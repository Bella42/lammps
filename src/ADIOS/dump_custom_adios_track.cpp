/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Norbert Podhorszki (ORNL) (DAI support: Kira Duwe (OVGU))
------------------------------------------------------------------------- */

#include "dump_custom_adios.h"
#include "dump_custom_adios_track.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

#include "adios2.h"
#include "adios_common.h"

//TODO: include JULEA header + DAI header

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

DumpCustomADIOSTrack::DumpCustomADIOSTrack(LAMMPS *lmp, int narg, char **arg) : DumpCustomADIOS(lmp, narg, arg)
{
  // create a default adios2_config.xml if it doesn't exist yet.
  FILE *cfgfp = fopen("adios2_config.xml", "r");
  if (!cfgfp) {
    cfgfp = fopen("adios2_config.xml", "w");
    if (cfgfp) fputs(default_config, cfgfp);
  }
  if (cfgfp) fclose(cfgfp);
}

/* ---------------------------------------------------------------------- */

DumpCustomADIOSTrack::~DumpCustomADIOSTrack()
{
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOSTrack::write()
{
  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // nme = # of dump lines this proc contributes to dump

  nme = count();

  // ntotal = total # of atoms in snapshot
  // atomOffset = sum of # of atoms up to this proc (exclusive prefix sum)

  bigint bnme = nme;
  MPI_Allreduce(&bnme, &ntotal, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  bigint atomOffset;    // sum of all atoms on processes 0..me-1
  MPI_Scan(&bnme, &atomOffset, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  atomOffset -= nme;    // exclusive prefix sum needed

  // Now we know the global size and the local subset size and offset
  // of the atoms table
  auto nAtomsGlobal = static_cast<size_t>(ntotal);
  auto startRow = static_cast<size_t>(atomOffset);
  auto nAtomsLocal = static_cast<size_t>(nme);
  auto nColumns = static_cast<size_t>(size_one);
  internal->varAtoms.SetShape({nAtomsGlobal, nColumns});
  internal->varAtoms.SetSelection({{startRow, 0}, {nAtomsLocal, nColumns}});

  // insure filewriter proc can receive everyone's info
  // limit nmax*size_one to int since used as arg in MPI_Rsend() below
  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (nme > maxbuf) {
    if ((bigint) nme * size_one > MAXSMALLINT) error->all(FLERR, "Too much per-proc info for dump");
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf, (maxbuf * size_one), "dump:buf");
  }
  if (sort_flag && sortcol == 0 && nme > maxids) {
    maxids = nme;
    memory->destroy(ids);
    memory->create(ids, maxids, "dump:ids");
  }

  if (sort_flag && sortcol == 0)
    pack(ids);
  else
    pack(nullptr);
  if (sort_flag) sort();

  //TODO: declare how ADIOS2 variable is used -> DAI Call which features in atoms.data should have separate columns in table
  //void j_dai_track_feature(gchar const* file_name, gchar const* var_name, gchar const* feature_name, size_t index_in_data_vector)

  openfile();

  internal->fh.BeginStep();
  // write info on data as scalars (by me==0)
  if (me == 0) {
    internal->fh.Put<uint64_t>("ntimestep", update->ntimestep);
    internal->fh.Put<int>("nprocs", nprocs);

    internal->fh.Put<double>("boxxlo", boxxlo);
    internal->fh.Put<double>("boxxhi", boxxhi);
    internal->fh.Put<double>("boxylo", boxylo);
    internal->fh.Put<double>("boxyhi", boxyhi);
    internal->fh.Put<double>("boxzlo", boxzlo);
    internal->fh.Put<double>("boxzhi", boxzhi);

    if (domain->triclinic) {
      internal->fh.Put<double>("boxxy", boxxy);
      internal->fh.Put<double>("boxxz", boxxz);
      internal->fh.Put<double>("boxyz", boxyz);
    }
  }
  // Everyone needs to write scalar variables that are used as dimensions and
  // offsets of arrays
  internal->fh.Put<uint64_t>("natoms", ntotal);
  internal->fh.Put<int>("ncolumns", size_one);
  internal->fh.Put<uint64_t>("nme", bnme);
  internal->fh.Put<uint64_t>("offset", atomOffset);
  // now write the atoms
  internal->fh.Put<double>("atoms", buf);
  internal->fh.EndStep();    // I/O will happen now...

  if (multifile) { internal->fh.Close(); }
}
