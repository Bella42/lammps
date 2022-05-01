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
   Contributing author: Norbert Podhorszki (ORNL)
------------------------------------------------------------------------- */

#include "dump_custom_adios.h"

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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpCustomADIOS::DumpCustomADIOS(LAMMPS *lmp, int narg, char **arg) : DumpCustom(lmp, narg, arg)
{
  // create a default adios2_config.xml if it doesn't exist yet.
  FILE *cfgfp = fopen("adios2_config.xml", "r");
  if (!cfgfp) {
    cfgfp = fopen("adios2_config.xml", "w");
    if (cfgfp) fputs(default_config, cfgfp);
  }
  if (cfgfp) fclose(cfgfp);

  internal = new DumpCustomADIOSInternal();
  try {
    internal->ad = new adios2::ADIOS("adios2_config.xml", world, adios2::DebugON);
  } catch (std::ios_base::failure &e) {
    error->all(FLERR, "ADIOS initialization failed with error: {}", e.what());
  }

  internal->columnNames.reserve(nfield);
  for (int i = 0; i < nfield; ++i) { internal->columnNames.push_back(earg[i]); }
}

/* ---------------------------------------------------------------------- */

DumpCustomADIOS::~DumpCustomADIOS()
{
  internal->columnNames.clear();
  if (internal->fh) { internal->fh.Close(); }
  delete internal->ad;
  delete internal;
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::openfile()
{
  if (multifile) {
    // if one file per timestep, replace '*' with current timestep
    auto filecurrent = utils::star_subst(filename, update->ntimestep, padflag);
    internal->fh = internal->io.Open(filecurrent, adios2::Mode::Write, world);
    if (!internal->fh) error->one(FLERR, "Cannot open dump file {}", filecurrent);
  } else {
    if (!singlefile_opened) {
      internal->fh = internal->io.Open(filename, adios2::Mode::Write, world);
      if (!internal->fh) error->one(FLERR, "Cannot open dump file {}", filename);
      singlefile_opened = 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::init_style()
{
  // assemble column string from defaults and user values

  delete[] columns;
  std::string combined;
  int icol = 0;
  for (auto item : utils::split_words(columns_default)) {
    if (combined.size()) combined += " ";
    if (keyword_user[icol].size())
      combined += keyword_user[icol];
    else
      combined += item;
    ++icol;
  }
  columns = utils::strdup(combined);

  // setup boundary string

  domain->boundary_string(boundstr);

  // remove % from filename since ADIOS always writes a global file with
  // data/metadata
  char *ptr = strchr(filename, '%');
  if (ptr) {
    while (*ptr) {
      ptr[0] = ptr[1];
      ++ptr;
    }
  }

  /* The next four loops are copied from dump_custom_mpiio, but nothing is
   * done with them.
   * It is unclear why we need them here.
   * For metadata, variable[] will be written out as an ADIOS attribute if
   * nvariable>0
   */
  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable
  for (int i = 0; i < ncompute; i++) {
    compute[i] = modify->get_compute_by_id(id_compute[i]);
    if (!compute[i])
      error->all(FLERR, "Could not find dump custom/adios compute ID {}", id_compute[i]);
  }

  for (int i = 0; i < nfix; i++) {
    fix[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fix[i]) error->all(FLERR, "Could not find dump custom/adios fix ID {}", id_fix[i]);
    if (nevery % fix[i]->peratom_freq)
      error->all(FLERR, "dump custom/adios and fix {} with ID {} not computed at compatible times",
                 fix[i]->style, id_fix[i]);
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) error->all(FLERR, "Could not find dump custom/adios variable name");
    variable[i] = ivariable;
  }

  // set index and check validity of region
  if (idregion && !domain->get_region_by_id(idregion))
    error->all(FLERR, "Region {} for dump custom/adios does not exist", idregion);

  /* Define the group of variables for the atom style here since it's a fixed
   * set */
  internal->io = internal->ad->DeclareIO(internal->ioName);
  if (!internal->io.InConfigFile()) {
    // if not defined by user, we can change the default settings
    // BPFile is the default writer
    internal->io.SetEngine("BPFile");
    int num_aggregators = multiproc;
    if (num_aggregators == 0) num_aggregators = 1;
    auto nstreams = std::to_string(num_aggregators);
    internal->io.SetParameters({{"substreams", nstreams}});
    if (me == 0)
      utils::logmesg(lmp, "ADIOS method for {} is n-to-m (aggregation with {} writers)\n", filename,
                     nstreams);
  }

  internal->io.DefineVariable<uint64_t>("ntimestep");
  internal->io.DefineVariable<uint64_t>("natoms");

  internal->io.DefineVariable<int>("nprocs");
  internal->io.DefineVariable<int>("ncolumns");

  internal->io.DefineVariable<double>("boxxlo");
  internal->io.DefineVariable<double>("boxxhi");
  internal->io.DefineVariable<double>("boxylo");
  internal->io.DefineVariable<double>("boxyhi");
  internal->io.DefineVariable<double>("boxzlo");
  internal->io.DefineVariable<double>("boxzhi");

  internal->io.DefineVariable<double>("boxxy");
  internal->io.DefineVariable<double>("boxxz");
  internal->io.DefineVariable<double>("boxyz");

  internal->io.DefineAttribute<int>("triclinic", domain->triclinic);

  int *boundaryptr = reinterpret_cast<int *>(domain->boundary);
  internal->io.DefineAttribute<int>("boundary", boundaryptr, 6);

  auto nColumns = static_cast<size_t>(size_one);
  internal->io.DefineAttribute<std::string>("columns", internal->columnNames.data(), nColumns);
  internal->io.DefineAttribute<std::string>("columnstr", columns);
  internal->io.DefineAttribute<std::string>("boundarystr", boundstr);
  internal->io.DefineAttribute<std::string>("LAMMPS/dump_style", "custom");
  internal->io.DefineAttribute<std::string>("LAMMPS/version", lmp->version);
  internal->io.DefineAttribute<std::string>("LAMMPS/num_ver", std::to_string(lmp->num_ver));

  internal->io.DefineVariable<uint64_t>("nme",
                                        {adios2::LocalValueDim});    // local dimension variable
  internal->io.DefineVariable<uint64_t>("offset",
                                        {adios2::LocalValueDim});    // local dimension variable

  // atom table size is not known at the moment
  // it will be correctly defined at the moment of write
  size_t UnknownSizeYet = 1;
  internal->varAtoms = internal->io.DefineVariable<double>(
      "atoms", {UnknownSizeYet, nColumns}, {UnknownSizeYet, 0}, {UnknownSizeYet, nColumns});
}
