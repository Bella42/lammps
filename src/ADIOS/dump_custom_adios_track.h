/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(custom/adios-track, DumpCustomADIOSTrack);
// clang-format on
#else

#ifndef LMP_DUMP_CUSTOM_ADIOS_TRACK_H
#define LMP_DUMP_CUSTOM_ADIOS_TRACK_H

#include "dump_custom.h"
#include "dump_custom_adios.h"


namespace LAMMPS_NS {

// class DumpCustomTrackADIOSInternal;

/**
*  Class using dai_track_feature to mitigate unusual variable usage
*/
class DumpCustomADIOSTrack : public DumpCustomADIOS {
 public:
  DumpCustomADIOSTrack(class LAMMPS *, int, char **);
  ~DumpCustomADIOSTrack() override;

 protected:
//   void openfile() override;
  void write() override;
//   void init_style() override;

 private:
//   DumpCustomADIOSInternal *internal;
};
}    // namespace LAMMPS_NS

#endif
#endif
