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
DumpStyle(custom/adios, DumpCustomADIOS);
// clang-format on
#else

#ifndef LMP_DUMP_CUSTOM_ADIOS_H
#define LMP_DUMP_CUSTOM_ADIOS_H

#include "dump_custom.h"
#include "adios2.h"
#include "adios_common.h"

namespace LAMMPS_NS {

// class DumpCustomADIOSInternal;

class DumpCustomADIOSInternal {

 public:
  DumpCustomADIOSInternal(){};
  ~DumpCustomADIOSInternal() = default;

  // name of adios group, referrable in adios2_config.xml
  const std::string ioName = "custom";
  adios2::ADIOS *ad = nullptr;    // adios object
  adios2::IO io;                  // adios group of variables and attributes in this dump
  adios2::Engine fh;              // adios file/stream handle object
  // one ADIOS output variable we need to change every step
  adios2::Variable<double> varAtoms;
  // list of column names for the atom table
  // (individual list of 'columns' string)
  std::vector<std::string> columnNames;
};


// class DumpCustomADIOS : virtual public DumpCustom {
class DumpCustomADIOS : public DumpCustom {
 public:
  DumpCustomADIOS(class LAMMPS *, int, char **);
  ~DumpCustomADIOS() override;

 protected:
  void openfile() override;
//   void write() override;
  void init_style() override;
  DumpCustomADIOSInternal *internal;

 private:
};
}    // namespace LAMMPS_NS

#endif
#endif
