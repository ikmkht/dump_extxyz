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
DumpStyle(extxyz,DumpEXTXYZ);
// clang-format on
#else

#ifndef LMP_DUMP_EXTXYZ_H
#define LMP_DUMP_EXTXYZ_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpEXTXYZ : public Dump {
 public:
  DumpEXTXYZ(class LAMMPS *, int, char **);
  ~DumpEXTXYZ() override;

 protected:
  int ntypes;
  char **typenames;

  void init_style() override;
  void write_header(bigint) override;
  typedef void (DumpEXTXYZ::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;
  void header_binary(bigint);
  void header_binary_triclinic(bigint);
  typedef void (DumpEXTXYZ::*FnPtrPack)(tagint *);
  FnPtrPack pack_choice;    // ptr to pack functions
  void pack(tagint *);
  void pack_triclinic(tagint *);
  int convert_string(int, double *) override;
  void write_data(int, double *) override;
  int modify_param(int, char **) override;

  typedef void (DumpEXTXYZ::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;    // ptr to write data functions
  void write_string(int, double *);
  void write_lines(int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
