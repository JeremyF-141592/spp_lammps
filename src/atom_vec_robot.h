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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(robot,AtomVecRobot);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_ROBOT_H
#define LMP_ATOM_VEC_ROBOT_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecRobot : public AtomVec {
 public:
  AtomVecRobot(class LAMMPS *);

  void grow_pointers() override;
  void data_atom_post(int) override;

 private:
  double **mu;
  double *rt_state;
  int *lr_state;
};

}    // namespace LAMMPS_NS

#endif
#endif
