/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(base/rna,PairBaseRNA);
// clang-format on
#else

#ifndef LMP_PAIR_BASE_RNA_H
#define LMP_PAIR_BASE_RNA_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBaseRNA : public Pair {
 public:
  PairBaseRNA(class LAMMPS *);
  ~PairBaseRNA() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void ev_tally6(int , int , int , int , int , int , double,
  double , double , double , double *,
  double *, double *, double *, double *, 
  double *, double *, double *, double *, 
  double *, double *, double *, double *, 
  double *, double *, double *, double *, 
  double *, double *, double *, 
  double *, double *, double *, 
  double *, double *, double *, 
  double *, double *, double *);
  
 protected:
  double cut_global;
  double **cut;
  double **ubp0, **kr, **r0, **k_theta, **k_phi;
  double **theta1, **theta2, **phi1, **phi2;
  double **offset;

  virtual void allocate();
};
}    // namespace LAMMPS_NS
#endif
#endif