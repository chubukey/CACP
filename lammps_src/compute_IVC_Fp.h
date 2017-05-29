/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(IVCFp,ComputeIVCFp)

#else

#ifndef LMP_COMPUTE_IVC_FP_H
#define LMP_COMPUTE_IVC_FP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeIVCFp : public Compute {
 public:
  ComputeIVCFp(class LAMMPS *, int, char **);
  ~ComputeIVCFp();
  double compute_scalar();
  void compute_vector();
  void init();
  double memory_usage();

 protected:
  int me;
  int nvalues;
  char **ids;
  int maxatom;

 private:
  int nnode_inter;
  int invoked_scalar;
  int nsystem,nplane,nFp;
  void decompose_into_plane(double * , double **, double *,double,int);
  void plane_summation(double **, int, double *);
  void Frotate(double *, double *, double *);
  double *compute_one();
  double **Fp_tot_plane;
  double *Fp_tot;
  double *natoms_plane_one,*natoms_plane_all;
  double natoms_allslip_one,natoms_allslip_all;
  int whichplane;             // index of slip plane atom belong to
  int flag_Fp_yes;            // if this atom has non-trivial Fp
  char *compute_id;
  double onecount;
  double allcount;
};

}

#endif
#endif

