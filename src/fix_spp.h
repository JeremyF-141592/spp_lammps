#ifdef FIX_CLASS

FixStyle(spp, SelfAlignment2D)

#else
#ifndef LMP_selfAlignment2D_H
#define LMP_selfAlignment2D_H

#include "fix.h"

namespace LAMMPS_NS {

class SelfAlignment2D : public Fix {
 public:
  SelfAlignment2D(class LAMMPS *, int, char **);
  virtual ~SelfAlignment2D();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);


 private: 
 double dt, sqrtdt;

 protected:
  class RanMars *random;
  int seed;
  double t_start,t_stop,t_period,t_target,tsqrt;
  double gamma1,gamma2,gamma3;
  double D,cosphi,sinphi;
  double tau_v, tau_n;
  char *id_temp;
  void compute_target();

};

}

#endif
#endif
