#ifdef FIX_CLASS

FixStyle(spprt, SelfAlignment2DRT)

#else
#ifndef LMP_selfAlignment2DRT_H
#define LMP_selfAlignment2DRT_H

#include "fix.h"

namespace LAMMPS_NS {

class SelfAlignment2DRT : public Fix {
 public:
  SelfAlignment2DRT(class LAMMPS *, int, char **);
  virtual ~SelfAlignment2DRT();
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
  double turn_rate;
  double r_time, t_time_max;
  double speed_modifier, effect_radius;
  char *id_temp;
  void compute_target();

};

}

#endif
#endif
