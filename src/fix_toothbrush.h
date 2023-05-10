#ifdef FIX_CLASS

FixStyle(toothbrush, Toothbrush2D)

#else
#ifndef LMP_toothbrush2D_H
#define LMP_toothbrush2D_H

#include "fix.h"

namespace LAMMPS_NS {

class Toothbrush2D : public Fix {
 public:
  Toothbrush2D(class LAMMPS *, int, char **);
  virtual ~Toothbrush2D();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);


 private: 
 double dt, sqrtdt;

 protected:
  class RanMars *random;
  int seed;
  double D;
  int epsilon;
  double tau_v, tau_n;
  double J, bias_std;
  char *id_temp;

};

}

#endif
#endif
