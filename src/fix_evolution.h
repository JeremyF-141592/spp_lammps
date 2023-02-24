#ifdef FIX_CLASS

FixStyle(evo2D, Evolution2D)

#else
#ifndef LMP_evolution2D_H
#define LMP_evolution2D_H

#include "fix.h"

namespace LAMMPS_NS {

class Evolution2D : public Fix {
 public:
  Evolution2D(class LAMMPS *, int, char **);
  virtual ~Evolution2D();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);


 private: 
 double dt, sqrtdt;

 protected:
  class RanMars *random;
  int seed;
  double eta;
  double tau_v, tau_n;
  int epsilon;
  double alpha, alphaq, betaq;
  double comm_radius;
  char *id_temp;

  char* idregion;
  class Region *region;

};

}

#endif
#endif