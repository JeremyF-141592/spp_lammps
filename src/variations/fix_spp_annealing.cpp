#include "fix_spp.h"
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

SelfAlignment2D::SelfAlignment2D(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
  
{

  if (narg != 8) error->all(FLERR,"Illegal selfAlignment2D command");
  D = utils::numeric(FLERR,arg[3],false,lmp);
  if (D < 0.0) error->all(FLERR,"Diffusion coefficient must be >= 0.0");
  tau_v = utils::numeric(FLERR,arg[4],false,lmp);
  tau_n = utils::numeric(FLERR,arg[5],false,lmp);
  epsilon = utils::numeric(FLERR,arg[6],false,lmp);
  seed = utils::numeric(FLERR,arg[7],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal selfAlignment2D command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
}
/* ---------------------------------------------------------------------- */

SelfAlignment2D::~SelfAlignment2D()
{

  delete random;

}

/* ---------------------------------------------------------------------- */


int SelfAlignment2D::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void SelfAlignment2D::init()
{

}


void SelfAlignment2D::initial_integrate(int vflag)
{
 
  // function to update x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **mu = atom->mu;
  double mag;


  // Using 'radius' to store particle orentiation angle
  //double *phi = atom->radius;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int step = update->ntimestep;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dt = update->dt;
  sqrtdt = sqrt(dt);

  // Set initial particle orientation
 if (step <= 1) {
    for (int i = 0; i < nlocal; i++){       
       // Initialise Active vector with random direction at beginining of simulation
       mu[i][0] = random->gaussian();       
       mu[i][1] = random->gaussian(); 
       mu[i][2] = 0; // 2D

       // Normailise active vector
       mag = sqrt(pow(mu[i][0],2)+pow(mu[i][1],2));
       double mag_inv = 1.0 / mag;
       mu[i][0] *= mag_inv;
       mu[i][1] *= mag_inv;
    }
  }
  double true_D = D;
  double fifth = (update->laststep - update->firststep) / 5.;
  if(step > 3*fifth){
    true_D *= (1. - (step - 3*fifth)/(2*fifth));
  } else if(step > fifth){
    true_D *= (step - fifth) / (2*fifth);
  } else {
    true_D = 0;
  }
  // Integrator 
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      
      // Initialise angular noise
      double ang_noise = random->gaussian();

      
      ang_noise *= sqrt(2*true_D*dt)/ tau_n;
      
       // Update Active Vector
      double mux, muy;

      mux = epsilon * (mu[i][1]*mu[i][1]*v[i][0] - mu[i][0]*mu[i][1]*v[i][1]) / tau_n;
      muy = epsilon * (mu[i][0]*mu[i][0]*v[i][1] - mu[i][0]*mu[i][1]*v[i][0]) / tau_n;
           
      // Update velocities
      v[i][0] += dt*(mu[i][0] - v[i][0] + f[i][0])/tau_v;
      v[i][1] += dt*(mu[i][1] - v[i][1] + f[i][1])/tau_v;
      v[i][2] = 0;
      
      // Apply angular noise
      double tmp_mu0 = mu[i][0];
      mu[i][0] += dt * mux - ang_noise * mu[i][1];
      mu[i][1] += dt * muy + ang_noise * tmp_mu0;


      // Normalise updated Active vector
      mag = sqrt(pow(mu[i][0],2)+pow(mu[i][1],2)+pow(mu[i][2],2));
      double mu_mag_inv = 1.0/mag;
      mu[i][0] *= mu_mag_inv;
      mu[i][1] *= mu_mag_inv;
      
      // Update positions
      x[i][0] +=  v[i][0]*dt;
      x[i][1] +=  v[i][1]*dt;
      x[i][2] = 0;                                                                  
      }
}

/* ---------------------------------------------------------------------- */
