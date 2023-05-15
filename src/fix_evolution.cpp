#include "fix_evolution.h"
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
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"

#define PI 3.14159265


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

Evolution2D::Evolution2D(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
  
{

  if (narg != 14) error->all(FLERR,"Illegal evolution2D command");
  eta = utils::numeric(FLERR,arg[3],false,lmp);
  if (eta < 0.0) error->all(FLERR,"Diffusion coefficient must be >= 0.0");
  tau_v = utils::numeric(FLERR,arg[4],false,lmp);
  tau_n = utils::numeric(FLERR,arg[5],false,lmp);
  epsilon = utils::numeric(FLERR,arg[6],false,lmp);
  alpha = utils::numeric(FLERR,arg[7],false,lmp);
  alphaq = utils::numeric(FLERR,arg[8],false,lmp);
  betaq = utils::numeric(FLERR,arg[9],false,lmp);
  comm_radius = utils::numeric(FLERR,arg[10],false,lmp);
  region = domain->get_region_by_id(arg[11]);
  idregion = utils::strdup(arg[11]);
  controller_flag = utils::numeric(FLERR,arg[12],false,lmp);
  seed = utils::numeric(FLERR,arg[13],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal evolution2D command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
}
/* ---------------------------------------------------------------------- */

Evolution2D::~Evolution2D()
{
  
  delete random;

}

/* ---------------------------------------------------------------------- */


int Evolution2D::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void Evolution2D::init()
{

}

float wca(float r){
  float a = 8;
  bool leaky = true;
  
  if(leaky) r += (pow(2, 1./a)-1);
  if(r > pow(2, 1./a)){
    return 0.0;
  }
  return 4 * (pow(r, -2*a) - pow(r, -a)) + 1;
}


void Evolution2D::initial_integrate(int vflag)
{
 
  // function to update x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **mu = atom->mu;
  double **phi = atom->phi;
  double *q_reward = atom->q_reward;
  double mag;

  double F_min = 0;
  double D_min = 0;


  // Using 'radius' to store particle orentiation angle
  //double *phi = atom->radius;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int step = update->ntimestep;
  
  double q_received = 0.0;
  double comm_sq = comm_radius * comm_radius;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dt = update->dt;
  sqrtdt = sqrt(dt);
  
  // Set initial particle orientation
  if (step <= 1) {
    for (int i = 0; i < nlocal; i++){       
       // Initialize Active vector with random direction at beginining of simulation
       mu[i][0] = random->gaussian();       
       mu[i][1] = random->gaussian(); 
       mu[i][2] = 0; // 2D

       // Normalize active vector
       mag = sqrt(pow(mu[i][0],2)+pow(mu[i][1],2));
       double mag_inv = 1.0 / mag;
       mu[i][0] *= mag_inv;
       mu[i][1] *= mag_inv;
    
       // Initialize parameters with random direction at beginining of simulation
       phi[i][0] = random->uniform();
       phi[i][1] = random->uniform();
       phi[i][2] = random->uniform();
       q_reward[i] = 0.0;
       }
  }
  int *ilist,*jlist,*numneigh,**firstneigh;
  int i, j, ii, jj, inum, jnum;
  double rsq, xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  NeighList *list = neighbor->lists[0];  // Works but might break
 
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  // Update parameters and reward
  if(step > 1){
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely;
      if(rsq == 0) continue;

      if (rsq < comm_sq) {
        if (q_reward[j]  >= q_reward[i]){
	        phi[i][0] += alpha* (phi[j][0] - phi[i][0]) * dt;
	        phi[i][1] += alpha* (phi[j][1] - phi[i][1]) * dt;
	        phi[i][2] += alpha* (phi[j][2] - phi[i][2]) * dt;
	        q_reward[i] += alpha* (q_reward[j] - q_reward[i]) * dt;
	      }
	      if (q_reward[i]  >= q_reward[j]){
	        phi[j][0] += alpha* (phi[i][0] - phi[j][0]) * dt;
	        phi[j][1] += alpha* (phi[i][1] - phi[j][1]) * dt;
	        phi[j][2] += alpha* (phi[i][2] - phi[j][2]) * dt;
		q_reward[j] += alpha* (q_reward[i] - q_reward[j]) * dt;
	      }
      }
      /*if (q_reward[j]  >= q_reward[i]){
	        phi[i][0] += alpha* sin(phi[j][0] - phi[i][0]) * dt;
	        phi[i][1] += alpha* sin(phi[j][1] - phi[i][1]) * dt;
	        phi[i][2] += alpha* sin(phi[j][2] - phi[i][2]) * dt;
	        q_reward[i] += alpha* (q_reward[j] - q_reward[i]) * dt;
	      }
	      if (q_reward[i]  >= q_reward[j]){
	        phi[j][0] += alpha* sin(phi[i][0] - phi[j][0]) * dt;
	        phi[j][1] += alpha* sin(phi[i][1] - phi[j][1]) * dt;
	        phi[j][2] += alpha* sin(phi[i][2] - phi[j][2]) * dt;
		q_reward[j] += alpha* (q_reward[i] - q_reward[j]) * dt;
      }*/
    }

  }
  }
  // Integrator 
  for (int i = 0; i < nlocal; i++)
	  
    if (mask[i] & groupbit) {
    
    
      
      // Mutation noise
      phi[i][0] += random->gaussian() * sqrt(2*dt*eta) + wca(phi[i][0] + 1) - wca(2 - phi[i][0]);
      phi[i][1] += random->gaussian() * sqrt(2*dt*eta) + wca(phi[i][1] + 1) - wca(2 - phi[i][1]);
      phi[i][2] += random->gaussian() * sqrt(2*dt*eta) + wca(phi[i][2] + 1) - wca(2 - phi[i][2]);
      
      
      /*phi[i][0] += random->gaussian() * sqrt(2*dt*eta);
      phi[i][1] += random->gaussian() * sqrt(2*dt*eta);
      phi[i][2] += random->gaussian() * sqrt(2*dt*eta);*/
      
      
      if(region->match(x[i][0], x[i][1], x[i][2])){
        q_received = 0.75;
      } else {
        q_received = 0.25;
      }
      // Update reward
      q_reward[i] += alphaq*(q_received - betaq * q_reward[i]) *dt;

      
      //double Fa = 0.5 * (1 + cos(phi[i][1] * 2 * 3.141592));
      double Fa = 0;
      switch(controller_flag){
        case 0:
        Fa = phi[i][0];
        break;
        
        case 1:
        Fa = 0.5*(1+cos(PI * (phi[i][0]+1)));
        break;
        
        case 2:
        Fa = 0.5*(1+cos(2*PI*phi[i][0]));
        break;
        
        case 3:
        Fa = -2*phi[i][0] + 1;
        if(phi[i][0] >= 0.5) Fa = 2*phi[i][0] - 1;
        break;
        
        case 4:
        Fa = 0.5*(1+cos(2*PI*phi[i][0] + PI));
        break;
        
        case 5:
        Fa = 2*phi[i][0];
        if(phi[i][0] >= 0.5) Fa = -2*phi[i][0] + 2;
        break;
        
        default:
        Fa = phi[i][0];
        break;
      }
      if( phi[i][1] >= 0.5) Fa = 2*phi[i][1] - 1;
      
      
      double D = 0.01;
    

      // Update velocities
      v[i][0] += dt*(Fa*mu[i][0] - v[i][0] + f[i][0])/tau_v;
      v[i][1] += dt*(Fa*mu[i][1] - v[i][1] + f[i][1])/tau_v;
      v[i][2] = 0;

      // overdamped
      /*v[i][0] = mu[i][0] + f[i][0];
      v[i][1] = mu[i][1] + f[i][1];
      v[i][2] = 0;*/

       // Update Active Vector
      double mux, muy;
      mux = epsilon * (mu[i][1]*mu[i][1]*v[i][0] - mu[i][0]*mu[i][1]*v[i][1]) / tau_n;
      muy = epsilon * (mu[i][0]*mu[i][0]*v[i][1] - mu[i][0]*mu[i][1]*v[i][0]) / tau_n;
      mu[i][0] += dt * mux;
      mu[i][1] += dt * muy; 
      
      // Initialise angular noise
      double ang_noise = random->gaussian();
      ang_noise *= sqrt(2*D*dt);
      
      mux = mu[i][0];
      muy = mu[i][1];

      mu[i][0] = cos(ang_noise) * mux - sin(ang_noise) * muy;
      mu[i][1] = sin(ang_noise) * mux + cos(ang_noise) * muy;

     /* double vx, vy;
      vx = v[i][0];
      vy = v[i][1];

      v[i][0] = cos(ang_noise) * vx - sin(ang_noise) * vy;
      v[i][1] = sin(ang_noise) * vx + cos(ang_noise) * vy;
	*/
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
