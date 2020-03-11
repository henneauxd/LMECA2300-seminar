#ifndef __CHECKDERIVATIES__
#define __CHECKDERIVATIES__

#include "BOV.h"
#include "kernel.h"
#include "neighborhood_search.h"
#include <time.h>
#include <math.h>

enum fieldNames {
  Density = 1,
  Velocity = 2,
  Pressure = 3,
  Temperature = 4
};

typedef struct mySingleParticle {
  double* coordinates;
  double* values;
  int size_values;
  neighborhood* particle_neighbours;
} mySingleParticle;

typedef struct singleParticleDerivatives {
  double* divergence;
  double* gradient;
  double* laplacian;
} singleParticleDerivatives;


void init1DSegmentWithParticles(GLfloat(* data)[14], double* x_lim, int nb_particles, double (*myFun)(double));

double myFunctionToDerive(double* x);


#endif