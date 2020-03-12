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

typedef struct singleParticleDerivatives {
  double* divergence;
  double* gradient;
  double* laplacian;
} singleParticleDerivatives;

typedef struct mySingleParticle {
  double* coordinates; // x-y
  double* values; // density or velocity (vector!) or pressure or temperature or anything else
  int size_values;
  double mass;
  double density;
  neighborhood_options* particle_neighbours; // all the neighbouring particles infos related to the particle of interest
  singleParticleDerivatives* particle_derivatives; // divergence, gradient and laplacian (depending if the quantity carried by the particle is a scalar or vector)
} mySingleParticle;

mySingleParticle* create_array_of_particles(int nbParticles, int size_values, neighborhood* nh);

singleParticleDerivatives* initialize_particle_derivatives(int size_values);

void init1DSegmentWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nb_particles, int index_part, int size_values);
// void init1DSegmentWithParticles(GLfloat(* data)[14], double* x_lim, int nb_particles, double (*myFun)(double));

double myFunctionToDerive(double* x);

void computeDerivativesOfParticleQuantity(mySingleParticle* myPart, int index_part);

void computeDerivatiesAllParticles(mySingleParticle* myPart, int nbParticles);

#endif