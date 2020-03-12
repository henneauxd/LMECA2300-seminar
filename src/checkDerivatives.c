
#include "checkDerivatives.h"

singleParticleDerivatives* initialize_particle_derivatives(int size_values)
{
  singleParticleDerivatives* my_deriv = malloc(sizeof(singleParticleDerivatives));
  
  if (size_values == 1) { // scalar quantity
      my_deriv->divergence = NULL; //  divergence of a scalar has no meaning (= gradient)
      my_deriv->gradient = (double*) calloc(2,sizeof(double));
      my_deriv->laplacian = (double*) calloc(1, sizeof(double));
  }
  else if (size_values == 2) { // vector quantity
      my_deriv->divergence = (double*) calloc(1, sizeof(double));
      my_deriv->gradient = (double*) calloc(4,sizeof(double));
      my_deriv->laplacian = (double*) calloc(2, sizeof(double));
  }
  else {
    printf("---------------  The size of the quantity for which you try to compute the derivatives is not correct ------------------------");
    return 0;
  }
  
  return my_deriv;
}

mySingleParticle* create_array_of_particles(int nbParticles, int size_values, neighborhood* nh) 
{
  double x_lim[2] = {-1.0, 1.0};
  mySingleParticle* my_particles = malloc(nbParticles*sizeof(mySingleParticle));
  for (int i = 0; i < nbParticles; i++) {
    my_particles[i].coordinates = calloc(2,sizeof(double));
    my_particles[i].values = calloc(size_values,sizeof(double));
    my_particles[i].mass = 0.0;
    my_particles[i].density = 0.0;
    my_particles[i].size_values = size_values;
    my_particles[i].particle_neighbours = neighborhood_options_init(0.0, 0.0);
    my_particles[i].particle_neighbours->nh = &(nh[i]);
    my_particles[i].particle_derivatives = initialize_particle_derivatives(size_values);
    
    init1DSegmentWithParticles(x_lim, my_particles[i].coordinates, my_particles[i].values, &(my_particles[i].mass), &(my_particles[i].density), nbParticles, i, size_values);  
  }
  return my_particles;
}

void initSquareWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int* nb_particles, int index_part, int size_values)// double (*myFun)(double*)) 
{
  double x_min = x_lim[0];
  double x_max = x_lim[1];
  double y_min = x_lim[2];
  double y_max = x_lim[3];
  if (x_min > x_max) {
      double x_temp = x_min;
      x_min = x_max;
      x_max = x_temp;
  }
  if (y_min > y_max) {
      double y_temp = y_min;
      y_min = y_max;
      y_max = y_temp;
  }
  double length_x = fabs(x_max-x_min);
  double delta_x = length_x / ((double)nb_particles[0] - 1.0);
  double length_y = fabs(y_max-y_min);
  double delta_y = length_y / ((double)nb_particles[1] - 1.0);
  for (int i=0; i<nb_particles[1]; i++) {
      for (int j=0; j<nb_particles[0]; j++) {
	x_coord
	
      }
  }
  
  double x_coord = x_min + index_part*delta_x; // uniformly distributed points on the segment
  // coordinates
  coord[0] = x_coord ; // X-dim
  coord[1] = 0.0; // 1-D segment so don't care about the second dimension
  // function values
  values[0] = myFunctionToDerive(coord);//(myFun)(&x_coord);
  if (size_values > 1) values[0] = 0.0; // 1-D segment so don't care about the second dimension
  *mass = 1.0;
  *density = (double)nb_particles / length;//1.0; 
}

void init1DSegmentWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nb_particles, int index_part, int size_values)// double (*myFun)(double*)) 
{
  double x_min = x_lim[0];
  double x_max = x_lim[1];
  if (x_min > x_max) {
      double x_temp = x_min;
      x_min = x_max;
      x_max = x_temp;
  }
  double length = fabs(x_max-x_min);
  double delta_x = length / ((double)nb_particles - 1.0);
  double x_coord = x_min + index_part*delta_x; // uniformly distributed points on the segment
  // coordinates
  coord[0] = x_coord ; // X-dim
  coord[1] = 0.0; // 1-D segment so don't care about the second dimension
  // function values
  values[0] = myFunctionToDerive(coord);//(myFun)(&x_coord);
  if (size_values > 1) values[0] = 0.0; // 1-D segment so don't care about the second dimension
  *mass = 1.0;
  *density = (double)nb_particles / length;//1.0; 
}

double myFunctionToDerive(double* x) 
{
   return x[0]*x[0]; 
}


// void init1DSegmentWithParticles(GLfloat(* data)[14], double* x_lim, int nb_particles, double (*myFun)(double*)) 
// {
//   double x_min = x_lim[0];
//   double x_max = x_lim[1];
//   if (x_min > x_max) {
//       double x_temp = x_min;
//       x_min = x_max;
//       x_max = x_temp;
//   }
//   double length = fabs(x_max-x_min);
//   double delta_x = length / ((double)nb_particles - 1.0);
//   double x_coord;
//   for (int i = 0; i < nb_particles; i++) {
//       x_coord = x_min + i*delta_x; // uniformly distributed points on the segment
//       // coordinates
//       data[i][0] = x_coord ; // X-dim
//       data[i][1] = 0.0; // 1-D segment so don't care about the second dimension
//       // function values
//       data[i][8] = (myFun)(&x_coord);
//       data[i][9] = 0.0; // 1-D segment so don't care about the second dimension
//   }
//   
// }


void computeDerivativesOfParticleQuantity(mySingleParticle* myPart, int index_part) {
   
  mySingleParticle* local_part = &(myPart[index_part]);
  double * divergence_part = local_part->particle_derivatives->divergence;
  double * grad_part = local_part->particle_derivatives->gradient;
  double * laplacian_part = local_part->particle_derivatives->laplacian;

//   double * divergence_part;
//   double * grad_part;
//   double * laplacian_part;
//   
//   if (myPart->size_values == 1) { // scalar quantity
//       divergence_part = NULL; //  divergence of a scalar has no meaning (= gradient)
//       grad_part = (double*) calloc(2,sizeof(double));
//       laplacian_part = (double*) calloc(1, sizeof(double));
//   }
//   else if (myPart->size_values == 2) { // vector quantity
//       divergence_part = (double*) calloc(1, sizeof(double));
//       grad_part = (double*) calloc(4,sizeof(double));
//       laplacian_part = (double*) calloc(2, sizeof(double));
//   }
//   else {
//     printf("---------------  The size of the quantity for which you try to compute the derivatives is not correct ------------------------");
//     return 0;
//   }

  int nNeigh = local_part->particle_neighbours->nh->nNeighbours;
  neighbours* List = local_part->particle_neighbours->nh->list;
  double kh = local_part->particle_neighbours->kh;
  for (int j = 0; j < nNeigh; j++) {
      int index_j = List->index;
      double d_ij = List->distance;
      double d_x_ij =  myPart[index_j].coordinates[0] - local_part->coordinates[0];//  data[index_node2][0] - data[i][0];
//       double d_ij = fabs(d_x_ij);
      double d_y_ij = myPart[index_j].coordinates[1] - local_part->coordinates[1];//data[index_node2][1] - data[i][1];
      
      /*
	You can choose here the desired kernel function for your code.
	*/
      
//       double weight_x = grad_w_cubic(d_ij, kh, d_x_ij);
//       double weight_y = grad_w_cubic(d_ij, kh, d_y_ij);
//       if (index_part == 1) printf("value = %d\n",index_j);
      double weight_x = grad_w_lucy(d_ij, kh, d_x_ij);
      double weight_y = grad_w_lucy(d_ij, kh, d_y_ij);
//       if (index_part == 1) printf("value = %2.15f \n",kh);
      
//       double weight_x = grad_w_newquartic(d_ij, kh, d_x_ij);
//       double weight_y = grad_w_newquartic(d_ij, kh, d_y_ij);
      
      //double weight_x = grad_w_quinticspline(distance, kh, d_x);
      //double weight_y = grad_w_quinticspline(distance, kh, d_y);
      
      // --- Divergence ---
      if (myPart->size_values == 2) divergence_part[0] += (1.0 / local_part->density) * ( (myPart[index_j].values[0] - local_part->values[0]) * weight_x 
							  + (myPart[index_j].values[1] - local_part->values[1]) * weight_y ) * myPart[index_j].mass;
      // --- Gradient & Laplacian ---
      for (int k=0; k < myPart->size_values; k++) {
	// d/dx value_k
	grad_part[2*k] += -local_part->density * myPart[index_j].mass * ((local_part->values[k] / (local_part->density*local_part->density)) 
								      + (myPart[index_j].values[k] / (myPart[index_j].density*myPart[index_j].density))) * weight_x;
// 	if (index_part == 1) printf("value = %2.15f \n",grad_part[2*k]);							      
	// d/dy value_k
	grad_part[2*k+1] += -local_part->density * myPart[index_j].mass * ((local_part->values[k] / (local_part->density*local_part->density)) 
								      + (myPart[index_j].values[k] / (myPart[index_j].density*myPart[index_j].density))) * weight_y;
	// d^2/dx^2 value_k + d^2/dy^2 value_k
	laplacian_part[k] += 2.0 * (myPart[index_j].mass / myPart[index_j].density) * (local_part->values[k] - myPart[index_j].values[k]) * (d_x_ij * weight_x + d_y_ij * weight_y) /(d_ij*d_ij);
      }
      
      List = List->next;
  }

//   free(divergence_part);
//   free(grad_part);
//   free(laplacian_part);

}


void computeDerivatiesAllParticles(mySingleParticle* myPart, int nbParticles) {
  
  printf("i  (X,Y)	f(x)			Grad_x		Grad_y		Laplacian\n");
  for (int i=0; i<nbParticles; i++) {
//   int i = 75;
      computeDerivativesOfParticleQuantity(myPart, i);
      double* x = myPart[i].coordinates;
      double* grad = myPart[i].particle_derivatives->gradient;
      double* lapl = myPart[i].particle_derivatives->laplacian;
      double value = myPart[i].values[0];
      printf("%d  (%2.3f,%2.3f)		%2.6f        %2.6f		%2.6f		%2.6f\n",i,x[0],x[1],value,grad[0],grad[1],lapl[0]);
  }
  
  
  
}  
  
  