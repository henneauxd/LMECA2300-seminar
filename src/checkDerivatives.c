
#include "checkDerivatives.h"



void init1DSegmentWithParticles(GLfloat(* data)[14], double* x_lim, int nb_particles, double (*myFun)(double*)) 
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
  double x_coord;
  for (int i = 0; i < nb_particles; i++) {
      x_coord = x_min + i*delta_x; // uniformly distributed points on the segment
      // coordinates
      data[i][0] = x_coord ; // X-dim
      data[i][1] = 0.0; // 1-D segment so don't care about the second dimension
      // function values
      data[i][8] = (myFun)(&x_coord);
      data[i][9] = 0.0; // 1-D segment so don't care about the second dimension
  }
  
}

double myFunctionToDerive(double* x) 
{
   return x[0]*x[0]; 
}



void computeDerivativesOfParticleQuantity(mySingleParticle* myPart, double kh, double* divergence, double* gradient) {
   double DENSITY = 0.0; // WARNING
   double MASS = 0.0; // WARNING
    for (int i = 0; i < NPTS; i++) {
        // Let's imagine we want only to compute the derivatives of the velocity field for the moment
        double val_node_x = data[i][2]; // velocity component in the x-direction
        double val_node_y = data[i][3]; // velocity component in the y-direction
        int nNeigh = nh[i].nNeighbours;
        double val_div = 0;
        double val_grad_x = 0;
        double val_grad_y = 0;
        double val_lapl = 0;
        double dens2 = pow(DENSITY, 2);
        neighbours* List = nh[i].list;
        if (nNeigh > 0) {
            for (int j = 0; j < nNeigh; j++) {
                int index_node2 = List->index;
                double distance = List->distance;
                double d_x = data[index_node2][0] - data[i][0];
                double d_y = data[index_node2][1] - data[i][1];
                
                /*
                 You can choose here the desired kernel function for your code.
                 */
                
                //double weight_x = grad_w_cubic(distance, kh, d_x);
                //double weight_y = grad_w_cubic(distance, kh, d_y);
                
                double weight_x = grad_w_lucy(distance, kh, d_x);
                double weight_y = grad_w_lucy(distance, kh, d_y);
                
                //double weight_x = grad_w_newquartic(distance, kh, d_x);
                //double weight_y = grad_w_newquartic(distance, kh, d_y);
                
                //double weight_x = grad_w_quinticspline(distance, kh, d_x);
                //double weight_y = grad_w_quinticspline(distance, kh, d_y);
                
                val_div += -MASS / DENSITY * ((data[index_node2][2] - val_node_x) * weight_x + (data[index_node2][3] - val_node_y) * weight_y);
                val_grad_x += -DENSITY * MASS * ((val_node_x / dens2) + (data[index_node2][2] / dens2)) * weight_x;
                val_grad_y += -DENSITY * MASS * ((val_node_y / dens2) + (data[index_node2][3] / dens2)) * weight_y;
                val_lapl += 2.0 * MASS / DENSITY * (val_node_x - data[index_node2][8]) * (d_x * weight_x + d_y * weight_y) / pow(distance,2);
                
                List = List->next;
            }
        }
        // All the values of the divergent gradient and laplacien are stored in the data table
        data[i][10] = val_div;
        data[i][11] = val_grad_x;
        data[i][12] = val_grad_y;
        data[i][13] = val_lapl;
    }
    
    //Computation of the error based on the already know function.
    for (int j = 0; j < NPTS; j++) {
        double exact = 3 * pow(data[j][0], 2);
        double error = exact - data[j][10];
    }
}