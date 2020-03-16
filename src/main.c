#include "neighborhood_search.h"
#include "kernel.h"
#include "checkDerivatives.h"
#include <math.h>

static int NPTS_BOUNDARIES = 10;
static int NPTS_DOMAIN = 10*10;
int NPTS = 10*10 + 4*10;//NPTS_DOMAIN + 4*NPTS_BOUNDARIES;

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void colormap(float v, float color[3])
{
	float v1 = 3.5 * (v - 0.7);
	float v2 = 1.25 * v;
	float v3 = fminf(0.5, v) * 2.0;

	color[0] = -v1 * v1 + 1.0f;
	color[1] = 6.0f * v2 * v2 * (1.0f - v2);
	color[2] = 5.5f * v3 * (1.0f - v3) * (1.0f - v3);

	// alternative: classical jet colormap
	// color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	// color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	// color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

static void myColormap(float v, float color[3], float v_max, float v_min)
{
	v -= v_min;
	v /= (v_max - v_min);

	// alternative: classical jet colormap
	color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

// function to fill the data table of the nPoints particles positions, speeds, colors and transparency and the coord table with the nPoints particles positions used to draw;
// data[i][0] == coord[i][0] && data[i][1] == coord[i][1]
void fillData(GLfloat(* data)[8])
{
// 	float rmax = 100.0 * sqrtf(2.0f);
  float rmax = 1.0 * sqrtf(2.0f);
// 	for (int i = 0; i < NPTS; i++) {
// 		data[i][0] = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
// 		printf("data[i][0] = %2.6f \n",data[i][0]);
// 		data[i][1] = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
// 		double r = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
// 		data[i][2] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
// 		data[i][3] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
// 		colormap(r / rmax, &data[i][4]); // fill color
// 		data[i][7] = 0.8f; // transparency
// 	}
// 	double x_lim[2] = {-1.0, 1.0};
// 	for (int i = 0; i < NPTS; i++) {
// // 	  double coord[2] = {data[i][0], data[i][1]};
// 	  double* coord = calloc(2,sizeof(double));
// // 	  double values[2] = {0.0,0.0};
// 	  double* values = calloc(2,sizeof(double));
// 	  double mass = 1.0;
// 	  double density = 1.0;
// 	  
// 	  init1DSegmentWithParticles(x_lim, coord, values, &mass, &density, NPTS, i, 1);
// 	  data[i][0] = coord[0];
// 	  data[i][1] = coord[1];
// 	  data[i][2] = 0.0; 
// 	  data[i][3] = 0.0;
// 	  double r = sqrt(values[0] * values[0]+ values[1] * values[1]);
// 	    myColormap(r, &data[i][4], 1.0); // fill color
// 	  data[i][7] = 0.8f; // transparency
// // 	  printf("data[i][0] = %2.6f \n",data[i][0]);
// 	}
	
	double x_lim[4] = {-1.0, 1.0, -1.0, 1.0};
	int i=0;
	int nbPart_x = sqrt(NPTS);
	int nbPart_y = nbPart_x;
	for (int ind_y = 0; ind_y < nbPart_y; ind_y++) {
	  for (int ind_x = 0; ind_x < nbPart_x; ind_x++) {
// 	  double coord[2] = {data[i][0], data[i][1]};
	  double* coord = calloc(2,sizeof(double));
// 	  double values[2] = {0.0,0.0};
	  double* values = calloc(2,sizeof(double));
	    double mass = 1.0;
	    double density = 1.0;
	    
// 	    init1DSegmentWithParticles(x_lim, coord, values, &mass, &density, NPTS, i, 1);
	    initSquareWithParticles(x_lim, coord, values, &mass, &density, nbPart_x, nbPart_y, ind_x, ind_y, 1);
	    data[i][0] = coord[0];
	    data[i][1] = coord[1];
	    
	    data[i][2] = 0.0; 
	    data[i][3] = 0.0;
	    double r = sqrt(values[0] * values[0]+ values[1] * values[1]);
	    myColormap(r, &data[i][4], 2.0, -2.0); // fill color
	    data[i][7] = 0.8f; // transparency
	    
	    i++;
// 	    printf("(values[0], values[1]) = (%2.6f,%2.6f) \n", values[0], values[1]);
	}
      }
}

// void myFillData(GLfloat(* data)[8], mySingleParticle* my_array_of_particles)
// {
// 	float rmax = 1.0 * sqrtf(2.0f);
// 	
// 	double x_lim[4] = {-1.0, 1.0, -1.0, 1.0};
// 	int i=0;
// 	int nbPart_x = sqrt(NPTS);
// 	int nbPart_y = nbPart_x;
// 	double max_grad = 2.0;
// 	double min_grad = -2.0;
// // 	for (int ind_y = 0; ind_y < nbPart_y; ind_y++) {
// // 	  for (int ind_x = 0; ind_x < nbPart_x; ind_x++) {
// // 	    if (my_array_of_particles[i].particle_derivatives->gradient[0] > max_grad) max_grad = my_array_of_particles[i].particle_derivatives->gradient[0];
// // 	  }
// // 	}
// 	for (int ind_y = 0; ind_y < nbPart_y; ind_y++) {
// 	  for (int ind_x = 0; ind_x < nbPart_x; ind_x++) {
// 	  double* coord = calloc(2,sizeof(double));
// 	  double* values = calloc(2,sizeof(double));
// 	  double mass = 1.0;
// 	  double density = 1.0;
// 	  
// // 	    init1DSegmentWithParticles(x_lim, coord, values, &mass, &density, NPTS, i, 1);
// 	  initSquareWithParticles(x_lim, coord, values, &mass, &density, nbPart_x, nbPart_y, ind_x, ind_y, 1);
// 	  data[i][0] = coord[0];
// 	  data[i][1] = coord[1];
// 	  
// 	  values[0] = my_array_of_particles[i].particle_derivatives->gradient[0];
// 	  double r = values[0];//sqrt(values[0] * values[0]+ values[1] * values[1]);
// 	  myColormap(r, &data[i][4], max_grad, min_grad); // fill color
// 	  data[i][7] = 0.8f; // transparency
// 	  
// 	  i++;
// 	}
//       }
// }

void myFillData(GLfloat(* data)[8], allParticles* my_array_of_particles, double* extrema)
{
	int nbPart = my_array_of_particles->nb_particles;
	for (int i=0; i<nbPart; i++) {
	    data[i][0] = my_array_of_particles->array_of_particles[i].coordinates[0];
	    data[i][1] = my_array_of_particles->array_of_particles[i].coordinates[1];
	    double r = my_array_of_particles->array_of_particles[i].values[0];
	    myColormap(r, &data[i][4], extrema[1], extrema[0]);
	    data[i][7] = 0.8f; // transparency
	}
}

void draw_particles(double* x_lim, GLfloat(* data)[8]) {
      bov_window_t* window = bov_window_new(1024, 780, "ANM Project: SPH");
      bov_window_set_color(window, (GLfloat[]){0.9f, 0.85f, 0.8f, 0.0f});

      /* send data to GPU, and receive reference to those data in a points object */
      bov_points_t *particles = bov_particles_new(data, NPTS, GL_STATIC_DRAW);

      /* setting particles appearance */
      bov_points_set_width(particles, 0.01);
      bov_points_set_outline_width(particles, 0.0025);
//       bov_points_set_outline_color(particles, ());


      /* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
      bov_points_scale(particles, (GLfloat[2]) {1.0, 1.0});
      bov_points_set_pos(particles, (GLfloat[2]) {0.0, -0.25});

      /* we got 0.2 at the top to write something. The screen goes from -1 to 1 */
      bov_text_t* msg =  bov_text_new(
	      (GLubyte[]) {"df/dx"},
	      GL_STATIC_DRAW);
      bov_text_set_pos(msg, (GLfloat[2]){-0.95, 0.82});
      bov_text_set_fontsize(msg, 0.1);

      while(!bov_window_should_close(window)){
// 	      fillData(data);
// 	      bov_particles_update(particles, data, NPTS);
	      
	      bov_particles_draw(window, particles, 0, BOV_TILL_END);
	      // bov_points_draw(window, particles, 0, BOV_TILL_END);
	      // bov_lines_draw(window, particles, 0, BOV_TILL_END);
	      // bov_triangles_draw(window, particles, 0, BOV_TILL_END);

	      bov_text_draw(window, msg);

	      // In your actual project, don't wait for events => bov_window_update(window)
	      bov_window_update_and_wait_events(window);
// 	      bov_window_update(window);
      }
      

      bov_text_delete(msg);
      bov_points_delete(particles);
      bov_window_delete(window);
}

void create_window_animation(GLfloat(* data)[8], bov_window_t* window, bov_points_t *particles) {
      window = bov_window_new(1024, 780, "ANM Project: SPH");
      bov_window_set_color(window, (GLfloat[]){0.9f, 0.85f, 0.8f, 0.0f});

      /* send data to GPU, and receive reference to those data in a points object */
      particles = bov_particles_new(data, NPTS, GL_STATIC_DRAW);

      /* setting particles appearance */
      bov_points_set_width(particles, 0.01);
      bov_points_set_outline_width(particles, 0.0025);
//       bov_points_set_outline_color(particles, ());


      /* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
      bov_points_scale(particles, (GLfloat[2]) {1.0, 1.0});
      bov_points_set_pos(particles, (GLfloat[2]) {0.0, -0.25});

      /* we got 0.2 at the top to write something. The screen goes from -1 to 1 */
      bov_text_t* msg =  bov_text_new(
	      (GLubyte[]) {"df/dx"},
	      GL_STATIC_DRAW);
      bov_text_set_pos(msg, (GLfloat[2]){-0.95, 0.82});
      bov_text_set_fontsize(msg, 0.1);
}

void display_particles(GLfloat(* data)[8], bov_window_t* window, bov_points_t *particles) {
  float transition_time = 1.0; 
  particles = bov_particles_update(particles,data,NPTS);
  bov_window_t* new_window = window;
  double tbegin = bov_window_get_time(window);
  while (bov_window_get_time(window) - tbegin < transition_time) {
    bov_particles_draw(window, particles, 0, BOV_TILL_END);
    bov_window_update(window);
  }
}

int main()
{
// 	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
// 	CHECK_MALLOC(data);
// 	// Seed the random
// 	time_t seed = time(NULL);
// 	//printf(" %u \n", seed);
// 	srand(seed);
// 	
// 	int nbParticles[2] = {(int) sqrt(NPTS), (int) sqrt(NPTS)};
// 	double domain_lim[4] = {-1.0,1.0,-1.0,1.0};
// 	allParticles* my_first_array_of_particles =  create_all_particles(nbParticles, domain_lim, 1, NULL, IN_DOMAIN, NULL, initFunction);
// 	
// 	// Creation of neighborhoods
// 	double timestep = 0.0;
// 	double maxspeed = 0.0;
// 	neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
// 	neighborhood* nh = options->nh;
// // 	neighborhood_update(options, nh, data, 0);
// 	neighborhood_update_new(options, nh, my_first_array_of_particles, 0);
// 	
// // 	// Creation of set of particles, for which to each of them are associated their coordinates, values of quantity they carry and their neighbourhood
// // 	mySingleParticle* my_array_of_particles = create_array_of_particles(NPTS, 1, nh);
// // 	// Compute the derivatives of the quantity carried by the particles 
// // 	computeDerivatiesAllParticles(my_array_of_particles, NPTS);
// 	
// 	
// 	allParticles* my_array_of_particles =  create_all_particles(nbParticles, domain_lim, 1, nh, IN_DOMAIN, NULL, initFunction);
// 	computeDerivativesAllParticles(my_array_of_particles, CUBIC);
// 	
// 	
// 	// Draw particles with color related to the derivatives of the quantity they carry
// 	myFillData(data, my_array_of_particles->array_of_particles);
// 	draw_particles(domain_lim, data);
// 
// 	neighborhood_options_delete(options,nh);
// 
// // 	double timestep = 0.5;
// // 	double maxspeed = 1;
// // 	neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
// // 	neighborhood* nh = options->nh;
// // 	int number_of_iterations = 10;
// // 	for (int iterations = 0; iterations < number_of_iterations;iterations++) {
// // 		if(iterations)
// // 			bouncyrandomupdate(data, timestep, options->half_length, maxspeed);
// // 		neighborhood_update(options, nh, data, iterations);
// // 		//kernel(data, nh, kh);
// // 	}
// // 	neighborhood_options_delete(options,nh);
// 
// 	free(data);
// 	delete_all_particles(my_array_of_particles);
// 	return EXIT_SUCCESS;
  
  
  // ***********************************************************************************************
  // ************************** HEAT EQUATION IN A SQUARE WITH DIRICHLET B.C. **********************
  // ***********************************************************************************************
  
    // *** 1 *** ASSIGN POSITIONS, DENSITY, AND MASS TO EVERY PARTICLES (neighbourhoods not yet defined) + BOUNDARY CONDITIONS
  
    // Creation of particles in the domain
    int nbParticles_domain[2] = {(int) sqrt(NPTS_DOMAIN), (int) sqrt(NPTS_DOMAIN)}; // number of particles in each direction (WARNING: same number of particles should be chosen for the moment)
    double domain_lim[4] = {0.0,1.0,0.0,1.0}; // limits of the computational domain
    int starting_index_part_in_domain = 0;
    double args_init_function[1] = {0.0}; // value of the temperature to be imposed everywhere in the domain. This argument is passed to the "initFunction" routine which will specify the values of each particle quantity
    int size_values = 1; // scalar temperature field with a single component
    
    allParticles* particles_in_domain = create_particles_in_domain(nbParticles_domain, domain_lim, size_values, args_init_function, initFunction, starting_index_part_in_domain);
    
    // Creation of particles on the boundaries
    int nbBoundaries = 4;
    int nbParticles_boundaries[4] = {NPTS_BOUNDARIES, NPTS_BOUNDARIES, NPTS_BOUNDARIES, NPTS_BOUNDARIES}; // number of particles on each boundary (WARNING: same number of particles should be chosen for the moment)
    double boundaries_lim[4][4] = {{0.0, 0.0, 0.0, 1.0},
				   {0.0, 1.0, 1.0, 1.0},
				   {1.0, 1.0, 0.0, 1.0},
				   {0.0, 1.0, 0.0, 0.0}};
    int starting_index_part_on_bound = NPTS_DOMAIN;
    double values_Dirichlet[4] = {0.0, 0.0, 0.0, 100.0}; 
    
    allParticles* particles_on_boundaries = create_particles_on_boundaries(nbParticles_boundaries, (double *)boundaries_lim, size_values, nbBoundaries, values_Dirichlet, starting_index_part_on_bound);
    
    // Assemble the particles in the domain and on the boundaries in a single structure, for the computation of the neighbourhoods just after
    allParticles* particles_everywhere = combine_two_particles_sets(particles_in_domain, particles_on_boundaries, IN_DOMAIN, ON_BOUNDARY);
    
    GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
    CHECK_MALLOC(data);
    double extrema[2] = {0.0, 100.0};
    myFillData(data, particles_everywhere, extrema);
//     draw_particles(domain_lim, data);
    
    
    // *** 2 *** CREATION OF THE NEIGHBOURHOODS OF EACH PARTICLE (in the domain and on the boundaries)
    
    // Creation of neighborhoods
    double timestep = 0.0;
    double maxspeed = 0.0;
    neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
    neighborhood* nh = options->nh;
    neighborhood_update_new(options, nh, particles_everywhere, 0);
    
    int test = particles_everywhere->array_of_particles[25].particle_neighbours->nh->nNeighbours;
    for(int i=0; i<test; i++) {
	int index = particles_everywhere->array_of_particles[25].particle_neighbours->nh[i].index;
// 	printf("index = %d \n", index);
    }
//     // Associate the computed neighbourhoods to each existing particle
//     associate_neighborhood_to_particles(particles_everywhere, nh);
    
    
    // *** 3 *** TIME LOOP TO RESOLVE THE EQUATIONS
//      bov_window_t* window = NULL;
//      bov_points_t *particles = NULL;
//      create_window_animation(data, window, particles);
    
    double alpha = 1.0;
    double time = 0.0;
    double time_max = 0.3;
    double dt = 1.0E-5;
    int nb_time_step = 0;
    int nb_time_step_max = time_max / dt;
    int print_every_time_step = 1000;
    int size_solution_vector = (int)(nb_time_step_max / print_every_time_step) + 1;
    double* solution_vector = malloc(size_values*size_solution_vector*sizeof(double));
    
    while (time < time_max) {
      // Store the solution at every n-time steps in a vector for visualization afterwards
      if (nb_time_step == nb_time_step_max) { //% print_every_time_step == 0) {
// 	store_solution(solution_vector, particles_everywhere);
	printf("------ Time step %d over %d -------\n", nb_time_step, nb_time_step_max+1);
	myFillData(data, particles_everywhere, extrema);
// 	double test = particles_everywhere->array_of_particles[25].values[0];
// 	double test = particles_everywhere->array_of_particles[25].particle_derivatives->laplacian[0];
	int test = particles_everywhere->array_of_particles[25].particle_neighbours->nh->nNeighbours;
// 	printf("test = %d\n",test);
//         display_particles(data, window, particles);
	draw_particles(domain_lim, data);
      }
      // Compute derivatives of all the particles inside the domain, in particular the Laplacian here
      computeDerivativesAllParticles(particles_everywhere, CUBIC);
      // Update the temperature of every particles with a chosen time integration scheme
      integrate_equation(particles_everywhere, dt, alpha);
//       // 
//       
      time += dt;
      nb_time_step++;
      // NOTE: No need to update the positions of the particles or to update the neighbourhoods since the particles are fixed
    }
    
    
    
    // Draw particles with color related to the derivatives of the quantity they carry
//     GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
//     CHECK_MALLOC(data);
//     myFillData(data, particles_everywhere->array_of_particles);
//     draw_particles(domain_lim, data);

    neighborhood_options_delete(options,nh);

// 	double timestep = 0.5;
// 	double maxspeed = 1;
// 	neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
// 	neighborhood* nh = options->nh;
// 	int number_of_iterations = 10;
// 	for (int iterations = 0; iterations < number_of_iterations;iterations++) {
// 		if(iterations)
// 			bouncyrandomupdate(data, timestep, options->half_length, maxspeed);
// 		neighborhood_update(options, nh, data, iterations);
// 		//kernel(data, nh, kh);
// 	}
// 	neighborhood_options_delete(options,nh);

    free(data);
    delete_all_particles(particles_in_domain);
    delete_all_particles(particles_on_boundaries);
    delete_all_particles(particles_everywhere);
    return EXIT_SUCCESS;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
