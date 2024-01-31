#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "hodgkin_huxley.hpp"

//#include "tentusscher_panfilov_2006_M_cell.hpp"
//#include "FHN.h"

#include <assert.h>
#include <stdint.h>
#include "lut/LUT.hpp"
#include "lut/LinearInterpolationLUT.hpp"

#include "lut/LinearInterpolationLUT.cpp"
#include "lut/ExpressionLUT.hpp"
#include "lut/ExpressionLUT.cpp"


// #include "bench.h"
// #include "bench_util.h"
// #include "rng.h"
// #include "schemes.h"


// Gotran generated C/C++ code for the "base_model" model
uint64_t ceil_div_uint64(uint64_t a, uint64_t b);
uint64_t ceil_to_multiple_uint64(uint64_t a, uint64_t b);

uint64_t ceil_div_uint64(uint64_t a, uint64_t b)
{
    return (a + b - 1) / b;
}

uint64_t ceil_to_multiple_uint64(uint64_t a, uint64_t b)
{
    return ceil_div_uint64(a, b) * b;
}

extern "C"
int detect_num_threads()
{
    int count = 0;
    #pragma omp parallel
    {
        #pragma omp critical
        {
            count++;
        }
    }
    return count;
}

extern "C"
void ode_solve_forward_euler(double* u, cellmodel_float_t* parameters,
                             double* u_values, double* t_values,
                             int num_timesteps, double dt)
{
  double t;
  int save_it = 1;
  int it, j;


  int num_cells = 10;
  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));


  unsigned int num_states = NUM_STATES;
  unsigned int num_parameters = NUM_PARAMS;
  size_t states_size = num_states * sizeof(double) *padded_num_cells;
  //size_t parameters_size = num_parameters * sizeof(double);


  cellmodel_float_t *states = (cellmodel_float_t *) aligned_alloc(alignment_bytes, states_size);
  //cellmodel_float_t *parameters = (cellmodel_float_t *) malloc(parameters_size);


  
  const std::vector<univariate_func> expressions_V = *model_lut.expressions_V;

  printf("\n");
  printf("Num expressions (V): %lu\n\n", expressions_V.size());


  double T_end = 40;
  double solve_dt = 0.01;
  int store_period = 1;
  double V_min = -100;
  double V_max = 100;
  double V_step = 0.05;  

  default_LUT_type lut_V = LinearInterpolationLUT(V_min, V_max, V_step, expressions_V, expressions_V.size(), solve_dt, parameters);


  /*

  double lookup_value;
  for(int V = 0; V < 100; V++){
    auto lut_V_state = lut_V.compute_input_state(V);
    for(int i = 0; i<100; i++){
      lookup_value = lut_V.lookup(i, lut_V_state);
      printf("%g\n", lookup_value);}}
  */



  for (it = 1; it <= num_timesteps; it++) {
    t = t_values[it-1];
    //forward_explicit_euler(u, t, dt, parameters, num_cells, padded_num_cells);
    forward_explicit_euler(u, t, dt, parameters, num_cells, padded_num_cells, lut_V);

    //printf("u: %f t: %f dt: %f param: %f\n", *u, t, dt, *parameters);

    //std::cout << u[0] << u[1] << u[2] << u[3] << std::endl;
    
    for (j=0; j < NUM_STATES; j++) {
      u_values[save_it*NUM_STATES + j] = u[j*padded_num_cells];
      //std::cout << u[j] << std::endl;
    }
    //printf("\n");
    save_it++;
  }
}

#if 0
void ode_solve_rush_larsen(double* u, const double* parameters,
                             double* u_values, double* t_values,
                             int num_timesteps, double dt)
{
  double t;
  int save_it = 1;
  int it, j;

  int num_cells = 8;
  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  unsigned int num_states = NUM_STATES;
  size_t states_size = num_states * sizeof(double) * padded_num_cells;

  unsigned int num_parameters = NUM_PARAMS;
  size_t parameters_size = num_parameters * sizeof(double);

  for (it = 1; it <= num_timesteps; it++) {
    t = t_values[it-1];
    forward_rush_larsen(u, t, dt, parameters, num_cells, padded_num_cells);
    for (j=0; j < NUM_STATES; j++) {
      u_values[save_it*NUM_STATES + j] = u[j];
    }
    save_it++;
  }
}
#endif

extern "C"
int state_count()
{
  return NUM_STATES;
}

extern "C"
int parameter_count()
{
  return NUM_PARAMS;
}

extern "C"
int main(int argc, char *argv[])
{
  double t_start = 0;
  //double dt = 0.0001;
  double dt = 0.01;
  
  
  int num_timesteps = (int) 5;
  if (argc > 1) {
    num_timesteps = atoi(argv[1]);
    printf("num_timesteps set to %d\n", num_timesteps);
    if(num_timesteps <= 0) {
        exit(EXIT_FAILURE);
    }
  }

  unsigned int num_states = NUM_STATES;
  unsigned int num_parameters = NUM_PARAMS;

  int num_cells = 11500000;
  //int num_cells = 50;
  
  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  size_t states_size = num_states * sizeof(double) *padded_num_cells;
  size_t parameters_size = num_parameters * sizeof(double);

  //cellmodel_float_t *states = aligned_alloc(alignment_bytes, states_size);
  //cellmodel_float_t *parameters = malloc(parameters_size);

  cellmodel_float_t *states = (cellmodel_float_t *) aligned_alloc(alignment_bytes, states_size);
  cellmodel_float_t *parameters = (cellmodel_float_t *) malloc(parameters_size);


  init_parameters_values(parameters);
  init_state_values(states, num_cells, padded_num_cells);

  double t = t_start;

  struct timespec timestamp_start, timestamp_now;
  double time_elapsed;


  //printf("num_timesteps: %ld\n", num_timesteps);
  int num_threads = detect_num_threads();
  printf("Using %d thread(s)\n", num_threads);

  size_t total_size = states_size + parameters_size;
  printf("Memory footprint: %lu bytes\n\n", total_size);
  


  //const std::vector<univariate_func> expressions_V = expressions_V;

  /*
  for (const auto& expr : expressions_V) {
      printf("Name: %s\n", expr.str);
      printf("Expression: %f\n", expr.f(1.0, 1.0, nullptr));
      printf("----------------------\n");
  }
  */

 
  const std::vector<univariate_func> expressions_V = *model_lut.expressions_V;

  printf("\n");
  printf("Num expressions (V) : %lu\n\n", expressions_V.size());

  


  double T_end = 1000;
  double solve_dt = 1E-3;
  int store_period = 1;

  double V_min = -100;
  double V_max = 100;
  double V_step = 0.05;
  
  default_LUT_type lut_V = LinearInterpolationLUT(V_min, V_max, V_step, expressions_V, expressions_V.size(), solve_dt, parameters);

 

  // forward euler
  printf("Scheme: Forward Euler\n");
  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_start);
  int it;
  for (it = 0; it < num_timesteps; it++) {
    //printf("%ld\n", it);

    //forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells);
    forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells, lut_V);
    
    t += dt;
  }

  

  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_now);
  time_elapsed = timestamp_now.tv_sec - timestamp_start.tv_sec + 1E-9 * (timestamp_now.tv_nsec - timestamp_start.tv_nsec);
  printf("Computed %d time steps in %g s. Computed cell steps per second: %g\n",
      num_timesteps, time_elapsed, num_cells*num_timesteps/time_elapsed);
  printf("\n");



  #if 0
  // Rush Larsen
  printf("Scheme: Rush Larsen (exp integrator on all gates)\n");
  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_start);
  init_state_values(states, num_cells, padded_num_cells);
  for (it = 0; it < num_timesteps; it++) {
    forward_rush_larsen(states, t, dt, parameters, num_cells, padded_num_cells);
    t += dt;
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_now);
  time_elapsed = timestamp_now.tv_sec - timestamp_start.tv_sec + 1E-9 * (timestamp_now.tv_nsec - timestamp_start.tv_nsec);
   printf("Computed %d time steps in %g s. Computed cell steps per second: %g\n",
      num_timesteps, time_elapsed, num_cells*num_timesteps/time_elapsed);
  printf("\n");

  printf("V= %f at t=%f\n", states[17*padded_num_cells],t);

  free(states);
  free(parameters);

  return 0;
  #endif
}