#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "tentusscher_panfilov_2006_M_cell.h"

#include <assert.h>
#include <stdint.h>

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


void ode_solve_forward_euler(double* u, const double* parameters,
                             double* u_values, double* t_values,
                             int num_timesteps, double dt)
{
  double t;
  int save_it = 1;
  int it, j;

  // unsigned int num_states = NUM_STATES;
  // size_t states_size = num_states * sizeof(double);

  // unsigned int num_parameters = NUM_PARAMS;
  // size_t parameters_size = num_parameters * sizeof(double);

  int num_cells = 1;
  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  for (it = 1; it <= num_timesteps; it++) {
    t = t_values[it-1];
    forward_explicit_euler(u, t, dt, parameters, num_cells, padded_num_cells);
    //printf("u: %f t: %f dt: %f param: %f\n", *u, t, dt, *parameters);

    for (j=0; j < NUM_STATES; j++) {
      u_values[save_it*NUM_STATES + j] = u[j];
    }
    save_it++;
  }
}

void ode_solve_rush_larsen(double* u, const double* parameters,
                             double* u_values, double* t_values,
                             int num_timesteps, double dt)
{
  double t;
  int save_it = 1;
  int it, j;

  // unsigned int num_states = NUM_STATES;
  // size_t states_size = num_states * sizeof(double);

  // unsigned int num_parameters = NUM_PARAMS;
  // size_t parameters_size = num_parameters * sizeof(double);

  int num_cells = 1;
  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  for (it = 1; it <= num_timesteps; it++) {
    t = t_values[it-1];
    forward_rush_larsen(u, t, dt, parameters, num_cells, padded_num_cells);
    for (j=0; j < NUM_STATES; j++) {
      u_values[save_it*NUM_STATES + j] = u[j];
    }
    save_it++;
  }
}

int state_count()
{
  return NUM_STATES;
}

int parameter_count()
{
  return NUM_PARAMS;
}

int main(int argc, char *argv[])
{
  double t_start = 0;
  double dt = 0.02E-3;
  int num_timesteps = (int) 1000000;
  if (argc > 1) {
    num_timesteps = atoi(argv[1]);
    printf("num_timesteps set to %d\n", num_timesteps);
    if(num_timesteps <= 0) {
        exit(EXIT_FAILURE);
    }
  }

  unsigned int num_states = NUM_STATES;
  unsigned int num_parameters = NUM_PARAMS;

  int num_cells = 1;
  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  size_t states_size = num_states * sizeof(double);
  size_t parameters_size = num_parameters * sizeof(double);

  cellmodel_float_t *states = aligned_alloc(alignment_bytes, states_size);
  cellmodel_float_t *parameters = malloc(parameters_size);

  init_parameters_values(parameters);
  init_state_values(states, num_cells, padded_num_cells);

  double t = t_start;

  struct timespec timestamp_start, timestamp_now;
  double time_elapsed;


  //printf("num_timesteps: %ld\n", num_timesteps);


  // forward euler
  printf("Scheme: Forward Euler\n");
  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_start);
  int it;
 
  for (it = 0; it < num_timesteps; it++) {
    //printf("%ld\n", it);

    forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells);


    
    t += dt;
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_now);
  time_elapsed = timestamp_now.tv_sec - timestamp_start.tv_sec + 1E-9 * (timestamp_now.tv_nsec - timestamp_start.tv_nsec);
  printf("Computed %d time steps in %g s. Time steps per second: %g\n",
      num_timesteps, time_elapsed, num_timesteps/time_elapsed);
  printf("\n");


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
  printf("Computed %d time steps in %g s. Time steps per second: %g\n",
      num_timesteps, time_elapsed, num_timesteps/time_elapsed);
  printf("\n");


  free(states);
  free(parameters);

  return 0;
}
