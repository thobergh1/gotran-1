#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>


// #include "beeler_reuter_1977_version06.hpp"
// #include "FHN.hpp"
// #include "GPB.hpp"
// #include "grandi_pasqualini_bers_2010.hpp"
// #include "Hinch_et_al_2004.hpp"
// #include "hodgkin_huxley.hpp"
// #include "irvine_model_1999.hpp"
// #include "iyer_2004.hpp"
// #include "JT21.hpp"
// #include "maleckar_2009.hpp"
// #include "maltsev_2009.hpp"
// #include "niederer_2006.hpp"
// #include "Niederer_et_al_2006.hpp"
// #include "pandit_et_al_2001_endo.hpp"
// #include "Pandit_Hinch_Niederer_Model.hpp"
// #include "rice_model_2008.hpp"
// #include "shannon_2004.hpp"
// #include "tentusscher_2004_mcell.hpp"
#include "tentusscher_panfilov_2006_M_cell.hpp"
// #include "terkildsen_2008.hpp"
// #include "winslow_model_1999.hpp"



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

#include <iostream>
using namespace std;

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
                             int num_timesteps, int num_cells, double dt)
{
  double t;
  int save_it = 1;
  int it, j;


  //int num_cells = 1;
  //int num_cells = 10;
  

  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  std::cout << num_timesteps << " num steps" << std::endl;
  std::cout << num_cells << " num cells" << std::endl;
  std::cout << padded_num_cells << " padded num cells" << std::endl;


  //double solve_dt = 0.01;
  double V_min = -100;
  double V_max = 100;
  double V_step = 0.0005;
  
  
  // double Ca_min = 0.00001;
  // double Ca_max = 10;
  // double Ca_step = 0.00001;
  

  const std::vector<univariate_func> expressions_V = *model_lut.expressions_V;
  default_LUT_type lut_V = LinearInterpolationLUT(V_min, V_max, V_step, expressions_V, expressions_V.size(), dt, parameters);

  // const std::vector<univariate_func> expressions_Ca_ss = *model_lut.expressions_Ca_ss;
  // default_LUT_type lut_Ca_ss = LinearInterpolationLUT(Ca_min, Ca_max, Ca_step, expressions_Ca_ss, expressions_Ca_ss.size(), dt, parameters);

 
  printf("\n");
  // printf("Num expressions (V) : %lu\n\n", expressions_V.size());
  //printf("Num expressions (Ca_ss) : %lu\n\n", expressions_Ca_ss.size());

  

  for (it = 1; it <= num_timesteps; it++) {
    t = t_values[it-1];
    // forward_explicit_euler(u, t, dt, parameters, num_cells, padded_num_cells);
    forward_explicit_euler(u, t, dt, parameters, num_cells, padded_num_cells, lut_V);
    // forward_explicit_euler(u, t, dt, parameters, num_cells, padded_num_cells, lut_V, lut_Ca_ss);

    //printf("u: %f t: %f dt: %f param: %f\n", *u, t, dt, *parameters);
    
    for (j=0; j < NUM_STATES; j++) {
      //u_values[save_it*NUM_STATES + j] = u[j*padded_num_cells];
      u_values[save_it*NUM_STATES + j] = u[j*padded_num_cells];

      // std::cout << u[j*padded_num_cells] << " ";
    }
    // std::cout << std::endl;
    //printf("\n");
    save_it++;
  }
}


extern "C"
void ode_solve_rush_larsen(double* u, cellmodel_float_t* parameters,
                             double* u_values, double* t_values,
                             int num_timesteps, int num_cells, double dt)
{
 double t;
  int save_it = 1;
  int it, j;


  //int num_cells = 1;
  //int num_cells = 10;
  

  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  std::cout << num_cells << " num cells" << std::endl;
  std::cout << padded_num_cells << " padded num cells" << std::endl;

  //double solve_dt = 0.01;
  double V_min = -100;
  double V_max = 100;
  double V_step = 0.001;
  
  
  // double Ca_min = 0.00001;
  // double Ca_max = 10;
  // double Ca_step = 0.00001;
  

  const std::vector<univariate_func> expressions_V = *model_lut.expressions_V;
  default_LUT_type lut_V = LinearInterpolationLUT(V_min, V_max, V_step, expressions_V, expressions_V.size(), dt, parameters);

  // const std::vector<univariate_func> expressions_Ca_ss = *model_lut.expressions_Ca_ss;
  // default_LUT_type lut_Ca_ss = LinearInterpolationLUT(Ca_min, Ca_max, Ca_step, expressions_Ca_ss, expressions_Ca_ss.size(), dt, parameters);

 
  printf("\n");
  // printf("Num expressions (V) : %lu\n\n", expressions_V.size());
  //printf("Num expressions (Ca_ss) : %lu\n\n", expressions_Ca_ss.size());

  

  for (it = 1; it <= num_timesteps; it++) {
    t = t_values[it-1];
    // forward_rush_larsen(u, t, dt, parameters, num_cells, padded_num_cells);
    forward_rush_larsen(u, t, dt, parameters, num_cells, padded_num_cells, lut_V);
    // forward_rush_larsen(u, t, dt, parameters, num_cells, padded_num_cells, lut_V, lut_Ca_ss);

    //printf("u: %f t: %f dt: %f param: %f\n", *u, t, dt, *parameters);
    //std::cout << u[0] << u[1] << u[2] << u[3] << std::endl;
    
    for (j=0; j < NUM_STATES; j++) {
      //u_values[save_it*NUM_STATES + j] = u[j*padded_num_cells];
      u_values[save_it*NUM_STATES + j] = u[j*padded_num_cells];

      //std::cout << u[j] << std::endl;
    }
    //printf("\n");
    save_it++;
  }
}


extern "C"
void ode_solve_generalized_rush_larsen(double* u, cellmodel_float_t* parameters,
                             double* u_values, double* t_values,
                             int num_timesteps, int num_cells, double dt)
{
  double t;
  int save_it = 1;
  int it, j;


  //int num_cells = 1;
  //int num_cells = 10;
  

  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  std::cout << num_timesteps << " num steps" << std::endl;
  std::cout << num_cells << " num cells" << std::endl;
  std::cout << padded_num_cells << " padded num cells" << std::endl;


  //double solve_dt = 0.01;
  double V_min = -100;
  double V_max = 100;
  double V_step = 0.0005;
  
  
  // double Ca_min = 0.00001;
  // double Ca_max = 10;
  // double Ca_step = 0.00001;
  

  const std::vector<univariate_func> expressions_V = *model_lut.expressions_V;
  default_LUT_type lut_V = LinearInterpolationLUT(V_min, V_max, V_step, expressions_V, expressions_V.size(), dt, parameters);

  // const std::vector<univariate_func> expressions_Ca_ss = *model_lut.expressions_Ca_ss;
  // default_LUT_type lut_Ca_ss = LinearInterpolationLUT(Ca_min, Ca_max, Ca_step, expressions_Ca_ss, expressions_Ca_ss.size(), dt, parameters);

 
  printf("\n");
  // printf("Num expressions (V) : %lu\n\n", expressions_V.size());
  //printf("Num expressions (Ca_ss) : %lu\n\n", expressions_Ca_ss.size());

  

  for (it = 1; it <= num_timesteps; it++) {
    t = t_values[it-1];
    // forward_generalized_rush_larsen(u, t, dt, parameters, num_cells, padded_num_cells);
    forward_generalized_rush_larsen(u, t, dt, parameters, num_cells, padded_num_cells, lut_V);
    // forward_generalized_rush_larsen(u, t, dt, parameters, num_cells, padded_num_cells, lut_V, lut_Ca_ss);

    //printf("u: %f t: %f dt: %f param: %f\n", *u, t, dt, *parameters);
    
    for (j=0; j < NUM_STATES; j++) {
      //u_values[save_it*NUM_STATES + j] = u[j*padded_num_cells];
      u_values[save_it*NUM_STATES + j] = u[j*padded_num_cells];

      // std::cout << u[j*padded_num_cells] << " ";
    }
    // std::cout << std::endl;
    //printf("\n");
    save_it++;
  }
}


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
  double dt = 1E-3;
  // double dt = 1E-7;

  //double dt = 0.02E-3;
  
  
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

  //int num_cells = 11500000;
  int num_cells = 1;
  //int num_cells = 10;
  
  size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
  unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
          num_cells, alignment_bytes / sizeof(cellmodel_float_t));

  size_t states_size = num_states * sizeof(double) *padded_num_cells;
  size_t parameters_size = num_parameters * sizeof(double);

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
  
 
  double V_min = -100;
  double V_max = 100;
  double V_step = 0.001;

  
  // double Ca_min = 0.00001;
  // double Ca_max = 10;
  // double Ca_step = 0.00001;


  const std::vector<univariate_func> expressions_V = *model_lut.expressions_V;
  default_LUT_type lut_V = LinearInterpolationLUT(V_min, V_max, V_step, expressions_V, expressions_V.size(), dt, parameters);

  //const std::vector<univariate_func> expressions_Ca_ss = *model_lut.expressions_Ca_ss;
  //default_LUT_type lut_Ca_ss = LinearInterpolationLUT(Ca_min, Ca_max, Ca_step, expressions_Ca_ss, expressions_Ca_ss.size(), dt, parameters);

 
  printf("\n");
  // printf("Num expressions (V) : %lu\n\n", expressions_V.size());
  //printf("Num expressions (Ca_ss) : %lu\n\n", expressions_Ca_ss.size());
  


  // forward euler
  printf("Scheme: Forward Euler\n");
  int it;

  /*

  std::string folder_path = "/home/thomas/Dokumenter/gotran/venv/gotran/demo/gotran2c/results/time";
  std::string file_name = "time_lut_V_" + std::to_string(num_timesteps) + "_steps.txt";
  std::string file_path = folder_path + "/" + file_name;

  std::ofstream outputFile(file_path);


  for(int timing=0; timing < 10; timing++){
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_start);
    for (it = 0; it < num_timesteps; it++) {
      // forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells);
      forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells, lut_V);
      //forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells, lut_V, lut_Ca_ss);
      t += dt;
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_now);

    time_elapsed = timestamp_now.tv_sec - timestamp_start.tv_sec + 1E-9 * (timestamp_now.tv_nsec - timestamp_start.tv_nsec);
    outputFile << std::fixed << std::setprecision(8) << time_elapsed << std::endl;
    //std::cout << std::setprecision(8) << time_elapsed << std::endl;
  }
  outputFile.close();
  */

  // unsigned int half = int(num_states/2);

  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_start);
  for (it = 0; it < num_timesteps; it++) {
    // forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells);
    forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells, lut_V);
    //forward_explicit_euler(states, t, dt, parameters, num_cells, padded_num_cells, lut_V, lut_Ca_ss);

    // for(unsigned int j=0; j < num_states; j++){
    //   if (j == half) {
    //     std::cout << endl;
    //   }
    //   else{
    //     std::cout << std::fixed << std::setprecision(8) << states[j*padded_num_cells] << " ";
    //   }
    // }
    // std::cout << "\n" << std::endl;

    t += dt;
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_now);
  time_elapsed = timestamp_now.tv_sec - timestamp_start.tv_sec + 1E-9 * (timestamp_now.tv_nsec - timestamp_start.tv_nsec);



  printf("Computed %d time steps in %g s. Time steps per second: %g\n",
      num_timesteps, time_elapsed, num_timesteps/time_elapsed);
  printf("\n");
  // printf("V= %f at t=%f\n", states[7*padded_num_cells],t);
  
  /*
  time_elapsed = timestamp_now.tv_sec - timestamp_start.tv_sec + 1E-9 * (timestamp_now.tv_nsec - timestamp_start.tv_nsec);
  printf("Computed %d time steps in %g s. Computed cell steps per second: %g\n",
      num_timesteps, time_elapsed, num_cells*num_timesteps/time_elapsed);
  */

  printf("\n");



  // Rush Larsen
  printf("Scheme: Rush Larsen (exp integrator on all gates)\n");
  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_start);
  init_state_values(states, num_cells, padded_num_cells);
  for (it = 0; it < num_timesteps; it++) {
    // forward_rush_larsen(states, t, dt, parameters, num_cells, padded_num_cells);
    forward_rush_larsen(states, t, dt, parameters, num_cells, padded_num_cells, lut_V);
    t += dt;
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_now);
  time_elapsed = timestamp_now.tv_sec - timestamp_start.tv_sec + 1E-9 * (timestamp_now.tv_nsec - timestamp_start.tv_nsec);
   printf("Computed %d time steps in %g s. Computed cell steps per second: %g\n",
      num_timesteps, time_elapsed, num_cells*num_timesteps/time_elapsed);
  printf("\n");
  // printf("V= %f at t=%f\n", states[17*padded_num_cells],t);

  free(states);
  free(parameters);
 
  return 0;

}
