#ifndef _CELLMODEL_LUT_H
#define _CELLMODEL_LUT_H

#include <vector>

#include "cellmodel.h"
#include "lut/LUT.hpp"
#include "lut/LinearInterpolationLUT.hpp"

typedef LinearInterpolationLUT default_LUT_type;




typedef void (*step_lut_func)(cellmodel_float_t *states, const cellmodel_float_t t,
                              const cellmodel_float_t dt, const cellmodel_float_t *parameters,
                              default_LUT_type &lut_V);

// typedef void (*step_lut_func)(cellmodel_float_t *states, const cellmodel_float_t t,
//                               const cellmodel_float_t dt, const cellmodel_float_t *parameters,
//                               const long num_cells, const long padded_num_cells,
//                               default_LUT_type &lut_V);

// typedef void (*step_lut_func)(cellmodel_float_t *states, const cellmodel_float_t t,
//                               const cellmodel_float_t dt, const cellmodel_float_t *parameters,
//                               const long num_cells, const long padded_num_cells,
//                               default_LUT_type &lut_V, default_LUT_type &lut_Ca);


typedef void (*step_rhs_func)(const double *states,
                              double t,
                              const double *parameters,
                              double *values,
                              default_LUT_type &lut_V);

// typedef void (*step_rhs_func)(const double *states,
//                               double t,
//                               const double *parameters,
//                               double *values,
//                               long num_cells,
//                               long padded_num_cells,
//                               default_LUT_type &lut_V);

// typedef void (*step_rhs_func)(const double *states,
//                               double t,
//                               const double *parameters,
//                               double *values,
//                               long num_cells,
//                               long padded_num_cells,
//                               default_LUT_type &lut_V,
//                               default_LUT_type &lut_Ca);



// struct cellmodel_lut {
//     init_states_func init_states;
//     init_parameters_func init_parameters;
//     state_index_func state_index;
//     parameter_index_func parameter_index;
//     step_rhs_func step_rhs;
//     step_lut_func step_RL;
//     int num_states;
//     int num_parameters;
//     const std::vector<univariate_func> *expressions_V;
// };


struct cellmodel_lut {
    init_states_func init_states;
    init_parameters_func init_parameters;
    state_index_func state_index;
    parameter_index_func parameter_index;
    step_rhs_func step_rhs;
    step_lut_func step_FE;
    step_lut_func step_RL;
    step_lut_func step_GRL;
    int num_states;
    int num_parameters;
    const std::vector<univariate_func> *expressions_V;
};



// struct cellmodel_lut {
//     init_states_func init_states;
//     init_parameters_func init_parameters;
//     state_index_func state_index;
//     parameter_index_func parameter_index;
//     step_lut_func step_FE;
//     step_lut_func step_RL;
//     step_lut_func step_GRL1;
//     int num_states;
//     int num_parameters;
//     const std::vector<univariate_func> *expressions_V;
//     const std::vector<univariate_func> *expressions_Ca;
// };

#endif
 