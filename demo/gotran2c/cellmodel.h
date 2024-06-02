#ifndef _CELLMODEL_H
#define _CELLMODEL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

enum scheme { SCHEME_FE, SCHEME_RL, SCHEME_GRL1, NUM_SCHEMES };
enum memory_layout { LAYOUT_STRUCT_OF_ARRAYS, LAYOUT_ARRAY_OF_STRUCTS };
typedef enum scheme scheme_type;

//#define CELLMODEL_FLOAT_USE_SINGLE_PRECISION

#ifdef CELLMODEL_FLOAT_USE_SINGLE_PRECISION

#define Exp(x) expf(x)
#ifdef USE_EXPM1
#define Expm1(x) expm1f(x)
#else
#define Expm1(x) (expf(x) - 1.f)
#endif
#define Log(x) logf(x)
#define FP_LITERAL(num) num##f
typedef float cellmodel_float_t;

#else

#define Exp(x) exp(x)
#ifdef USE_EXPM1
#define Expm1(x) expm1(x)
#else
#define Expm1(x) (exp(x) - 1)
#endif
#define Log(x) log(x)
#define FP_LITERAL(num) num
typedef double cellmodel_float_t;

#endif

#ifdef VECTOR_LENGTH
#define MULTISTEP_VECTOR_LENGTH VECTOR_LENGTH
#else
#define MULTISTEP_VECTOR_LENGTH 8
#endif

//#define CELLMODEL_STATES_ALIGNMENT_BYTES 256
#if defined(VECTOR_LENGTH) && VECTOR_LENGTH > MULTISTEP_VECTOR_LENGTH
#define CELLMODEL_STATES_ALIGNMENT_BYTES (VECTOR_LENGTH * sizeof(cellmodel_float_t))
#else
#define CELLMODEL_STATES_ALIGNMENT_BYTES (MULTISTEP_VECTOR_LENGTH * sizeof(cellmodel_float_t))
#endif

struct state_colour_sets {
    const int num_sets;
    const uint8_t *
            set_map; // array of length num states, entries indicate which colour set the state belongs to
};
// typedef void (*init_states_func)(cellmodel_float_t *states, long num_cells, long padded_num_cells);
typedef void (*init_states_func)(cellmodel_float_t *states);


typedef void (*init_parameters_func)(cellmodel_float_t *parameters);
typedef int (*state_index_func)(const char name[]);
typedef int (*parameter_index_func)(const char name[]);
typedef void (*step_func)(cellmodel_float_t *states, const cellmodel_float_t t,
                          const cellmodel_float_t dt, const cellmodel_float_t *parameters,
                          long num_cells, long padded_num_cells);
typedef void (*multistep_func)(cellmodel_float_t *states, const cellmodel_float_t t_start,
                               const cellmodel_float_t dt, const cellmodel_float_t *parameters,
                               long num_cells, long padded_num_cells, int num_timesteps);

struct cellmodel {
    init_states_func init_states;
    init_parameters_func init_parameters;
    state_index_func state_index;
    parameter_index_func parameter_index;
    step_func step_FE;
    step_func step_RL;
    step_func step_GRL1;
    multistep_func multistep_FE;
    multistep_func multistep_RL;
    multistep_func multistep_GRL1;
    int num_states;
    int num_parameters;
    enum memory_layout layout;
    const struct state_colour_sets *colour_sets;
    const char **state_names;
};

#ifdef __cplusplus
} // extern "C"
#endif

#endif
