#ifndef _LUT__LINEARINTERPOLATIONLUT_HPP
#define _LUT__LINEARINTERPOLATIONLUT_HPP

#include <assert.h>

#include <vector>

#include "ExpressionLUT.hpp"

struct LinearInterpolationInputState {
    double weight_above;
    double weight_below;
    int row_above;
    int row_below;
};


class LinearInterpolationLUT : public ExpressionLUT
{
  private:
    int row_above;
    int row_below;
    double weight_above;
    double weight_below;

  public:
    const char* class_desc = "Linear interpolation (row-major)";

    LinearInterpolationLUT(double V_min, double V_max, double V_step,
                           const std::vector<univariate_func> expressions, int _num_expressions,
                           double dt, double *param)
        : ExpressionLUT(V_min, V_max, V_step, expressions, _num_expressions, dt, param)
    {
    }
    bool set_input(double V);

    LinearInterpolationInputState compute_input_state(double V);

    inline double lookup(int expression_id)
    {
#ifdef DEBUG
        assert(expression_id < this->num_expressions);
#endif
        double val = this->weight_above * this->table[table_index(this->row_above, expression_id)]
                     + this->weight_below * this->table[table_index(this->row_below, expression_id)];
        return val;
    }

    inline double lookup(int expression_id, const LinearInterpolationInputState &state)
    {
#ifdef DEBUG
        assert(expression_id < this->num_expressions);
#endif
        double val = state.weight_above * this->table[table_index(state.row_above, expression_id)]
                     + state.weight_below * this->table[table_index(state.row_below, expression_id)];
        return val;
    }
};

#endif
