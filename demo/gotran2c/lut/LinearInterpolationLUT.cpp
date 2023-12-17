#include <assert.h>
#include <stdio.h>

#include "LinearInterpolationLUT.hpp"


bool LinearInterpolationLUT::set_input(double V)
{
    if (V < this->V_min || this->V_max < V) {
        return false;
    }

    this->V = V;
    this->V_is_set = true;

    double row_fractional = (V - this->V_min) / this->V_step;
    this->row_above = (int) floor(row_fractional);
    this->row_below = (int) ceil(row_fractional);

    this->weight_above = 1 - (row_fractional - row_above);
    this->weight_below = 1 - weight_above;

#if 0
    printf("V=%g yields row_above = %d, row_below = %d, and weight_above = %g and weight_below=%g\n",
        V, row_above, row_below, weight_above, weight_below);
#endif

    return true;
}

LinearInterpolationInputState LinearInterpolationLUT::compute_input_state(double V)
{
    bool valid = input_is_valid(V);
    if (!valid) {
        printf("V = %g exceeds the LUT range [%g, %g]\n", V, this->V_min, this->V_max);
        assert(input_is_valid(V));
    }

    double row_fractional = (V - this->V_min) * this->V_step_inv;

    int row_above;
    int row_below;
    double weight_above;
    double weight_below;

    row_above = (int) floor(row_fractional);
    row_below = row_above+1;

    weight_above = 1 - (row_fractional - row_above);
    weight_below = 1 - weight_above;

    LinearInterpolationInputState s;
    s.row_above = row_above;
    s.row_below = row_below;
    s.weight_above = weight_above;
    s.weight_below = weight_below;

    return s;
}
