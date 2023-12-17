#ifndef _LUT__LUT_HPP
#define _LUT__LUT_HPP

#include <math.h>
#include <stddef.h>

typedef double (*univariate_func)(double x, double dt, double *param);

struct univariate_func_tuple {
    const char *str;
    const univariate_func f;
    constexpr univariate_func_tuple(const char *str, univariate_func f) : str(str), f(f)
    {
    }
};

class LUT
{
  protected:
    double V_min;
    double V_max;
    double V_step;
    double V_step_inv;
    int num_V_values;
    int num_expressions;
    double *table = nullptr;
    size_t size = 0;

    double V = 0;
    bool V_is_set = false;

    bool check_input_value(double V)
    {
        return this->V_is_set;
    }

    void setup_V_limits(double V_min, double V_max, double V_step)
    {
        this->V_min = V_min;
        this->V_step = V_step;
        this->V_step_inv = 1. / V_step;

        this->num_V_values = ((int) ceil((V_max - V_min) / V_step)) + 1;
        this->V_max = V_min + (this->num_V_values - 1) * V_step;
    }

    LUT()
    {
    }

    LUT(double V_min, double V_max, double V_step)
    {
        this->setup_V_limits(V_min, V_max, V_step);
    }

  public:
    inline size_t get_size()
    {
        return size;
    }

    inline bool input_is_valid(double V)
    {
        return this->V_min <= V && V <= this->V_max;
    }

    virtual bool set_input(double V) = 0;

    virtual double lookup(int expression_id) = 0;
};

#endif // _LUT__LUT_HPP
