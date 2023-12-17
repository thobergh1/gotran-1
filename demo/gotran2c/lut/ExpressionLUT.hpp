#ifndef _LUT__EXPRESSIONLUT_HPP
#define _LUT__EXPRESSIONLUT_HPP

#include <vector>

#include "LUT.hpp"

class ExpressionLUT : public LUT
{
  protected:
    const int num_expressions;
    inline int table_index(int V_index, int expr_index) {
      return V_index * num_expressions + expr_index;
    }

  public:
    ExpressionLUT(double V_min, double V_max, double V_step,
                  std::vector<univariate_func> expressions, int _num_expressions, double dt,
                  double *param);
    ~ExpressionLUT();
};

#endif
