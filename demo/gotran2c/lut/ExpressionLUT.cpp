#include "ExpressionLUT.hpp"

ExpressionLUT::ExpressionLUT(double V_min, double V_max, double V_step,
                             const std::vector<univariate_func> expressions, int _num_expressions,
                             double dt, double *param)
    : LUT(V_min, V_max, V_step), num_expressions(_num_expressions)
{
    //printf("Allocating LUT\n");
    int num_columns = num_expressions;
    int num_V_values = this->num_V_values;
    this->table = new double[num_columns * num_V_values];
    size = num_columns * num_V_values * sizeof(double);

    // populate table
    for (int j = 0; j < num_expressions; j++) {
        univariate_func expr = expressions[j];
        double V = V_min;
        for (int i = 0; i < num_V_values; i++) {
            double val = expr(V, dt, param);
            this->table[table_index(i,j)] = val;
            V += this->V_step;
        }
    }
}

ExpressionLUT::~ExpressionLUT()
{
    //printf("Deallocating LUT\n");
    delete[] table;
}
