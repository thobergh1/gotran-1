"""test for odecomponent module"""

# Imports for evaluation of generated code
from __future__ import division
import numpy as np
import math

import unittest
import gotran
import sys

from modelparameters.logger import suppress_logging
from modelparameters.codegeneration import sympycode
from modelparameters.sympytools import sp_namespace, symbols_from_expr
from modelparameters.parameters import *
from modelparameters.utils import list_timings

globals().update(sp_namespace)

from gotran.common import GotranException, parameters

from gotran.model.odeobjects2 import *
from gotran.model.ode2 import *
from gotran.model.loadmodel2 import *
from gotran.model.expressions2 import *
from gotran.codegeneration.algorithmcomponents import *
from gotran.codegeneration.codecomponent import *
from gotran.codegeneration.codegenerator2 import *


suppress_logging()

#default_params = parameters["code_generation"].copy()

def get_indexed(comp, name):
    check_arg(comp, CodeComponent)
    return [expr for expr in comp.body_expressions \
            if isinstance(expr, IndexedExpression) and expr.basename == name]

# Load reference CodeComponents
ode = load_ode("tentusscher_2004_mcell_updated.ode")

# Python code generator with default code generation paramters
codegen = PythonCodeGenerator()

# Options for code generation
default_params = parameters["code_generation"]

state_repr_opts = sorted(dict.__getitem__(default_params.states, "representation")._options)
param_repr_opts = dict.__getitem__(default_params.parameters, "representation")._options
body_repr_opts = dict.__getitem__(default_params.body, "representation")._options
body_optimize_opts = dict.__getitem__(default_params.body, "optimize_exprs")._options

module_code = codegen.module_code(ode)
exec(module_code)

states_values = init_state_values()
parameter_values = init_parameter_values()
rhs_ref_values = rhs(states_values, 0.0, parameter_values)
jac_ref_values = compute_jacobian(states_values, 0.0, parameter_values)

debug = 0

if debug:
    print "REF RHS:", rhs_ref_values
    print "REF JAC:", jac_ref_values

code_params = parameters["code_generation"]

test_map = {}

# A Function closure for the generation of test methods
def function_closure(body_repr, body_optimize, param_repr, \
                     state_repr, use_cse, float_precision, \
                     parameters_name="parameters", states_name="states", \
                     body_array_name="body", body_in_arg=False):

    # The test that will be attached to the TestCase class below
    def test(self):
        test_name = "body repr: {0}, body opt: {1}, state repr: {2}, "\
                    "parameter repr: {3}, use_cse: {4}, float_precision: "\
                    "{5}, state_name: {6}, param_name: {7}, body_name: {8}, "\
                    "body_in_arg: {9}".format(\
            body_repr, body_optimize, state_repr, param_repr, use_cse, \
            float_precision, states_name, parameters_name, body_array_name, body_in_arg)
        print "Testing code generation with parameters: " + test_name

        # Update main parameters
        code_params["body"]["optimize_exprs"] = body_optimize
        code_params["body"]["representation"] = body_repr
        code_params["body"]["array_name"] =  body_array_name
        code_params["parameters"]["representation"] = param_repr
        code_params["parameters"]["array_name"] = parameters_name
        code_params["states"]["representation"] = state_repr
        code_params["states"]["array_name"] = states_name
        code_params["body"]["use_cse"] = use_cse
        code_params["float_precision"] = float_precision
        code_params["default_arguments"] = "stp" + ("b" if body_in_arg else "")

        # Reload ODE for each test
        ode = load_ode("tentusscher_2004_mcell_updated.ode")
        codegen = PythonCodeGenerator()
        rhs_comp = rhs_expressions(ode)
        rhs_code = codegen.function_code(rhs_comp)
        exec rhs_code in globals(), locals()

        args = [states_values, 0.0]
        if param_repr != "numerals":
            args.append(parameter_values)

        if body_in_arg:
            body = np.zeros(rhs_comp.shapes[body_array_name])
            args.append(body)

        # Call the generated rhs function
        rhs_values = rhs(*args)
        rhs_norm = np.sqrt(np.sum(rhs_ref_values-rhs_values)**2)

        # DEBUG
        if debug:
            test_name = "_".join([body_repr, body_optimize, state_repr, param_repr, \
                                  float_precision, states_name, parameters_name, \
                                  body_array_name, str(body_in_arg)])
            test_name += "_use_cse_" + str(use_cse)

            open("rhs_code_{0}.py".format(test_name), "w").write(rhs_code)
            print "rhs norm:", rhs_norm

        eps = 1e-8 if float_precision == "double" else 1e-6
        self.assertTrue(rhs_norm<eps)

        # Only evaluate jacobian if not using cse
        #if use_cse:
        #    return

        jac_comp = jacobian_expressions(ode)
        jac_code = codegen.function_code(jac_comp)
        exec jac_code in globals(), locals()

        args = [states_values, 0.0]
        if param_repr != "numerals":
            args.append(parameter_values)

        if body_in_arg:
            body = np.zeros(jac_comp.shapes[body_array_name])
            args.append(body)

        jac_values = compute_jacobian(*args)
        
        jac_norm = np.sqrt(np.sum(jac_ref_values-jac_values)**2)
        if debug:
            test_name = "_".join([body_repr, body_optimize, state_repr, param_repr, \
                                  float_precision, states_name, parameters_name, \
                                  body_array_name, str(body_in_arg)])
            test_name += "_use_cse_" + str(use_cse)

            open("jac_code_{0}.py".format(test_name), "w").write(jac_code)
            print "jac norm:", jac_norm
            print "JAC", jac_values

        eps = 1e-8 if float_precision == "double" else 1e-3
        self.assertTrue(jac_norm<eps)

    return test

# Generate tests using function closure
for use_cse in [False, True]:
#for use_cse in [False]:
    for state_repr in state_repr_opts:
    #for state_repr in ["named"]:
        for param_repr in param_repr_opts:
        #for param_repr in ["named"]:
            for body_repr in body_repr_opts:
            #for body_repr in ["named"]:
                for body_optimize in body_optimize_opts:
                #for body_optimize in ["none"]:
                    for float_precision in ["double", "single"]:
                    #for float_precision in ["double"]:
                        
                        test_name = "_".join([body_repr, body_optimize, \
                                              state_repr, param_repr, float_precision])
                        test_name += "_use_cse_" + str(use_cse)

                        test_map["test_"+test_name] = function_closure(\
                            body_repr, body_optimize, param_repr, state_repr, use_cse, \
                        float_precision)

# Generate an empy test class
class TestCodeComponent(unittest.TestCase):pass

for param_name, state_name, body_name in [["PARAMETERS", "STATES", "ALGEBRAIC"], \
                                          ["params", "states", "algebraic"]]:
    for use_cse in [False, True]:
        for body_repr in ["array", "reused_array"]:
            for body_in_arg in [False, True]:

                test_name = "_".join([param_name, body_name, body_repr, \
                                      str(body_in_arg)])
                test_name += "_use_cse_" + str(use_cse)
                setattr(TestCodeComponent, "test_"+test_name, function_closure(\
                    body_repr, "none", "array", "array", use_cse, \
                        "double", param_name, state_name, body_name, body_in_arg))

# Populate the test class with methods
for test_name, test_function in test_map.items():
    setattr(TestCodeComponent, test_name, test_function)
        
if __name__ == "__main__":
    unittest.main()
