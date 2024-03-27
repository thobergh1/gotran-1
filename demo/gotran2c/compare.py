import os
from ctypes import c_double
from ctypes import c_int
from ctypes import cdll


import matplotlib.pyplot as plt
import numpy as np
# import seaborn as sns

import time

path = "~/Dokumenter/gotran/venv/gotran/demo/gotran2c"
expanded_path = os.path.expanduser(path)
os.chdir(expanded_path)

libname = "libdemo.so"
libdir = "."
libpath = os.path.join(libdir, libname)


assert os.path.isfile(libpath)

"""
nm -gC libdemo.so
"""



#libbase_model = cdll[libpath]
libbase_model = np.ctypeslib.load_library(libname, libdir)


# Get number of states and parameters from the C library
num_states = libbase_model.state_count()
num_parameters = libbase_model.parameter_count()

#print(num_states)
#print(num_parameters)


def init_lib():
    """
    Make sure that arrays passed to C is of the correct types.
    """

    float64_array = np.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="contiguous")
    float64_array_2d = np.ctypeslib.ndpointer(
        dtype=c_double,
        ndim=2,
        flags="contiguous",
    )

    libbase_model.init_state_values.restype = None  # void
    libbase_model.init_state_values.argtypes = [float64_array]

    libbase_model.init_parameters_values.restype = None  # void
    libbase_model.init_parameters_values.argtypes = [float64_array]

    solve_functions = [
        libbase_model.ode_solve_forward_euler,
        libbase_model.ode_solve_rush_larsen,
        libbase_model.ode_solve_generalized_rush_larsen
    ]

    for func in solve_functions:
        func.restype = None  # void
        func.argtypes = [  
            float64_array,  # u
            float64_array,  # parameters
            float64_array_2d,  # u_values
            float64_array,  # t_values
            c_int,  # num_timesteps
            c_int,  # num_cells
            c_double,  # v_step
            c_double,  # dt
        ]


def init_parameters():
    parameters = np.zeros(num_parameters, dtype=np.float64)
    libbase_model.init_parameters_values(parameters)
    return parameters


def solve(t_start, t_end, dt, V_step, num_steps=None, method="fe"):
    parameters = init_parameters()
    
    #print("Here", type(dt))

    
    if type(dt) is not float:
        dt = float(dt)
    if num_steps is not None:
        assert type(num_steps) is int
        t_end = dt * num_steps
    else:
        num_steps = round((t_end - t_start) / dt)

    t_values = np.linspace(t_start, t_end, num_steps + 1)

    num_cells = 1
    padded_num_cells = 8
    
    u = np.zeros(num_states*padded_num_cells, dtype=np.float64)
    libbase_model.init_state_values(u, num_cells, padded_num_cells)

    # print(u)
    # print(u.shape[0])
    

    u_values = np.zeros((num_steps + 1, int(u.shape[0]/padded_num_cells)), dtype=np.float64)
    
    for i in range(int(u.shape[0]/padded_num_cells)):
        u_values[0, i] = u[i*padded_num_cells]


    # print(u_values)
    # print(u_values.shape[0])
    # print(u_values.shape)


    if method == "fe":
        libbase_model.ode_solve_forward_euler(
            u,
            parameters,
            u_values,
            t_values,
            num_steps,
            num_cells,
            V_step,
            dt,
        )
    elif method == "rl":
        libbase_model.ode_solve_rush_larsen(
            u,
            parameters,
            u_values,
            t_values,
            num_steps,
            num_cells,
            V_step,
            dt,
        )
    elif method == "grl1":
        libbase_model.ode_solve_rush_larsen(
            u,
            parameters,
            u_values,
            t_values,
            num_steps,
            num_cells,
            V_step,
            dt,
        )
    else:
        raise ValueError(f"Invalid method {method}")

    return t_values, u_values

def main():
    """   
    t_start = 0.0
    t_end = 200.0
    dt = 1e-4

    total = 0.0

    t_values, u_values = solve(t_start, t_end, dt, method="fe")
    #t_values, u_values = solve(t_start, t_end, dt, method="rush_larsen")

    V_idx = libbase_model.state_index("V")

    fig, ax = plt.subplots()
    ax.plot(t_values[:], u_values[:, V_idx])
    ax.set_title("Membrane potential")
    ax.set_xlabel("Time (ms)")
    plt.show()
    """ 


    t_start = 0.0
    t_end = 100.0
    dt = 1e-5
    
    # V_step = 10
    
    # V_steps = [10,1,1e-1,1e-2,1e-3,1e-4,1e-5]
    V_steps = [10, 5 , 1 , 0.1, 0.5 , 0.01, 0.05, 0.001, 0.005, 0.0001, 0.0005, 0.00001]



    V_idx = libbase_model.state_index("V")

    infilename = f"original_tentusscher_dt_{dt}.txt"
    infilepath = "/home/thomas/Dokumenter/gotran/venv/gotran/demo/gotran2c/results/rmse/TP06/"+infilename


    outpath = "results/rmse/TP06"
    outfilename = f"TP06_dt_{dt}_changing_V_steps.txt"
    outfilepath = os.path.join(outpath, outfilename)
    os.makedirs(os.path.dirname(outfilepath), exist_ok=True)

    # Open the file in read mode
    with open(infilepath, 'r') as file:
        lines = file.readlines()

    original_v = []

    for line in lines:
        original_v.append(float(line))

    with open(outfilepath, 'w') as file:
        methods = ["fe", "rl", "grl1"]
        for mtd in methods:
            file.write(f"Method: {mtd}\n")
            state_values = []
            rmse_values = []

            for V_step in V_steps:
                t_values, u_values = solve(t_start, t_end, dt, V_step, method=mtd)
                state_values.append(u_values[:, V_idx])
                
            for i, state_v in enumerate(state_values):
                mse = np.mean((original_v - state_v)**2, axis=0)
                rmse = np.sqrt(mse)
                rmse_values.append(rmse)
                print(f"RMSE: {rmse:10.8f} | Step size {V_steps[i]:10.6f}")
                file.write(f"RMSE: {rmse:10.8f} | Step size {V_steps[i]:10.6f}\n")
            file.write("\n\n")

    file.close()

if __name__ == "__main__":
    init_lib()
    main()
