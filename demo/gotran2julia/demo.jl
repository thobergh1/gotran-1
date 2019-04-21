# Checkout http://docs.juliadiffeq.org/latest/
using DifferentialEquations

# Include file generated by gotran
include("tentusscher_panfilov_2006_M_cell.jl")

# Initial states
y0 = init_state_values()

# Parameters
parameters = init_parameter_values()

# Time steps
tspan = (0.0,100.0)
prob = ODEProblem(rhs,y0,tspan, parameters)
sol = solve(prob)

V_idx = state_indices("V")
using Plots
p1 = plot(sol,vars=(V_idx), title="State V")

# Monitored values
monitored = zeros((90, length(sol.t)))
for i = 1:length(sol.t)
    monitored[:, i] = monitor(monitored[:, i], sol.u[i], parameters, sol.t[i])
end

i_Kr_idx= monitor_indices("i_Kr")
i_Kr = monitored[i_Kr_idx, :]
p2 = plot(sol.t, i_Kr, title="iKr")

plot(p1, p2)
png("results_julia")

#plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
#     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
#plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")
