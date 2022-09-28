using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

include("../src/QStirlingEngine.jl")
using .QStirlingEngine
using LinearAlgebra
using Random
using DifferentialEquations
using NPZ

tend = 100_000.
tspan = (0., tend)
rΔp, ωm, G, Th, Tc, γh, γc = rand(Float64, (7))
while Th<Tc
    global Th, Tc = rand(Float64, (2))
end
γh *= 10^(-3)
γc *= 10^(-3)
ρ0 = [4.08869024190911 0.3020135004112786 -2.6885011832628827 -0.4528955715802489
      0.3020135004112787 3.0426737415267766 0.21373970846112722 4.438387229442656
      -2.6885011832628827 0.21373970846112744 1.8910786497021073 0.9534619295961205
      -0.45289557158024896 4.438387229442656 0.9534619295961203 6.981967317786721]

##################################
# Copy parameters from specific run
##################################
vars = npzread("./data/run_03.npz")
rΔp, ωm, G, Th, Tc, γh, γc, ρ0 = vars["rDeltap"], vars["omegam"], vars["G"], vars["Th"], vars["Tc"], vars["gammah"], vars["gammac"], vars["rho0"]
##################################

Δp = t -> rΔp + 10^(-4)*rΔp*t
# Δp = t -> 2*ωm - (2*ωm - 10^(-2)*ωm)/tspan[end] * t

P = evalP(Δp, ωm, G)
u0matrix = transpose(P(0))*ρ0*P(0)
u0 = [u0matrix[1,1], u0matrix[2,2], u0matrix[1,2], u0matrix[3,3], u0matrix[4,4], u0matrix[3,4], u0matrix[1,4], u0matrix[2,3], u0matrix[1,3], u0matrix[2,4]]

sol = eqSolver(u0, tspan, Δp, ωm, G, Th, Tc, γh, γc)

cov = zeros(Float64, (length(sol.t), 10))
for j in 1:10
    cov[:, j] = [sol.u[i][j] for i in eachindex(sol.t)]
end

npzwrite("./data/run_05.npz", Dict("t" => sol.t,
                                   "cov" => cov,
                                   "rDeltap" => rΔp,
                                   "omegam" => ωm,
                                   "G" => G,
                                   "Th" => Th,
                                   "Tc" => Tc,
                                   "gammah" => γh,
                                   "gammac" => γc,
                                   "rho0" => ρ0,
                                   "tend" => tend))