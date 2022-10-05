using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

include("../src/QStirlingEngine.jl")
using .QStirlingEngine
using LinearAlgebra
using Random
using DifferentialEquations
using NPZ

# Copy parameters from specific run
vars = npzread("./data/run_03_slow_driving.npz")
rΔp, ωm, G, Th, Tc, γh, γc, ρ0 = vars["rDeltap"], vars["omegam"], 0.001*vars["G"], vars["Th"], vars["Tc"], vars["gammah"], vars["gammac"], vars["rho0"]

# Generate data to plot eigenvalues(Δp(t))
tilt = 2 * 10^(-1)
Δp = t -> 2 * ωm * (1-tilt) * t/tend + tilt*ωm

tend = 100_000.
tspan = LinRange(0., tend, 1000)
Δpspan = Δp.(tspan)

eigplus = eigScan.(true, Δpspan, ωm, G) 
eigminus = eigScan.(false, Δpspan, ωm, G)

npzwrite("./data/eig_scan.npz", Dict("tspan" => tspan,
                                     "omegam" => ωm,
                                     "G" => G,
                                     "deltap" => Δpspan,
                                     "eigplus" => eigplus,
                                     "eigminus" => eigminus))