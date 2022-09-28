module QStirlingEngine

using LinearAlgebra
using Random
using DifferentialEquations

include("Decay.jl");
include("EigenSystem.jl");
include("Matrices.jl")
include("EqSolver.jl")

export n, decay
export evalΩ1, evalΩ2, evalP11, evalP12, evalP21, evalP22, evalP
export evalΔ, evalΣ, evalξ, evalM, evalC
export eqSolver

end