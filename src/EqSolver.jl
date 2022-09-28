function eqSolver(u0, tspan, Δp::Function, ωm, G, Th, Tc, γh, γc)
    M = evalM(Δp, ωm, G, Th, Tc, γh, γc)
    C = evalC(Δp, ωm, G, Th, Tc, γh, γc)
    function f(du, u, par, t)
        du .= M(t)*u + C(t)
    end
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Rosenbrock23(), abstol=1e-8, reltol=1e-8, maxiters=Int(1e7), progress=true, progress_steps=20000)
    return sol
end