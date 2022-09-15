function solver(tspan, ρ0, Δp::Function, ωm, G, Th, Tc, γh, γc; saveat=[])
    P = evalP(Δp, ωm, G)
    M = evalM(Δp, ωm, G, Th, Tc, γh, γc)
    C = evalC(Δp, ωm, G, Th, Tc, γh, γc)
    u0matrix = transpose(P(0))*ρ0*P(0)
    u0 = [u0matrix[1,1], u0matrix[2,2], u0matrix[1,2], u0matrix[3,3], u0matrix[4,4], u0matrix[3,4], u0matrix[1,4], u0matrix[2,3], u0matrix[1,3], u0matrix[2,4]]
    function f(du, u, par, t)
        du .= M(t)*u+C(t)
    end
    prob = ODEProblem(f, u0, tspan)
    return solve(prob, Rosenbrock23(), abstol=1e-8, reltol=1e-8, maxiters=Int(1e7), saveat=saveat)
end