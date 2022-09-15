
function evalΔ(Ω::Function, T, γ)
    t -> (decay(-Ω(t), T, γ) - decay(Ω(t), T, γ))/(4*Ω(t))
end

function evalΣ(Ω::Function, T, γ)
    t -> (decay(-Ω(t), T, γ) + decay(Ω(t), T, γ))/(4*Ω(t))
end

function evalξ(Pa::Function, Pb::Function, Δa::Function, Δb::Function)
    t -> Pa(t)^2*Δa(t) + Pb(t)^2*Δb(t)
end

function evalM(Δp::Function, ωm, G, Th, Tc, γh, γc)
    Ω1 = evalΩ1(Δp, ωm, G)
    Ω2 = evalΩ2(Δp, ωm, G)
    P11 = evalP11(Δp, ωm, G)
    P12 = evalP12(Δp, ωm, G)
    P21 = evalP21(Δp, ωm, G)
    P22 = evalP22(Δp, ωm, G)
    Δh1 = evalΔ(Ω1, Th, γh)
    Δc1 = evalΔ(Ω1, Tc, γc)
    Δh2 = evalΔ(Ω2, Th, γh)
    Δc2 = evalΔ(Ω2, Tc, γc)
    ξ1 = evalξ(P11, P21, Δh1, Δc1)
    ξ2 = evalξ(P12, P22, Δh2, Δc2)
    t -> [2*ξ1(t) 0 2 0 0 0 0 0 0 0
          0 2*ξ1(t) -2*Ω1(t)^2 0 0 0 0 0 0 0
          -Ω1(t)^2 1 2*ξ1(t) 0 0 0 0 0 0 0
          0 0 0 2*ξ2(t) 0 2 0 0 0 0
          0 0 0 0 2*ξ2(t) -2*Ω2(t)^2 0 0 0 0
          0 0 0 -Ω2(t)^2 1 2*ξ2(t) 0 0 0 0
          0 0 0 0 0 0 ξ1(t)+ξ2(t) 0 -Ω2(t)^2 1
          0 0 0 0 0 0 0 ξ1(t)+ξ2(t) -Ω1(t)^2 1
          0 0 0 0 0 0 1 1 ξ1(t)+ξ2(t) 0
          0 0 0 0 0 0 -Ω1(t)^2 -Ω2(t)^2 0 ξ1(t)+ξ2(t)]
end

function evalC(Δp::Function, ωm, G, Th, Tc, γh, γc)
    Ω1 = evalΩ1(Δp, ωm, G)
    Ω2 = evalΩ2(Δp, ωm, G)
    P11 = evalP11(Δp, ωm, G)
    P12 = evalP12(Δp, ωm, G)
    P21 = evalP21(Δp, ωm, G)
    P22 = evalP22(Δp, ωm, G)
    Σh1 = evalΣ(Ω1, Th, γh)
    Σc1 = evalΣ(Ω1, Tc, γc)
    Σh2 = evalΣ(Ω2, Th, γh)
    Σc2 = evalΣ(Ω2, Tc, γc)
    t -> [(P11(t)^2*Σh1(t))/Ω1(t)+(P21(t)^2*Σc1(t))/Ω1(t),
           P11(t)^2*Ω1(t)*Σh1(t)+P21(t)^2*Ω1(t)*Σc1(t),
           0,
           (P12(t)^2*Σh2(t))/Ω2(t)+(P22(t)^2*Σc2(t))/Ω2(t),
           P12(t)^2*Ω2(t)*Σh2(t)+P22(t)^2*Ω2(t)*Σc2(t),
           0, 0, 0, 0, 0]
end