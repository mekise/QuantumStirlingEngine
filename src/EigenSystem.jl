function evalΩ1(Δp::Function, ωm, G)
    t -> (1/4)*(ωm^2+Δp(t)^2-sqrt(ωm^4+16*G^2*ωm*Δp(t)-2*ωm^2*Δp(t)^2+Δp(t)^4))
end

function evalΩ2(Δp::Function, ωm, G)
    t -> (1/4)*(ωm^2+Δp(t)^2+sqrt(ωm^4+16*G^2*ωm*Δp(t)-2*ωm^2*Δp(t)^2+Δp(t)^4))
end

function evalP11(Δp::Function, ωm, G)
    t -> -((ωm^2-Δp(t)^2+sqrt(ωm^4+16*G^2*ωm*Δp(t)-2*ωm^2*Δp(t)^2+Δp(t)^4))/(4*G*sqrt(1+(1/16)*abs((ωm^2-Δp(t)^2+sqrt(ωm^4+16*G^2*ωm*Δp(t)-2*ωm^2*Δp(t)^2+Δp(t)^4))/(G*sqrt(ωm*Δp(t))))^2)*sqrt(ωm*Δp(t))))
end

function evalP12(Δp::Function, ωm, G)
    t -> return -((ωm^2-Δp(t)^2-sqrt(ωm^4+16*G^2*ωm*Δp(t)-2*ωm^2*Δp(t)^2+Δp(t)^4))/(4*G*sqrt(1+(1/16)*abs((ωm^2-Δp(t)^2-sqrt(ωm^4+16*G^2*ωm*Δp(t)-2*ωm^2*Δp(t)^2+Δp(t)^4))/(G*sqrt(ωm*Δp(t))))^2)*sqrt(ωm*Δp(t))))
end

function evalP21(Δp::Function, ωm, G)
    t -> 1/sqrt(1+(1/16)*abs((ωm^2-Δp(t)^2+sqrt(ωm^4+16*G^2*ωm*Δp(t)-2*ωm^2*Δp(t)^2+Δp(t)^4))/(G*sqrt(ωm*Δp(t))))^2)
end

function evalP22(Δp::Function, ωm, G)
    t -> 1/sqrt(1+(1/16)*abs((ωm^2-Δp(t)^2-sqrt(ωm^4+16*G^2*ωm*Δp(t)-2*ωm^2*Δp(t)^2+Δp(t)^4))/(G*sqrt(ωm*Δp(t))))^2)
end

function evalP(Δp::Function, ωm, G)
    P11 = evalP11(Δp, ωm, G)
    P12 = evalP12(Δp, ωm, G)
    P21 = evalP21(Δp, ωm, G)
    P22 = evalP22(Δp, ωm, G)
    t -> [P11(t) 0 P12(t) 0
          0 P11(t) 0 P12(t)
          P21(t) 0 P22(t) 0
          0 P21(t) 0 P22(t)]
end