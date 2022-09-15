n = (ω, T) -> (exp(ω/T)-1)^(-1)

function decay(ω, T, γ; Ωc=1000, d=1)
    return 2*γ*ω^d*(Ωc^2/(ω^2+Ωc^2))*(1+n(ω, T))
end