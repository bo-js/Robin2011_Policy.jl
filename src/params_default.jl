function params_default(; path = "x1.txt")
    b = parse.(Float64, readlines(path))
    
    logit = x -> 1/(1 + exp(x))

    ν = exp(b[1])
    μ = exp(b[2])
    δ = logit(b[3])
    λ0 = logit(b[4])
    ρ = tanh(b[5])
    σ = exp(b[6])
    z0 = exp(b[7])
    C = exp(b[8])

    r = 0.05/4
    τ = 0.5
    β = (1-δ)/(1+r)
    α = 0.64
    λ1 = 0.12 * λ0
    B = 1

    return (
        r = r,
        δ = δ,
        λ0 = λ0,
        λ1 = λ1,
        ρ = ρ,
        σ = σ,
        ν = ν,
        μ = μ,
        τ = τ,
        β = β,
        α = α,
        z0 = z0,
        C = C,
        B = B
    )
end
