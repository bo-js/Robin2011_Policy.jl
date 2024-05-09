using Robin2011_RepPackage, NLsolve, Roots, Random, LinearAlgebra

function deficit(tax, uizx, zsub; wmin = 0, M = 50, N =50,T = 6000, burn = 1000, draw = rand(burn+T), b = params_default(), grid = grids(b; M = M, N = N), tol = 0.001)
    r = b.r
    δ = b.δ
    λ0 = b.λ0
    λ1 = b.λ1
    ρ = b.ρ
    σ = b.σ
    ν = b.ν
    μ = b.μ
    τ = b.τ
    α = b.α
    β = b.β
    z0 = b.z0
    C = b.C
    B = b.B

    x = grid[:x]
    y = grid[:y]
    Π = grid[:Π]
    l = grid[:l]

    p = matchprod(x, y; B = B, C = C)
    z = homeprod(x, y; B = B, C = C, α = α, z0 = z0)

    sub = [(y[i] < zsub[1])*zsub[2] for i in 1:N, j in 1:M]

    Stup = SminVFI(wmin, sub + (1-tax) .* p, z + uizx, Π; β = β, λ0 = λ0, λ1 = λ1, r = r, tol = tol)
    if Stup.success == false
        return Inf
    end
    
    Sx = Stup.S
    Smin = Stup.Smin

    ux = (δ/(δ + λ0)).*(Sx .> max.(Smin, 0)) .+ (Sx .≤ max.(Smin,0))

    statet = zeros(Int64, T+burn)
    i = Integer(round(N/2))
    statet[1] = i
    yt = y[i] * ones(T+burn)
    uxt = ones(T+burn+1, M)
    uxt[1, :] = ux[i, :]
    Sxt = repeat(Sx[i, :]', T+burn+1, 1)

    for t = 1:T+burn
        i = min(1 + sum(draw[t] .> cumsum(Π[i, :])), N)
        statet[t] = i
        yt[t] = y[i]
        Sxt[t, :] = Sx[i, :]
        uxt[t+1, :] = 1 .- (Sxt[t, :] .> max.(Smin[i, :],0)) .* ((1 - δ) .* (1 .- uxt[t, :]) + λ0 .* uxt[t, :])
    end

    deficit = sum(dot(uizx[statet[t], :], uxt[t, :] .* l)  + dot(sub[statet[t], :] .- tax .* p[statet[t], :], (1 .- uxt[t, :]) .* l) for t in (burn+1):(T+burn))/T

    return deficit
end

export deficit