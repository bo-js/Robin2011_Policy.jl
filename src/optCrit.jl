function optCrit(zsub, uizx, wmin, tax; M = 500, N =100,T = 5000, burn = 1000, draw = rand(burn+T), b = params_default(), grid = grids(b; M = M, N = N))
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

    Stup = SminVFI(wmin, sub + (1-tax) .* p, z + uizx, Π; β = β, λ0 = λ0, λ1 = λ1, r = r)
    Sx = Stup.S
    Smin = Stup.Smin

    Ux = (I(N) - Π./(1 + r))\(z + uizx + (λ1/(1+r)) .* Π * ((Sx .> max.(Smin, 0)) .* max.(Smin, 0)))

    ux = (δ/(δ + λ0)).*(Sx .> max.(Smin, 0)) .+ (Sx .≤ max.(Smin,0))

    statet = zeros(Int64, T+burn)
    i = Integer(round(N/2))
    statet[1] = i
    yt = y[i] * ones(T+burn)
    uxt = repeat(ux[i, :]', T+burn+1, 1)
    Sxt = repeat(Sx[i, :]', T+burn+1, 1)

    for t = 1:T+burn
        i = min(1 + sum(draw[t] .> cumsum(Π[i, :])), N)
        statet[t] = i
        yt[t] = y[i]
        Sxt[t, :] = Sx[i, :]
        uxt[t+1, :] = 1 .- (Sxt[t, :] .> max.(Smin[i, :],0)) .* ((1 - δ) .* (1 .- uxt[t, :]) + λ0 .* uxt[t, :])
    end
    ut = uxt * l

    wd = WageVFI(Sx, Smin, Π, z; λ1 = λ1, β = β)

    wdt = wage_dens_path(Sx, Smin, uxt, wd, l, Ux, statet, T+burn; λ0 = λ0, λ1 = λ1, δ = δ)

    swf = sum(sum(wdt[t, :, :, 1] .* wd[:wmin]) + sum(wdt[t, :, :, 2] .* wd[:wmin]) + dot(z[statet[t], :] + uizx[statet[t], :], uxt[t, :] .* l) for t in (burn+1):(T+burn))/T

    deficit = sum(sub[statet[t], 1] * (1 - ut[t]) + dot(uizx[statet[t], :], uxt[t, :] .* l) - tax * dot(p[statet[t], :], (1 .- uxt[t, :]) .* l) for t in (burn+1):(T+burn))/T

    if deficit > 0.1
        return -Inf
    else
        return swf
    end
end

function optCrit(zsub, uizx, wmin, tax, prod; M = 500, N =100, b = params_default(), grid = grids(b; M = M, N = N))
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

    T = length(prod)
    burn = 0
    T1 = burn+T;


    x = grid[:x]
    y = grid[:y]
    Π = grid[:Π]
    l = grid[:l]

    p = matchprod(x, y; B = B, C = C)
    z = homeprod(x, y; B = B, C = C, α = α, z0 = z0)

    sub = [(y[i] < zsub[1])*zsub[2] for i in 1:N, j in 1:M]

    Stup = SminVFI(wmin, sub + (1-tax) .* p, z + uizx, Π; β = β, λ0 = λ0, λ1 = λ1, r = r)
    Sx = Stup.S
    Smin = Stup.Smin

    Ux = (I(N) - Π./(1 + r))\(z + uizx + (λ1/(1+r)) .* Π * ((Sx .> max.(Smin, 0)) .* max.(Smin, 0)))

    ux = (δ/(δ + λ0)).*(Sx .> max.(Smin, 0)) .+ (Sx .≤ max.(Smin,0))

    uxt = ones(T1+1, M)
    yt = ones(T1)
    statet = zeros(Int, T1)
    
    i = min(sum(y.-prod[1] .<= 0), N);
    yt[1] = y[i];
    statet[1] = i;
    uxt[1, :] = ux[i,:]

    ## Productivity and Unemployment Dynamics
    for t in 1:(burn+T)

        uxt_1  = 1 .- (Sx .> max.(Smin, 0)) .* repeat((1 - δ) * (1 .- uxt[t, :]) + λ0 * uxt[t, :], 1, N)'
        ext = repeat(l', N, 1) .* (1 .- uxt_1)
        et = sum(ext, dims = 2)
        xt = vec(sum(p .* ext, dims = 2)./et)
            
        global i = argmin(abs.(xt .- prod[t]))
        statet[t] = i
        yt[t] = y[i]

        uxt[t+1, :] = [1 - (Sx[i, m] > max(Smin[i,m], 0)) * ((1 - δ) * (1 - uxt[t, m]) + λ0 * uxt[t, m]) for m in 1:M]
    end
    ut = uxt * l

    wd = WageVFI(Sx, Smin, Π, z; λ1 = λ1, β = β)

    wdt = wage_dens_path(Sx, Smin, uxt, wd, l, Ux, statet, T+burn; λ0 = λ0, λ1 = λ1, δ = δ)

    swf = sum(sum(wdt[t, :, :, 1] .* wd[:wmin]) + sum(wdt[t, :, :, 2] .* wd[:wmin]) + dot(z[statet[t], :] + uizx[statet[t], :], uxt[t, :] .* l) for t in (burn+1):(T+burn))/T

    deficit = sum(sub[statet[t], 1] * (1 - ut[t]) + dot(uizx[statet[t], :], uxt[t, :] .* l) - tax * dot(p[statet[t], :], (1 .- uxt[t, :]) .* l) for t in (burn+1):(T+burn))/T

    if deficit > 0.1
        return -Inf
    else
        return swf
    end
end