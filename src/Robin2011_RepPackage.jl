module Robin2011_RepPackage

using CSV
using DataFrames
using LinearAlgebra
using Distributions
using Base.Threads

function matchprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686)
    return y * (B * x .+ C)'
end

function homeprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686, α::Number = 0.64, z0::Number = 0.766752794650811)
    return z0 .+ α * (matchprod(x, y; B, C) .- z0)
end

function SurplusVFI(p::Matrix, z::Matrix, Π::Matrix;tol = 0.00001, β::Number = 0.946603693905558)
    S = (I(length(Π[1,:])) - β * Π )\(p - z)
    e = norm(S - max.(S, 0), 2)

    while e > tol
        S1 = S
        S = p - z + β * Π * max.(S, 0)
        e = norm(S - S1, 2)
    end

    return S
end

function states(draw, Π)
    T1 = length(draw)
    statet = zeros(Int64, T1)
    i = Integer(round(length(Π[1, :])/2))
    statet[1] = i

    for t = 1:T1
        i = min(1 + sum(draw[t] .> cumsum(Π[i, :])), length(Π[1, :]))
        statet[t] = i
    end

    return statet
end

include("deficit.jl")

include("wages.jl")

include("dynamics.jl")

include("grids.jl")

include("estCrit.jl")

include("min_surplus.jl")

include("optCrit.jl")

include("params_default.jl")

include("uiyx.jl")

export states
export uiyx
export params_default
export optCrit
export SminVFI
export estCrit
export grids
export wage_dens_path
export SurplusVFI
export WageVFI
export matchprod
export homeprod
export unemp_path

end
