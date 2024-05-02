module Robin2011_RepPackage

using LinearAlgebra
using Distributions
using Base.Threads

function matchprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686)
    return y * (B * x .+ C)'
end

function homeprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686, α::Number = 0.75, z0::Number = 0.766752794650811)
    return z0 .+ α * (matchprod(x, y; B, C) .- z0)
end

function SurplusVFI(p::Matrix, z::Matrix, Π::Matrix; β::Number = 0.946603693905558)
    S = (I(length(Π[1,:])) - β * Π )\(p - z)
    e = norm(S - max.(S, 0), 2)

    while e > 0.00001
        S1 = S
        S = p - z + β * Π * max.(S, 0)
        e = norm(S - S1, 2)
    end

    return S
end

include("wages.jl")

include("dynamics.jl")

include("grids.jl")

include("estCrit.jl")

include("min_surplus.jl")

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
