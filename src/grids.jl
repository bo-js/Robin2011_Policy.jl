using Distributions
using Copulas

function grids(; M = 500, N = 100, ρ = 0.913702234286476, σ = 0.0257,ν::Number = 2.019365636076711, μ::Number = 5.786082109731152 )
     trim = 0.002

     x = collect(LinRange(trim, 1 - trim, M))
     F = collect(LinRange(trim, 1- trim, N))
     y = exp.(σ * quantile(Normal(), F))
     l = pdf(Beta(ν, μ), x)
     l = l ./ sum(l)

     Σ = [1 ρ
          ρ 1]

     cop = GaussianCopula(Σ)

     P = [pdf(cop, [F[i], F[j]]) for i in 1:N, j in 1:N ]

     Π = [P[i, j]/sum(P[i, :]) for i in 1:N, j in 1:N]

     return Dict(:x => x, :y => y, :F => F, :l => l,:Π => Π)
end

function grids(params::NamedTuple; M = 500, N = 100, trim = 0.002)
     
     ρ = params.ρ
     σ = params.σ
     ν = params.ν
     μ = params.μ

     x = collect(LinRange(trim, 1 - trim, M))
     F = collect(LinRange(trim, 1- trim, N))
     y = exp.(σ * quantile(Normal(), F))
     l = pdf(Beta(ν, μ), x)
     l = l ./ sum(l)

     Σ = [1 ρ
          ρ 1]

     cop = GaussianCopula(Σ)

     P = [pdf(cop, [F[i], F[j]]) for i in 1:N, j in 1:N ]

     Π = [P[i, j]/sum(P[i, :]) for i in 1:N, j in 1:N]

     return Dict(:x => x, :y => y, :F => F, :l => l,:Π => Π)
end
