# Packages
using Robin2011_RepPackage
using Optimization, OptimizationNLopt
using CSV
using DataFrames
using Distributions, Statistics, LinearAlgebra
using DelimitedFiles
using Random

# Parameters and Grids
M = 50
N = 50


b = params_default()

g = grids(b; M = M, N = N)

x = g[:x]
y = g[:y]
l = g[:l]
Π = g[:Π]
F = g[:F]

# Draw
T = 5000
burn = 1000
Random.seed!(2077)
draw = rand(T+burn)

# Data
# dta = CSV.read("/Users/bojs/Desktop/Robin 2011 Rep Files/matlab/USquarterly.csv", DataFrame, header = false)

# filter!(row -> !isnan(row.Column1) && !isnan(row.Column4), dta);

# prod = collect(dta.Column1)
# prod = prod/mean(prod)

# Val Added and Home Production
p = matchprod(x, y; B = b.B, C = b.C)
z = homeprod(x, y; B = b.B, C = b.C, α = b.α, z0 = b.z0)


# Callback Function
function callback(state, loss_val)
    display(loss_val)
    return false
end

# Baseline - No Policy
baseline = optCrit([0,0], zeros(N, M), 0.0, 0.0; T = T, burn = burn, M = M, N = N, b = b, grid = g, fullinfo = true)

# Min Wage
fmw = OptimizationFunction((u, _) -> optCrit([0,0], zeros(N, M), u[1], 0.0;T = T, burn = burn, draw = draw, M = M, N = N, b = b, grid = g, fullinfo = false))

probmw = OptimizationProblem(fmw, [0.75])

sol = solve(probmw, NLopt.LN_SBPLX(); maxiters = 100000, reltol = 1e-8)

mw = 0.75

writedlm("output/wminwage.txt", sol.u)

# opt UI

fui = OptimizationFunction((u, minwage) -> optCrit([0,0], uiyx(u[1], p), minwage[1], u[2];T = T, burn = burn, draw = draw, M = M, N = N, b =b, grid = g, fullinfo = false))

#  No MW
probui = OptimizationProblem(fui, [0.0, 0.0], [0.0]; lb = [0.0, 0.0], ub = [1.0, 1.0])

sol = solve(probui, NLopt.LN_SBPLX(); reltol = 1e-8, callback =callback)

ui = sol.u

writedlm("output/ui.txt", ui)

# MW
probui_mw = OptimizationProblem(fui, [0.0,0.0], mw; lb = [0.0, 0.0], ub = [1.0, 1.0])

sol = solve(probui_mw, NLopt.LN_SBPLX(); reltol = 1e-4, callback =callback)

ui_mw = sol.u

writedlm("output/ui_mw.txt", ui_mw)

# opt Subsidy

fsub = OptimizationFunction((u, minwage) -> optCrit([u[1],u[2]], zeros(N, M), minwage[1], u[3];T = T, burn = burn, draw = draw, M = M, N = N, b =b, grid = g, fullinfo = false))

# No MW
probsub = OptimizationProblem(fsub, [1.0, 0.0, 0.0], 0; lb = [minimum(y), 0.0, 0.0], ub = [maximum(y), 10.0, 1.0])

sol = solve(probsub, NLopt.LN_SBPLX(); reltol = 1e-4, callback =callback)

sub = sol.u

writedlm("output/sub.txt", sub)

# MW
probsub_mw = OptimizationProblem(fsub, [1.0, 0.0, 0.0], mw; lb = [minimum(y), 0.0, 0.0], ub = [maximum(y), 1.0, 1.0])

sol = solve(probsub_mw, NLopt.LN_SBPLX(); reltol = 1e-8, callback =callback)

sub_mw = sol.u

writedlm("output/sub_mw.txt", sub_mw)

# Both

fboth = OptimizationFunction((u, minwage) -> optCrit([u[1], u[2]],uiyx(u[3], p), minwage[1], u[4];T = T, burn = burn, draw = draw, M = M, N = N, b = b, grid = g, fullinfo = false))

# No MW
probboth = OptimizationProblem(fboth, [1.0, 0.0, ui[1], ui[2]], 0.0; lb = [minimum(y), 0.0, 0.0, 0.0], ub = [maximum(y), 10.0, 1.0, 1.0])

sol = solve(probboth, NLopt.LN_SBPLX(); reltol = 1e-8, callback =callback)

both = sol.u

writedlm("output/both.txt", both)

probboth_mw = OptimizationProblem(fboth, [1.0, 0.0, ui_mw[1], ui_mw[2]], mw; lb = [minimum(y), 0.0, 0.0, 0.0], ub = [maximum(y), 10.0, 1.0, 1.0])

sol = solve(probboth_mw, NLopt.LN_SBPLX(); reltol = 1e-8, callback =callback)

both_mw = sol.u

writedlm("output/both_mw.txt", both_mw)
