using Robin2011_RepPackage, Roots, Random, StatsBase, DelimitedFiles

T = 5000
burn = 1000

N = 100
M = 50

b = params_default()
g = grids(b; M = M, N = N)

p = matchprod(g[:x], g[:y]; B = b.B, C = b.C)

Random.seed!(2077)
draw = rand(T+burn)

ui_range = LinRange(0.001, 0.04, 40)

function find_tax_ui(ui; wmin = 0)
    ui_matrix = uiyx(ui, p)
    tax = find_zero(t -> deficit(t, ui_matrix, [0,0]; wmin = wmin, M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g), 0)
    return tax
end


ui_taxes = find_tax_ui.(ui_range)

ui_taxes_mw = find_tax_ui.(ui_range; wmin = 0.76)

writedlm("output/ui_grid.txt", ui_range)
writedlm("output/ui_taxes.txt", ui_taxes)
writedlm("output/ui_taxes_mw.txt", ui_taxes_mw)

Π = g[:Π]

statet = states(draw, Π)
statet = statet[burn+1:burn+T]

statesort = sort(statet)

statelow = maximum(statesort[1:2000])
statehigh = maximum(statesort[2000:4000])

function find_tax_sub(sub; wmin = 0, tol = 0.0001)
    tax = find_zero(t -> deficit(t, zeros(N, M), [g[:y][statelow], sub];tol = tol, wmin = wmin, M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g), 0)
    return tax
end

sub_range = LinRange(0.0005, 0.05, 100)
writedlm("output/sub_grid.txt", sub_range)
sub_taxes = find_tax_sub.(sub_range)
writedlm("output/sub_taxes.txt", sub_taxes)
sub_range_mw = LinRange(0.0005, 0.01, 20)
writedlm("output/sub_grid_mw.txt", sub_range_mw)
sub_taxes_mw = find_tax_sub.(sub_range_mw; wmin = 0.76, tol = 0.01)
writedlm("output/sub_taxes_mw.txt", sub_taxes_mw)