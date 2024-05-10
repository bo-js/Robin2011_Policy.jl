using Robin2011_RepPackage, Roots, Random, StatsBase, DelimitedFiles

T = 5000
burn = 1000

N = 50
M = 50

b = params_default()
g = grids(b; M = M, N = N)

p = matchprod(g[:x], g[:y]; B = b.B, C = b.C)

Random.seed!(2077)
draw = rand(T+burn)

ui_range = LinRange(0.001, 0.04, 40)

function find_tax_ui(ui)
    ui_matrix = uiyx(ui, p)
    tax = find_zero(t -> deficit(t, ui_matrix, [0,0]; M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g), 0)
    return tax
end


ui_taxes = find_tax_ui.(ui_range)

writedlm("output/ui_grid.txt", ui_range)
writedlm("output/ui_taxes.txt", ui_taxes)

Π = g[:Π]

statet = states(draw, Π)
statet = statet[burn+1:burn+T]

statesort = sort(statet)

statelow = maximum(statesort[1:2000])
statehigh = maximum(statesort[2000:4000])

function find_tax_sub(sub)
    tax = find_zero(t -> deficit(t, zeros(N, M), [g[:y][statelow], sub]; M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g), 0)
    return tax
end

sub_range = LinRange(0.001, 0.1, 100)
writedlm("output/sub_grid.txt", sub_range)
sub_taxes = find_tax_sub.(sub_range)
writedlm("output/sub_taxes.txt", sub_taxes)