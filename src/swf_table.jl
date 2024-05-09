using Robin2011_RepPackage, Random, LinearAlgebra, Latexify, PrettyTables

T = 5000
burn = 1000

N = 50
M = 50

Random.seed!(2077)
draw = rand(T+burn)

b = params_default()
grid = grids(; N = 50, M = 50)

p = matchprod(grid[:x], grid[:y]; B = params.B, C = params.C)

minwage = 0.75

# Baseline
bsline = optCrit(zeros(2), zeros(N, M), 0, 0; M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = grid, fullinfo = true)

swf_emp_bsline = sum(sum(bsline.wdt[t, :, :, 1] .* bsline.wd[:wmin]) + sum(bsline.wdt[t, :, :, 2] .* bsline.wd[:wmin]) for t in (burn+1):(T+burn))/T
swf_unemp_bsline = sum(dot(z[bsline.statet[t], :], bsline.uxt[t, :] .* l) for t in (burn+1):(burn+T))/T

# Baseline - MW

bsline_mw = optCrit(zeros(2), zeros(N, M), minwage, 0; M = M, N = N, T = T, burn = burn, draw = draw, b = params, grid = grid, fullinfo = true)

swf_emp_bslinemw = sum(sum(bsline_mw.wdt[t, :, :, 1] .* bsline_mw.wd[:wmin]) + sum(bsline_mw.wdt[t, :, :, 2] .* bsline_mw.wd[:wmin]) for t in (draw+1):(T+draw))/T
swf_unemp_bslinemw = sum(dot(z[bsline_mw.statet[t], :], bsline_mw.uxt[t, :] .* l) for t in (burn+1):(burn+T))/T


# UI - No MW
ui = parse.(Float64, readlines("output/ui.txt"))

ui_sim = optCrit(zeros(2), uiyx(ui[1], p), 0, ui[2]; M = M, N = N, T = T, burn = burn, draw = draw, b = params, grid = grid, fullinfo = true)
swf_emp_ui = sum(sum(ui_sim.wdt[t, :, :, 1] .* ui_sim.wd[:wmin]) + sum(ui_sim.wdt[t, :, :, 2] .* ui_sim.wd[:wmin]) for t in (draw+1):(T+draw))/T
swf_unemp_ui = sum(dot(z[ui_sim.statet[t], :] .+ uiyx(ui[1], p)[ui_sim.statet[t], :], ui_sim.uxt[t, :] .* l) for t in (burn+1):(burn+T))/T

# UI - MW

uimw = parse.(Float64, readlines("output/ui_mw.txt"))

uimw_sim = optCrit(zeros(2), uiyx(uimw[1], p), minwage, uimw[2]; M = M, N = N, T = T, burn = burn, draw = draw, b = params, grid = grid, fullinfo = true)
swf_emp_uimw = sum(sum(uimw_sim.wdt[t, :, :, 1] .* uimw_sim.wd[:wmin]) + sum(uimw_sim.wdt[t, :, :, 2] .* uimw_sim.wd[:wmin]) for t in (draw+1):(T+draw))/T
swf_unemp_uimw = sum(dot(z[uimw_sim.statet[t], :] .+ uiyx(uimw[1], p), uimw_sim.uxt[t, :] .* l) for t in (burn+1):(burn+T))/T

# Sub - No MW
sub = parse.(Float64, readlines("output/sub.txt"))

sub_sim = optCrit([sub[1], sub[2]], zeros(N, M), 0.0, sub[3]; M = M, N = N, T = T, burn = burn, draw = draw, b = params, grid = grid, fullinfo = true)
swf_emp_sub = sum(sum(sub_sim.wdt[t, :, :, 1] .* sub_sim.wd[:wmin]) + sum(sub_sim.wdt[t, :, :, 2] .* sub_sim.wd[:wmin]) for t in (draw+1):(T+draw))/T
swf_unemp_sub = sum(dot(z[sub_sim.statet[t], :], sub_sim.uxt[t, :] .* l) for t in (burn+1):(burn+T))/T

#Sub - MW
submw = parse.(Float64, readlines("output/sub_mw.txt"))

submw_sim = optCrit([submw[1], submw[2]], zeros(N, M), minwage, submw[3]; M = M, N = N, T = T, burn = burn, draw = draw, b = params, grid = grid, fullinfo = true)
swf_emp_submw = sum(sum(submw_sim.wdt[t, :, :, 1] .* submw_sim.wd[:wmin]) + sum(submw_sim.wdt[t, :, :, 2] .* submw_sim.wd[:wmin]) for t in (draw+1):(T+draw))/T
swf_unemp_submw = sum(dot(z[submw_sim.statet[t], :], submw_sim.uxt[t, :] .* l) for t in (burn+1):(burn+T))/T


