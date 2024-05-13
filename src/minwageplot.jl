using Robin2011_RepPackage, Plots, Random, Statistics, StatsPlots, PrettyTables

N = 100
M = 50

T = 6000
burn = 1000

Random.seed!(2077)
draw = rand(T+burn)

b = params_default(; path = "x0.txt")
g = grids(b; N = N, M = M)

sim = optCrit([0,0], zeros(N, M), 0, 0; M = M, N = N, T = T, burn = burn, draw = draw, grid = g, b = b, fullinfo = true)
sim_mw = optCrit([0,0], zeros(N, M), 0.76, 0; M = M, N = N, T = T, burn = burn, draw = draw, grid = g, b = b, fullinfo = true)

wmin = sim.wd[:wmin]
wmax = sim.wd[:wmax]
wdt = sim.wdt[burn+1:T+burn, :, :, :]

wmin_mw = sim_mw.wd[:wmin]
wmax_mw = sim_mw.wd[:wmax]
wdt_mw = sim_mw.wdt[burn+1:T+burn, :, :, :]

wmin_long = []
wmin_dens_long = []
wmax_long = []
wmax_dens_long = []

for i in 1:N
    for j in 1:M
        push!(wmin_long, wmin[i, j])
        push!(wmin_dens_long, sum(wdt[:, i, j, 1]))
        push!(wmax_long, wmax[i, j])
        push!(wmax_dens_long, sum(wdt[:, i, j, 2]))
    end
end

wmin_dens_long = wmin_dens_long./T
wmax_dens_long = wmax_dens_long./T
wage_long = vcat(wmin_long, wmax_long)
wage_dens_long = vcat(wmin_dens_long, wmax_dens_long)


dens = density(wage_long, weights = wage_dens_long_mw, normalize = true, label = "No Minimum Wage", title = "Wage Distribution")
dens_min = density(wmin_long, weights = wmin_dens_long, normalize = true, label = "No Minimum Wage", title = "Wage Distribution (Monopsony Wages)")
dens_max = density(wmax_long, weights = wmax_dens_long, normalize = true, label = "No Minimum Wage", title = "Wage Distribution (Full Surplus Wages)")

wmin_long_mw = []
wmin_dens_long_mw = []
wmax_long_mw = []
wmax_dens_long_mw = []

for i in 1:N
    for j in 1:M
        push!(wmin_long_mw, wmin_mw[i, j])
        push!(wmin_dens_long_mw, sum(wdt_mw[:, i, j, 1]))
        push!(wmax_long_mw, wmax_mw[i, j])
        push!(wmax_dens_long_mw, sum(wdt_mw[:, i, j, 2]))
    end
end

wmin_dens_long_mw = wmin_dens_long_mw./T
wmax_dens_long_mw = wmax_dens_long_mw./T
wage_long_mw = vcat(wmin_long_mw, wmax_long_mw)
wage_dens_long_mw = vcat(wmin_dens_long_mw, wmax_dens_long_mw)

density!(dens, wage_long_mw, weights = wage_dens_long_mw, normalize = true, label = "Minimum Wage")
xlabel!("Wage")
ylabel!("Density")
density!(dens_min, wmin_long_mw, weights = wmin_dens_long_mw, normalize = true, label = "Minimum Wage")
xlabel!("Wage")
ylabel!("Density")
density!(dens_max, wmax_long_mw, weights = wmax_dens_long_mw, normalize = true, label = "Minimum Wage")
xlabel!("Wage")
ylabel!("Density")

ut = sim.ut[burn+1:burn+T]
ut_mw = sim_mw.ut[burn+1:burn+T]

u_mean = round(100* mean(ut); digits = 2)
u_min = round(100* minimum(ut); digits = 2)
u_max = round(100* maximum(ut); digits = 2)

u_mean_mw = round(100* mean(ut_mw); digits = 2)
u_min_mw = round(100* minimum(ut_mw); digits = 2)
u_max_mw = round(100* maximum(ut_mw); digits = 2)

c1 = ["Mean Unemployment Rate", "Minimum Unemployment Rate", "Maximum Unemployment Rate"]
c2 = ["$u_mean%", "$u_min%", "$u_max%"]
c3 = ["$u_mean_mw%", "$u_min_mw%", "$u_max_mw%"]
header = ["", "No Minimum Wage", "Minimum Wage"]

tab = pretty_table([c1 c2 c3]; header = header)
