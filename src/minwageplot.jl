using Robin2011_RepPackage, Plots, Random, Statistics, StatsPlots, PrettyTables, LaTeXStrings, LinearAlgebra

outdir = "output/figures"

N = 100
M = 50

T = 6000
burn = 1000

Random.seed!(2077)
draw = rand(T+burn)

statet = states(draw, Π)
statet = statet[burn+1:burn+T]
statesort = sort(statet)
statelow = maximum(statesort[1:2000])
statehigh = maximum(statesort[2000:4000])

low = (statet .≤ statelow)
med = (statelow .< statet .≤ statehigh)
high = (statet .> statehigh)

tlow = [(statet[t] ≤ statelow) for t in 1:T]
tmed = [(statelow < statet[t] ≤ statehigh) for t in 1:T]
thigh = [(statet[t] > statehigh) for t in 1:T]

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

k = 0.001
dens = density(wage_long, weights = wage_dens_long, normalize = true, label = "No Min Wage", title = "Wage Distribution", kernel = k)
dens_min = density(wmin_long, weights = wmin_dens_long, normalize = true, label = "No Min Wage", title = "Monopsony Wages", legend = false, kernel = k)
dens_max = density(wmax_long, weights = wmax_dens_long, normalize = true, label = "No Min Wage", title = "Full Surplus Wages", legend = false, kernel = k)

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

density!(dens, wage_long_mw, weights = wage_dens_long_mw, normalize = true, label = "Min Wage = 0.76", kernel = k)
xlabel!("Wage")
ylabel!("Density")
density!(dens_min, wmin_long_mw, weights = wmin_dens_long_mw, normalize = true, label = "Min Wage = 0.76", kernel = k)
xlabel!("Wage")
ylabel!("Density")
density!(dens_max, wmax_long_mw, weights = wmax_dens_long_mw, normalize = true, label = "Min Wage = 0.76", kernel = k)
xlabel!("Wage")
ylabel!("Density")



plot(dens, dens_min, dens_max; layout = @layout[a;b c])

savefig("$outdir/mwplot.svg")
savefig("$outdir/mwplot.png")


ut = sim.ut[burn+1:burn+T]
ut_mw = sim_mw.ut[burn+1:burn+T]

u_mean = round(100* mean(ut); digits = 2)
u_min = round(100* minimum(ut); digits = 2)
u_max = round(100* maximum(ut); digits = 2)
u_low = round(100* mean(ut[tlow]); digits = 2)
u_high = round(100* mean(ut[thigh]); digits = 2)
u_med = round(100* mean(ut[tmed]); digits = 2)


u_mean_mw = round(100* mean(ut_mw); digits = 2)
u_min_mw = round(100* minimum(ut_mw); digits = 2)
u_max_mw = round(100* maximum(ut_mw); digits = 2)
u_low_mw = round(100* mean(ut_mw[tlow]); digits = 2)
u_high_mw = round(100* mean(ut_mw[thigh]); digits = 2)
u_med_mw = round(100* mean(ut_mw[tmed]); digits = 2)


c1 = ["Mean Unemployment Rate", "Mean Unemp. Rate in Low State", "Mean Unemp. Rate in Medium State", "Mean Unemp. Rate in High State", "Minimum Unemployment Rate", "Maximum Unemployment Rate"]
c2 = ["$u_mean%", "$u_low%", "$u_med%", "$u_high%", "$u_min%", "$u_max%"]
c3 = ["$u_mean_mw%","$u_low_mw%", "$u_med_mw%", "$u_high_mw%", "$u_min_mw%", "$u_max_mw%"]
header = ["", "No Minimum Wage", "Minimum Wage = 0.76"]

tab = pretty_table(String, [c1 c2 c3]; header = header, backend = Val(:latex))
file = open("$outdir/minwagetab.tex", "w")
write(file, tab)
close(file)
