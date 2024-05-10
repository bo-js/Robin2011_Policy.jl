using Robin2011_RepPackage, Random, Plots, Base.Threads, Statistics

subs = parse.(Float64, readlines("output/sub_grid.txt"))
sub_taxes = parse.(Float64, readlines("output/sub_taxes.txt"))

N = 50
M = 50

b = params_default()
g = grids(b; N = N, M = M)

x = g[:x]
y = g[:y]
Π = g[:Π]
l = g[:l]

T = 6000
burn = 1000

Random.seed!(2077)
draw = rand(T+burn)

# Categorise States
statet = states(draw, Π)
statet = statet[burn+1:burn+T]
statesort = sort(statet)
statelow = maximum(statesort[1:2000])
statehigh = maximum(statesort[2000:4000])

# Categorise Workers
x_cdf = cumsum(l)
x_low = maximum(x[(x_cdf .≤ (1/3))])
x_high = maximum(x[(x_cdf .≤ (2/3))])

function policy_sims_sub(sub, subtax; wmin = 0)
    
    monop = zeros(length(sub))
    monop_low = zeros(length(sub))
    monop_med = zeros(length(sub))
    monop_high = zeros(length(sub))

    nmonop = zeros(length(sub))
    nmonop_low = zeros(length(sub))
    nmonop_med = zeros(length(sub))
    nmonop_high = zeros(length(sub))

    nfullsurp = zeros(length(sub))
    nfullsurp_low = zeros(length(sub))
    nfullsurp_med = zeros(length(sub))
    nfullsurp_high = zeros(length(sub))

    avwages = zeros((length(sub)))
    avwages_low = zeros((length(sub)))
    avwages_med = zeros((length(sub)))
    avwages_high = zeros((length(sub)))

    mean_u = zeros((length(sub)))
    mean_u_low = zeros((length(sub)))
    mean_u_med = zeros((length(sub)))
    mean_u_high = zeros((length(sub)))

    @threads for i in 1:lastindex(sub)
        sim = optCrit([y[statelow], sub[i]], zeros(N, M), wmin, subtax[i]; M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g, fullinfo = true)
        sel = burn+1:burn+T
        
        wdt = sim.wdt
        wdt = wdt[sel, :, :, :]
        statet = sim.statet
        statet = statet[sel]
        wagext = sim.wagext
        wagext = wagext[sel, :]
        ut = sim.ut
        ut = ut[sel]

        low = (statet .≤ statelow)
        med = (statelow .< statet .≤ statehigh)
        high = (statet .> statehigh)

        monop[i] = sum(wdt[:, :, :, 1])/sum(wdt[:, :, :, :])
        monop_low[i] = sum(wdt[low, :, :, 1])/sum(wdt[low, :, :, :])
        monop_med[i] = sum(wdt[med, :, :, 1])/sum(wdt[med, :, :, :])
        monop_high[i] = sum(wdt[high, :, :, 1])/sum(wdt[high, :, :, :])

        nmonop[i] = sum(wdt[:, :, :, 1])
        nmonop_low[i] = sum(wdt[low, :, :, 1])
        nmonop_med[i] = sum(wdt[med, :, :, 1])
        nmonop_high[i] = sum(wdt[high, :, :, 1])
        
        nfullsurp[i] = sum(wdt[:, :, :, 2])
        nfullsurp_low[i] = sum(wdt[low, :, :, 2])
        nfullsurp_med[i] = sum(wdt[med, :, :, 2])
        nfullsurp_high[i] = sum(wdt[high, :, :, 2])

        avwages[i] = mean(wagext * l)
        avwages_low[i] = mean(wagext[low, :] * l)
        avwages_med[i] = mean(wagext[med, :] * l)
        avwages_high[i] = mean(wagext[high, :] * l)
        
        mean_u[i] = mean(ut)
        mean_u_low[i] = mean(ut[low])
        mean_u_med[i] = mean(ut[med])
        mean_u_high[i] = mean(ut[high])


    end

    return Dict(
        :monop => monop, 
        :monop_low => monop_low, 
        :monop_med => monop_med, 
        :monop_high => monop_high,
        :nmonop => nmonop,
        :nmonop_low => nmonop_low, 
        :nmonop_med => nmonop_med, 
        :nmonop_high => nmonop_high,
        :nfullsurp => nfullsurp,
        :nfullsurp_low => nfullsurp_low, 
        :nfullsurp_med => nfullsurp_med, 
        :nfullsurp_high => nfullsurp_high,
        :avwages => avwages,
        :avwages_low => avwages_low,
        :avwages_med => avwages_med,
        :avwages_high => avwages_high,
        :mean_u => mean_u,
        :mean_u_low => mean_u_low,
        :mean_u_med => mean_u_med,
        :mean_u_high => mean_u_high
    )

end

d = policy_sims_sub(subs, sub_taxes)

# Proportion on Monopsony Wage
pmonop = plot(subs, d[:monop]; title = "Proportion on Monopsony Wage", xlabel = "Subsidy in Low State")
pmonop_low = plot(subs, d[:monop_low]; title = "Proportion on Monopsony Wage (Low State)", xlabel = "Subsidy in Low State")
pmonop_med = plot(subs, d[:monop_med]; title = "Proportion on Monopsony Wage (Medium State)", xlabel = "Subsidy in Low State")
pmonop_high = plot(subs, d[:monop_high]; title = "Proportion on Monopsony Wage (High State)", xlabel = "Subsidy in Low State")

# No Earning Monopsony Wage and No Earning Full Surplus
nmonop = d[:nmonop]
nmonop = nmonop./nmonop[1]
nfullsurp = d[:nfullsurp]
nfullsurp = nfullsurp/nfullsurp[1]
pwagedens = plot(subs, nmonop; label = "Monopsony Wage", title = "Wage Densities", xlabel = "Subsidy in Low State")
plot!(pwagedens, subs, nfullsurp; label = "Full Surplus Wage")

nmonop_low = d[:nmonop_low]
nmonop_low = nmonop_low./nmonop_low[1]
nfullsurp_low = d[:nfullsurp_low]
nfullsurp_low = nfullsurp_low/nfullsurp_low[1]
pwagedens_low = plot(subs, nmonop_low; label = "Monopsony Wage", title = "Wage Densities (Low State)", xlabel = "Subsidy in Low State")
plot!(pwagedens_low, subs, nfullsurp_low; label = "Full Surplus Wage")

nmonop_med = d[:nmonop_med]
nmonop_med = nmonop_med./nmonop_med[1]
nfullsurp_med = d[:nfullsurp_med]
nfullsurp_med = nfullsurp_med/nfullsurp_med[1]
pwagedens_med = plot(subs, nmonop_med; label = "Monopsony Wage", title = "Wage Densities (Medium State)", xlabel = "Subsidy in Low State")
plot!(pwagedens_med, subs, nfullsurp_med; label = "Full Surplus Wage")

nmonop_high = d[:nmonop_high]
nmonop_high = nmonop_high./nmonop_high[1]
nfullsurp_high = d[:nfullsurp_high]
nfullsurp_high = nfullsurp_high/nfullsurp_high[1]
pwagedens_high = plot(subs, nmonop_high; label = "Monopsony Wage", title = "Wage Densities (High State)", xlabel = "Subsidy in Low State")
plot!(pwagedens_high, subs, nfullsurp_high; label = "Full Surplus Wage")


# Average Wages Across States
pavwages = plot(subs, d[:avwages]; title = "Mean Wage", xlabel = "Subsidy in Low State")
pavwages_low = plot(subs, d[:avwages_low]; title = "Mean Wage (Low State)", xlabel = "Subsidy in Low State")
pavwages_med = plot(subs, d[:avwages_med]; title = "Mean Wage (Medium State)", xlabel = "Subsidy in Low State")
pavwages_high = plot(subs, d[:avwages_high]; title = "Mean Wage (High State)", xlabel = "Subsidy in Low State")

# Unemployment Rate
pu = plot(subs, d[:mean_u]; title = "Mean Unemployment Rate", xlabel = "Subsidy in Low State")
pu_low = plot(subs, d[:mean_u_low]; title = "Mean Unemployment Rate (Low State)", xlabel = "Subsidy in Low State")
pu_med = plot(subs, d[:mean_u_med]; title = "Mean Unemployment Rate (Medium State)", xlabel = "Subsidy in Low State")
pu_high = plot(subs, d[:mean_u_high]; title = "Mean Unemployment Rate (High State)", xlabel = "Subsidy in Low State")


