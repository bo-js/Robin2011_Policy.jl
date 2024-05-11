using Robin2011_RepPackage, Random, Plots, Base.Threads, Statistics, LinearAlgebra, Colors

subs = parse.(Float64, readlines("output/sub_grid_mw.txt"))
sub_taxes = parse.(Float64, readlines("output/sub_taxes_mw.txt"))
pushfirst!(subs, 0)
pushfirst!(sub_taxes, 0)

N = 100
M = 50

b = params_default(; path = "x0.txt")
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

# Production
p = matchprod(x, y; B = b.B, C = b.C)
z = homeprod(x, y; B = b.B, C = b.C, α = b.α, z0 = b.z0)

function policy_sims_sub(sub, subtax; wmin = 0, tol = 0.0001)
    
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

    vanet = zeros(length(sub))
    vanet_low = zeros(length(sub))
    vanet_med = zeros(length(sub))
    vanet_high = zeros(length(sub))

    homeprod = zeros(length(sub))
    homeprod_low = zeros(length(sub))
    homeprod_med = zeros(length(sub))
    homeprod_high = zeros(length(sub))

    wageincome = zeros(length(sub))
    wageincome_low = zeros(length(sub))
    wageincome_med = zeros(length(sub))
    wageincome_high = zeros(length(sub))


    @threads for i in 1:lastindex(sub)
        sim = optCrit([y[statelow+1], sub[i]], zeros(N, M), wmin, subtax[i];tol = tol, M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g, fullinfo = true)
        sel = burn+1:burn+T
        
        wd = sim.wd
        wdt = sim.wdt
        wdt = wdt[sel, :, :, :]
        statet = sim.statet
        statet = statet[sel]
        wagext = sim.wagext
        wagext = wagext[sel, :]
        ut = sim.ut
        ut = ut[sel]
        uxt = sim.uxt
        uxt = uxt[sel, :]

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

        avwages[i] = sum((wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/wdt[t, :, :, :] for t in 1:T)/T
        avwages_low[i] = sum((wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/wdt[t, :, :, :] for t in 1:T if statet[t] ≤ statelow)/sum(low)
        avwages_med[i] = sum((wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/wdt[t, :, :, :] for t in 1:T if statelow < statet[t] ≤ statehigh)/sum(med)
        avwages_high[i] = sum((wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/wdt[t, :, :, :] for t in 1:T if statehigh < statet[t])/sum(high)
        
        mean_u[i] = mean(ut)
        mean_u_low[i] = mean(ut[low])
        mean_u_med[i] = mean(ut[med])
        mean_u_high[i] = mean(ut[high])

        vanet[i] = sum((((1 - subtax[i]).*p[statet[t], :] .+ (statet[t] ≤ statelow)*sub[i]) .* (1 .- uxt[t, :])) ⋅ l for t in 1:T)/T
        vanet_low[i] = sum((((1 - subtax[i]).*p[statet[t], :] .+ sub[i]) .* (1 .- uxt[t, :])) ⋅ l for t in 1:T if statet[t] ≤ statelow)/sum(low)
        vanet_med[i] = sum((((1 - subtax[i]).*p[statet[t], :]) .* (1 .- uxt[t, :])) ⋅ l for t in 1:T if statelow < statet[t] ≤ statehigh)/sum(med)
        vanet_high[i] = sum((((1 - subtax[i]).*p[statet[t], :]) .* (1 .- uxt[t, :])) ⋅ l for t in 1:T if statehigh < statet[t])/sum(high)

        homeprod[i] = sum(z[statet[t], :] ⋅ (uxt[t, :] .* l) for t in 1:T)/T
        homeprod_low[i] = sum(z[statet[t], :] ⋅ (uxt[t, :] .* l) for t in 1:T if statet[t] ≤ statelow)/sum(low)
        homeprod_med[i] = sum(z[statet[t], :] ⋅ (uxt[t, :] .* l) for t in 1:T if statelow < statet[t] ≤ statehigh)/sum(med)
        homeprod_high[i] = sum(z[statet[t], :] ⋅ (uxt[t, :] .* l) for t in 1:T if statehigh < statet[t])/sum(high)

        wageincome[i] = sum(sum((wdt[t, :, :, 1] .* wd[:wmin]) ) + sum(wdt[t, :, :, 2] .* wd[:wmax]) for t in 1:T)/T
        wageincome_low[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin]) + sum(wdt[t, :, :, 2] .* wd[:wmax]) for t in 1:T if statet[t] ≤ statelow)/sum(low)
        wageincome_med[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin]) + sum(wdt[t, :, :, 2] .* wd[:wmax]) for t in 1:T if statelow < statet[t] ≤ statehigh)/sum(med)
        wageincome_high[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin]) + sum(wdt[t, :, :, 2] .* wd[:wmax]) for t in 1:T if statehigh < statet[t])/sum(high)



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
        :mean_u_high => mean_u_high,
        
        :vanet => vanet,
        :vanet_low => vanet_low,
        :vanet_med => vanet_med,
        :vanet_high => vanet_high,

        :homeprod => homeprod,
        :homeprod_low => homeprod_low,
        :homeprod_med => homeprod_med,
        :homeprod_high => homeprod_high,

        :wageincome => wageincome,
        :wageincome_low => wageincome_low,
        :wageincome_med => wageincome_med,
        :wageincome_high => wageincome_high

    )

end

d = policy_sims_sub(subs, sub_taxes; wmin = 0.76, tol = 0.008)

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

# Welfare and Income
wageinc = d[:wageincome]
wageinc = wageinc./wageinc[1]
valadd = d[:vanet]
valadd = valadd/valadd[1]
home = d[:homeprod]
home = home./home[1]
pinc = plot(subs, wageinc;label = "Wage Income", title = "Production and Income", xlabel = "Subsidy in Low State")
plot!(pinc, subs, valadd; label = "Net Market Production")
plot!(pinc, subs, home, label = "Home Production")

wageinc_low = d[:wageincome_low]
wageinc_low = wageinc_low./wageinc_low[1]
valadd_low = d[:vanet_low]
valadd_low = valadd_low/valadd_low[1]
home_low = d[:homeprod_low]
home_low = home_low./home_low[1]
pinc_low = plot(subs, wageinc_low;label = "Wage Income", title = "Production and Income (Low State)", xlabel = "Subsidy in Low State")
plot!(pinc_low, subs, valadd_low; label = "Net Market Production")
plot!(pinc_low, subs, home_low; label = "Home Production")

wageinc_med = d[:wageincome_med]
wageinc_med = wageinc_med./wageinc_med[1]
valadd_med = d[:vanet_med]
valadd_med = valadd_med/valadd_med[1]
home_med = d[:homeprod_med]
home_med = home_med./home_med[1]
pinc_med = plot(subs, wageinc_med;label = "Wage Income", title = "Production and Income (Medium State)", xlabel = "Subsidy in med State")
plot!(pinc_med, subs, valadd_med; label = "Net Market Production")
plot!(pinc_med, subs, home_med; label = "Home Production")

wageinc_high = d[:wageincome_high]
wageinc_high = wageinc_high./wageinc_high[1]
valadd_high = d[:vanet_high]
valadd_high = valadd_high/valadd_high[1]
home_high = d[:homeprod_high]
home_high = home_high./home_high[1]
pinc_high = plot(subs, wageinc_high;label = "Wage Income", title = "Production and Income (High State)", xlabel = "Subsidy in high State")
plot!(pinc_high, subs, valadd_high; label = "Net Market Production")
plot!(pinc_high, subs, home_high; label = "Home Production")

#SWF
swf = d[:vanet] + d[:homeprod]
swf = swf/swf[1]

swf_workers = d[:homeprod] + d[:wageincome]
swf_workers = swf_workers/swf_workers[1]

plot(subs, swf)
plot(subs, swf_workers)
plot!(twinx(), subs, d[:mean_u])

pswf_weighted = plot()
weights = LinRange(0, 1, 11)
c = colormap("Blues", length(weights))

for i in 1:lastindex(weights)
    w = weights[i]
    y = d[:homeprod] + d[:wageincome] + w*(d[:vanet] - d[:wageincome])
    y = y/y[1]
    plot!(subs, y, linecolor = c[i], label = "w = $w")
end

xlabel!("Subsidy in Low State")
ylabel!("Weighted Social Welfare (Normalised)")