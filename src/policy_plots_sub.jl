using Robin2011_RepPackage, Random, Plots, Base.Threads, Statistics, LinearAlgebra, Colors, LaTeXStrings

outdir = "output/figures/sub"

subs = parse.(Float64, readlines("output/sub_grid.txt"))
sub_taxes = parse.(Float64, readlines("output/sub_taxes.txt"))
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

low = (statet .≤ statelow)
med = (statelow .< statet .≤ statehigh)
high = (statet .> statehigh)

tlow = [(statet[t] ≤ statelow) for t in 1:T]
tmed = [(statelow < statet[t] ≤ statehigh) for t in 1:T]
thigh = [(statet[t] > statehigh) for t in 1:T]

# Categorise Workers
x_cdf = cumsum(l)
x_low = (x_cdf .≤ (1/20))
x_med = ((1/20) .< x_cdf .≤ (1/10))
x_high = ((1/20) .< x_cdf .≤ (3/20))
# x_high = (x_cdf .> (3/4))

# Production
p = matchprod(x, y; B = b.B, C = b.C)
z = homeprod(x, y; B = b.B, C = b.C, α = b.α, z0 = b.z0)

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

    avmonopwage = zeros(length(sub))
    avmonopwage_low = zeros(length(sub))
    avmonopwage_med = zeros(length(sub))
    avmonopwage_high = zeros(length(sub))

    avhighwage = zeros(length(sub))
    avhighwage_low = zeros(length(sub))
    avhighwage_med = zeros(length(sub))
    avhighwage_high = zeros(length(sub))

    wminsim = [zeros(N, M) for _ in 1:length(sub)]
    wmaxsim = [zeros(N, M) for _ in 1:length(sub)]
    wdtsim = [zeros(T, N, M, 2) for _ in 1:length(sub)]
    uxtsim = [zeros(T, M) for _ in 1:length(sub)]

    q = zeros(length(sub))
    q_low = zeros(length(sub))
    q_med = zeros(length(sub))
    q_high = zeros(length(sub))

    for i in 1:lastindex(sub)
        sim = optCrit([y[statelow+1], sub[i]], zeros(N, M), wmin, subtax[i]; M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g, fullinfo = true)
        sel = burn+1:burn+T
        
        wd = sim.wd
        wminsim[i] = wd[:wmin]
        wmaxsim[i] = wd[:wmax]
        wdt = sim.wdt
        wdt = wdt[sel, :, :, :]
        wdtsim[i] = wdt
        statet = sim.statet
        statet = statet[sel]
        wagext = sim.wagext
        wagext = wagext[sel, :]
        ut = sim.ut
        ut = ut[sel]
        uxt = sim.uxt
        uxt = uxt[sel, :]
        uxtsim[i] = uxt


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

        avwages[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, :]) for t in 1:T)/T
        avwages_low[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, :]) for t in 1:T if statet[t] ≤ statelow)/sum(low)
        avwages_med[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, :]) for t in 1:T if statelow < statet[t] ≤ statehigh)/sum(med)
        avwages_high[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, :]) for t in 1:T if statehigh < statet[t])/sum(high)
        
        avmonopwage[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin])/sum(wdt[t, :, :, 1]) for t in 1:T)/T
        avmonopwage_low[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin])/sum(wdt[t, :, :, 1]) for t in 1:T if statet[t] ≤ statelow)/sum(low)
        avmonopwage_med[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin])/sum(wdt[t, :, :, 1]) for t in 1:T if statelow < statet[t] ≤ statehigh)/sum(med)
        avmonopwage_high[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin])/sum(wdt[t, :, :, 1]) for t in 1:T if statehigh < statet[t])/sum(high)
        
        avhighwage[i] = sum(sum(wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, 2]) for t in 1:T)/T
        avhighwage_low[i] = sum(sum(wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, 2]) for t in 1:T if statet[t] ≤ statelow)/sum(low)
        avhighwage_med[i] = sum(sum(wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, 2]) for t in 1:T if statelow < statet[t] ≤ statehigh)/sum(med)
        avhighwage_high[i] = sum(sum(wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, 2]) for t in 1:T if statehigh < statet[t])/sum(high)

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

        qt = sim.qt[sel]
        q[i] = mean(qt)
        q_low[i] = mean(qt[low])
        q_med[i] = mean(qt[med])
        q_high[i] = mean(qt[high])

    end

    return Dict(
        :wminsim => wminsim,
        :wmaxsim => wmaxsim,
        :wdtsim => wdtsim,
        :uxtsim => uxtsim,

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

        :avmonopwage => avmonopwage,
        :avmonopwage_low => avmonopwage_low,
        :avmonopwage_med => avmonopwage_med,
        :avmonopwage_high => avmonopwage_high,

        :avhighwage => avhighwage,
        :avhighwage_low => avhighwage_low,
        :avhighwage_med => avhighwage_med,
        :avhighwage_high => avhighwage_high,
        
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
        :wageincome_high => wageincome_high,

        :q => q,
        :q_low => q_low,
        :q_med => q_med,
        :q_high => q_high


    )

end

d = policy_sims_sub(subs, sub_taxes)

# Quit Rate
pq = plot(subs, d[:q]; title = "All States", legend = false)
pq_low = plot(subs, d[:q_low], title = "Low State", legend = false)
pq_med = plot(subs, d[:q_med], title = "Medium State", legend = false)
pq_high = plot(subs, d[:q_high], title = "High State", legend = false)

pq_all = plot(pq, pq_low, pq_med, pq_high, plot_title = "Quit Rate", titlefontsize = 12, size = (700, 500))
plot!(xlabel = "Subsidy in Low State", guidefontsize = 8)
savefig("$outdir/quitrate.svg")
savefig("$outdir/quitrate.png")

# Productivity
yt = [y[statet[t]] for t in 1:T]
pt = matchprod(x, yt; B = b.B, C = b.C)

function productivity(uxt; timesel = nothing)
    if isnothing(timesel) == true
        ext = repeat(l', T, 1).*(1 .- uxt)
        et = sum(ext, dims = 2)
        xt = sum(pt.*ext, dims = 2)./et
        return mean(xt)
    else
        ext = repeat(l', sum(timesel), 1).*(1 .- uxt[timesel, :])
        et = sum(ext, dims = 2)
        xt = sum(pt[timesel, :].*ext, dims = 2)./et
        return mean(xt)
    end
end

pprod = plot(subs, productivity.(d[:uxtsim]))
pprod_low = plot(subs, productivity.(d[:uxtsim], timesel = low))
pprod_med = plot(subs, productivity.(d[:uxtsim], timesel = med))
pprod_high = plot(subs, productivity.(d[:uxtsim], timesel = high))
pprod_all = plot(pprod, pprod_low, pprod_med, pprod_high, title = "Labour Productivity", titlefontsize = 12, size = (700, 500))
plot!(xlabel = "Subsidy in Low State", guidefontsize = 8, legend = false)
savefig("$outdir/prod.svg")
savefig("$outdir/prod.png")


# Proportion on Monopsony Wage
pmonop = plot(subs, d[:monop]; title = "All States", xlabel = "Subsidy in Low State", titlefontsize = 12)
pmonop_low = plot(subs, d[:monop_low]; title = "Low State", xlabel = "Subsidy in Low State", titlefontsize = 12)
pmonop_med = plot(subs, d[:monop_med]; title = "Medium State", xlabel = "Subsidy in Low State", titlefontsize = 12)
pmonop_high = plot(subs, d[:monop_high]; title = "High State", xlabel = "Subsidy in Low State", titlefontsize = 12)
pmonop_all = plot(pmonop, pmonop_low, pmonop_med, pmonop_high, layout = 4, legend = false, plot_title = "Proportion of Employees on Monopsony Wage",plot_titlefontsize = 15, plot_titlevspan = 0.06)
plot!(guidefontsize = 8, titlefontsize = 12, plot_titlefontsize = 20, size = (700, 500)) 
savefig("$outdir/monopprop.svg")
savefig("$outdir/monopprop.png")

# No Earning Monopsony Wage and No Earning Full Surplus
nmonop = d[:nmonop]
nmonop = nmonop./nmonop[1]
nfullsurp = d[:nfullsurp]
nfullsurp = nfullsurp/nfullsurp[1]
pwagedens = plot(subs, nmonop; label = "Monopsony Wage", title = "All States", xlabel = "Subsidy in Low State", titlefontsize = 12)
plot!(pwagedens, subs, nfullsurp; label = "Full Surplus Wage")
plot!(legend = false)

nmonop_low = d[:nmonop_low]
nmonop_low = nmonop_low./nmonop_low[1]
nfullsurp_low = d[:nfullsurp_low]
nfullsurp_low = nfullsurp_low/nfullsurp_low[1]
pwagedens_low = plot(subs, nmonop_low; label = "Monopsony Wage", title = "Low State", xlabel = "Subsidy in Low State", titlefontsize = 12)
plot!(pwagedens_low, subs, nfullsurp_low; label = "Full Surplus Wage")
plot!(legend = false)

nmonop_med = d[:nmonop_med]
nmonop_med = nmonop_med./nmonop_med[1]
nfullsurp_med = d[:nfullsurp_med]
nfullsurp_med = nfullsurp_med/nfullsurp_med[1]
pwagedens_med = plot(subs, nmonop_med; label = "Monopsony Wage", title = "Medium State", xlabel = "Subsidy in Low State", titlefontsize = 12)
plot!(pwagedens_med, subs, nfullsurp_med; label = "Full Surplus Wage")
plot!(legend = false)

nmonop_high = d[:nmonop_high]
nmonop_high = nmonop_high./nmonop_high[1]
nfullsurp_high = d[:nfullsurp_high]
nfullsurp_high = nfullsurp_high/nfullsurp_high[1]
pwagedens_high = plot(subs, nmonop_high; label = "Monopsony Wage", title = "High State", xlabel = "Subsidy in Low State", titlefontsize = 12)
plot!(pwagedens_high, subs, nfullsurp_high; label = "Full Surplus Wage")
plot!(legend = :bottomleft)

pwagedens_all = plot(pwagedens, pwagedens_low, pwagedens_med, pwagedens_high, layout = 4,size = (700, 600), plot_title = "Wage Densities", titlefontsize = 12, plot_titlevspan = 0.06, guidefontsize = 8)
savefig("$outdir/wagedens.svg")
savefig("$outdir/wagedens.png")

# Average Wages Across States
pavwages = plot(subs, d[:avwages]; title = "All States", xlabel = "Subsidy in Low State")
pavwages_low = plot(subs, d[:avwages_low]; title = "Low State", xlabel = "Subsidy in Low State")
pavwages_med = plot(subs, d[:avwages_med]; title = "Medium State", xlabel = "Subsidy in Low State")
pavwages_high = plot(subs, d[:avwages_high]; title = "High State", xlabel = "Subsidy in Low State")
pavwages_all = plot(pavwages, pavwages_low, pavwages_med, pavwages_high, layout = 4)
plot!(size = (700, 500), guidefontsize = 8, titlefontsize = 12, plot_title = "Mean Wages", legend = false)
savefig("$outdir/avwage.svg")
savefig("$outdir/avwage.png")

pavlowwage = plot(subs, d[:avmonopwage];label = "Average Monopsony/Minimum Wage", title = "Average Wage by Wage Type", xlabel = "Subsidy in Low State" )
pavhighwage = plot(subs, d[:avhighwage], color = "orange", label = "Average Full Surplus Wage")
xlabel!("Subsidy in Low State")
pavwagebytype = plot(pavlowwage, pavhighwage, layout = [1;1], size = (500, 600))
savefig("$outdir/avwagetype.svg")
savefig("$outdir/avwagetype.png")

# Looks pretty much the same across States - so just do the above
# pavlowwagebytype_low = plot(subs, d[:avmonopwage_low];label = "Average Low Wage", title = "Low State", xlabel = "Subsidy in Low State")
# plot!(xlabel = "Subsidy in Low State",legend = false)
# pavhighwagebytype = plot(subs, d[:avhighwage_low];color = "orange", label = "Average High Wage", title  = "Low State")
# plot!(xlabel = "Subsidy in Low State")

# # plot!(legend = :outertop)

# pavlowwagebytype_med = plot(subs, d[:avmonopwage_med];label = "Average Low Wage", title = "Medium State", xlabel = "Subsidy in Low State")
# plot!(legend = false)
# pavhighwagebytype_ned = plot(subs, d[:avhighwage_med],color = "orange", label = "Average High Wage")
# # plot!(legend = false)

# pavlowwagebytype_high = plot(subs, d[:avmonopwage_high];label = "Average Low Wage", title = "High State", xlabel = "Subsidy in Low State")
# plot!(legend = false)
# pavlowwagebytype_high = plot( subs, d[:avhighwage_high],color = "orange" ,label = "Average High Wage")
# # plot!(legend = false)

# pavlowwagesbytype_all = plot(pavlowwagebytype, pavlowwagebytype_low, pavlowwagebytype_med, pavlowwagebytype_high, layout = 4, plot_titlefontsize = 15, plot_titlevspan = 0.06)
# plot!(guidefontsize = 8, plot_title = "Mean Low Wage")

# Unemployment Rate
pu = plot(subs, d[:mean_u]; title = "All States", xlabel = "Subsidy in Low State")
pu_low = plot(subs, d[:mean_u_low]; title = "Low State", xlabel = "Subsidy in Low State")
pu_med = plot(subs, d[:mean_u_med]; title = "Medium State", xlabel = "Subsidy in Low State")
pu_high = plot(subs, d[:mean_u_high]; title = "High State", xlabel = "Subsidy in Low State")
pu_all = plot(pu, pu_low, pu_med, pu_high, layout = 4)
plot!(legend = false, plot_title = "Mean Unemployment Rate", titlefontsize = 12, size = (700, 500))
savefig("$outdir/unemp.svg")
savefig("$outdir/unemp.png")

# Welfare and Income
wageinc = d[:wageincome]
wageinc = wageinc./wageinc[1]
valadd = d[:vanet]
valadd = valadd/valadd[1]
home = d[:homeprod]
home = home./home[1]
pinc = plot(subs, wageinc;label = "Wage Income", title = "All States", xlabel = "Subsidy in Low State")
plot!(pinc, subs, valadd; label = "Net Market Production")
plot!(legend = false)
# plot!(pinc, subs, home, label = "Home Production")

wageinc_low = d[:wageincome_low]
wageinc_low = wageinc_low./wageinc_low[1]
valadd_low = d[:vanet_low]
valadd_low = valadd_low/valadd_low[1]
home_low = d[:homeprod_low]
home_low = home_low./home_low[1]
pinc_low = plot(subs, wageinc_low;label = "Wage Income", title = "Low State", xlabel = "Subsidy in Low State")
plot!(pinc_low, subs, valadd_low; label = "Net Market Production")
plot!(legend = :bottomright)
# plot!(pinc_low, subs, home_low; label = "Home Production")

wageinc_med = d[:wageincome_med]
wageinc_med = wageinc_med./wageinc_med[1]
valadd_med = d[:vanet_med]
valadd_med = valadd_med/valadd_med[1]
home_med = d[:homeprod_med]
home_med = home_med./home_med[1]
pinc_med = plot(subs, wageinc_med;label = "Wage Income", title = "Medium State", xlabel = "Subsidy in Low State")
plot!(pinc_med, subs, valadd_med; label = "Net Market Production")
plot!(legend = false)
# plot!(pinc_med, subs, home_med; label = "Home Production")

wageinc_high = d[:wageincome_high]
wageinc_high = wageinc_high./wageinc_high[1]
valadd_high = d[:vanet_high]
valadd_high = valadd_high/valadd_high[1]
home_high = d[:homeprod_high]
home_high = home_high./home_high[1]
pinc_high = plot(subs, wageinc_high;label = "Wage Income", title = "High State", xlabel = "Subsidy in Low State")
plot!(pinc_high, subs, valadd_high; label = "Net Market Production")
plot!(legend = false)
# plot!(pinc_high, subs, home_high; label = "Home Production")

pinc_all = plot(pinc, pinc_low, pinc_med, pinc_high; plot_title = "Production and Income", size = (700, 500))
plot!(titlefontsize = 13, plot_titlevspan = 0.06, guidefontsize = 8)
savefig("$outdir/inc.svg")
savefig("$outdir/inc.png")

# Change in Unemployment Rate and Wages by Worker Type and State
uxtsim = d[:uxtsim]
wminsim = d[:wminsim]
wmaxsim = d[:wmaxsim]
wdtsim = d[:wdtsim]

u_l = [mean(uxtsim[s][:, x_low]*(l[x_low]./sum(l[x_low]))) for s in 1:length(subs)]
u_m = [mean(uxtsim[s][:, x_med]*(l[x_med]./sum(l[x_med]))) for s in 1:length(subs)]
u_h = [mean(uxtsim[s][:, x_high]*(l[x_high]./sum(l[x_high]))) for s in 1:length(subs)]

pu_disag = plot(subs, u_l./u_l[1], label = "First Ventile")
plot!(subs, u_m./u_m[1], label = "Second Ventile")
plot!(subs, u_h./u_h[1], label = "Third Ventile")
plot!(title = "All States", legend = false)

u_ll = [mean(uxtsim[s][low, x_low]*(l[x_low]./sum(l[x_low]))) for s in 1:length(subs)]
u_lm = [mean(uxtsim[s][low, x_med]*(l[x_med]./sum(l[x_med]))) for s in 1:length(subs)]
u_lh = [mean(uxtsim[s][low, x_high]*(l[x_high]./sum(l[x_high]))) for s in 1:length(subs)]

pu_disag_l = plot(subs, u_ll./u_ll[1], label = "First Ventile")
plot!(subs, u_lm./u_lm[1], label = "Second Ventile")
plot!(subs, u_lh./u_lh[1], label = "Third Ventile")
plot!(title = "Low State", legend = false)

u_ml = [mean(uxtsim[s][med, x_low]*(l[x_low]./sum(l[x_low]))) for s in 1:length(subs)]
u_mm = [mean(uxtsim[s][med, x_med]*(l[x_med]./sum(l[x_med]))) for s in 1:length(subs)]
u_mh = [mean(uxtsim[s][med, x_high]*(l[x_high]./sum(l[x_high]))) for s in 1:length(subs)]

pu_disag_m = plot(subs, u_ml./u_ml[1], label = "First Ventile")
plot!(subs, u_mm./u_mm[1], label = "Second Ventile")
plot!(subs, u_mh./u_mh[1], label = "Third Ventile")
plot!(title = "Medium State", legend = false)

u_hl = [mean(uxtsim[s][high, x_low]*(l[x_low]./sum(l[x_low]))) for s in 1:length(subs)]
u_hm = [mean(uxtsim[s][high, x_med]*(l[x_med]./sum(l[x_med]))) for s in 1:length(subs)]
u_hh = [mean(uxtsim[s][high, x_high]*(l[x_high]./sum(l[x_high]))) for s in 1:length(subs)]

pu_disag_h = plot(subs, u_hl./u_hl[1], label = "First Ventile")
plot!(subs, u_hm./u_hm[1], label = "Second Ventile")
plot!(subs, u_hh./u_hh[1], label = "Third Ventile")
plot!(title = "High State", legend = :bottomright)

pu_disag_all = plot(pu_disag, pu_disag_l, pu_disag_m, pu_disag_h, layout = 4)
plot!(plot_title = "Change in Unemployment Rate by Worker Skill", size = (700, 500), xlabel = "Subsidy in Low State")
savefig("$outdir/unemp_skill.svg")
savefig("$outdir/unemp_skill.png")

# Avg Wages by State x Worker Type
# Define Worker Quintiles
x_q1 = (x_cdf .≤ 1/5)
x_q2 = (1/5 .< x_cdf .≤ 2/5)
x_q3 = (2/5 .< x_cdf .≤ 3/5)
x_q4 = (3/5 .< x_cdf .≤ 4/5)
x_q5 = (4/5 .< x_cdf)

wagebyskill_q1 = [
    mean(
        (sum(wdtsim[s][t, :, x_q1, 1] .* wminsim[s][:, x_q1]) + sum(wdtsim[s][t, :, x_q1, 2] .* wmaxsim[s][:, x_q1]))/sum(wdtsim[s][t, :, x_q1, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q1, :]) != 0) for s in 1:length(subs)
]

wagebyskill_q2 = [
    mean(
        (sum(wdtsim[s][t, :, x_q2, 1] .* wminsim[s][:, x_q2]) + sum(wdtsim[s][t, :, x_q2, 2] .* wmaxsim[s][:, x_q2]))/sum(wdtsim[s][t, :, x_q2, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q2, :]) != 0) for s in 1:length(subs)
]

wagebyskill_q3 = [
    mean(
        (sum(wdtsim[s][t, :, x_q3, 1] .* wminsim[s][:, x_q3]) + sum(wdtsim[s][t, :, x_q3, 2] .* wmaxsim[s][:, x_q3]))/sum(wdtsim[s][t, :, x_q3, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q3, :]) != 0) for s in 1:length(subs)
]

wagebyskill_q4 = [
    mean(
        (sum(wdtsim[s][t, :, x_q4, 1] .* wminsim[s][:, x_q4]) + sum(wdtsim[s][t, :, x_q4, 2] .* wmaxsim[s][:, x_q4]))/sum(wdtsim[s][t, :, x_q4, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q4, :]) != 0) for s in 1:length(subs)
]

wagebyskill_q5 = [
    mean(
        (sum(wdtsim[s][t, :, x_q5, 1] .* wminsim[s][:, x_q5]) + sum(wdtsim[s][t, :, x_q5, 2] .* wmaxsim[s][:, x_q5]))/sum(wdtsim[s][t, :, x_q5, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q5, :]) != 0) for s in 1:length(subs)
]

pwagebyskill = plot(subs, wagebyskill_q1./wagebyskill_q1[1], label = "First Quintile")
plot!(subs, wagebyskill_q2./wagebyskill_q2[1], label = "Second Quintile")
plot!(subs, wagebyskill_q3./wagebyskill_q3[1], label = "Third Quintile")
plot!(subs, wagebyskill_q4./wagebyskill_q4[1], label = "Fourth Quintile")
plot!(subs, wagebyskill_q5./wagebyskill_q5[1], label = "Fifth Quintile")
plot!(title = "All States")

wagebyskill_lq1 = [
    mean(
        (sum(wdtsim[s][t, :, x_q1, 1] .* wminsim[s][:, x_q1]) + sum(wdtsim[s][t, :, x_q1, 2] .* wmaxsim[s][:, x_q1]))/sum(wdtsim[s][t, :, x_q1, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q1, :]) != 0 && tlow[t] == 1) for s in 1:length(subs)
]

wagebyskill_lq2 = [
    mean(
        (sum(wdtsim[s][t, :, x_q2, 1] .* wminsim[s][:, x_q2]) + sum(wdtsim[s][t, :, x_q2, 2] .* wmaxsim[s][:, x_q2]))/sum(wdtsim[s][t, :, x_q2, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q2, :]) != 0 && tlow[t] == 1) for s in 1:length(subs)
]

wagebyskill_lq3 = [
    mean(
        (sum(wdtsim[s][t, :, x_q3, 1] .* wminsim[s][:, x_q3]) + sum(wdtsim[s][t, :, x_q3, 2] .* wmaxsim[s][:, x_q3]))/sum(wdtsim[s][t, :, x_q3, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q3, :]) != 0 && tlow[t] == 1) for s in 1:length(subs)
]

wagebyskill_lq4 = [
    mean(
        (sum(wdtsim[s][t, :, x_q4, 1] .* wminsim[s][:, x_q4]) + sum(wdtsim[s][t, :, x_q4, 2] .* wmaxsim[s][:, x_q4]))/sum(wdtsim[s][t, :, x_q4, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q4, :]) != 0 && tlow[t] == 1) for s in 1:length(subs)
]

wagebyskill_lq5 = [
    mean(
        (sum(wdtsim[s][t, :, x_q5, 1] .* wminsim[s][:, x_q5]) + sum(wdtsim[s][t, :, x_q5, 2] .* wmaxsim[s][:, x_q5]))/sum(wdtsim[s][t, :, x_q5, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q5, :]) != 0 && tlow[t] == 1) for s in 1:length(subs)
]

pwagebyskill_l = plot(subs, wagebyskill_lq1./wagebyskill_lq1[1], label = "First Quintile")
plot!(subs, wagebyskill_lq2./wagebyskill_lq2[1], label = "Second Quintile")
plot!(subs, wagebyskill_lq3./wagebyskill_lq3[1], label = "Third Quintile")
plot!(subs, wagebyskill_lq4./wagebyskill_lq4[1], label = "Fourth Quintile")
plot!(subs, wagebyskill_lq5./wagebyskill_lq5[1], label = "Fifth Quintile")
plot!(title = "Low State", legend = false)


wagebyskill_mq1 = [
    mean(
        (sum(wdtsim[s][t, :, x_q1, 1] .* wminsim[s][:, x_q1]) + sum(wdtsim[s][t, :, x_q1, 2] .* wmaxsim[s][:, x_q1]))/sum(wdtsim[s][t, :, x_q1, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q1, :]) != 0 && tmed[t] == 1) for s in 1:length(subs)
]

wagebyskill_mq2 = [
    mean(
        (sum(wdtsim[s][t, :, x_q2, 1] .* wminsim[s][:, x_q2]) + sum(wdtsim[s][t, :, x_q2, 2] .* wmaxsim[s][:, x_q2]))/sum(wdtsim[s][t, :, x_q2, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q2, :]) != 0 && tmed[t] == 1) for s in 1:length(subs)
]

wagebyskill_mq3 = [
    mean(
        (sum(wdtsim[s][t, :, x_q3, 1] .* wminsim[s][:, x_q3]) + sum(wdtsim[s][t, :, x_q3, 2] .* wmaxsim[s][:, x_q3]))/sum(wdtsim[s][t, :, x_q3, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q3, :]) != 0 && tmed[t] == 1) for s in 1:length(subs)
]

wagebyskill_mq4 = [
    mean(
        (sum(wdtsim[s][t, :, x_q4, 1] .* wminsim[s][:, x_q4]) + sum(wdtsim[s][t, :, x_q4, 2] .* wmaxsim[s][:, x_q4]))/sum(wdtsim[s][t, :, x_q4, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q4, :]) != 0 && tmed[t] == 1) for s in 1:length(subs)
]

wagebyskill_mq5 = [
    mean(
        (sum(wdtsim[s][t, :, x_q5, 1] .* wminsim[s][:, x_q5]) + sum(wdtsim[s][t, :, x_q5, 2] .* wmaxsim[s][:, x_q5]))/sum(wdtsim[s][t, :, x_q5, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q5, :]) != 0 && tmed[t] == 1) for s in 1:length(subs)
]

pwagebyskill_m = plot(subs, wagebyskill_mq1./wagebyskill_mq1[1], label = "First Quintile")
plot!(subs, wagebyskill_mq2./wagebyskill_mq2[1], label = "Second Quintile")
plot!(subs, wagebyskill_mq3./wagebyskill_mq3[1], label = "Third Quintile")
plot!(subs, wagebyskill_mq4./wagebyskill_mq4[1], label = "Fourth Quintile")
plot!(subs, wagebyskill_mq5./wagebyskill_mq5[1], label = "Fifth Quintile")
plot!(title = "Medium State", legend = false)


wagebyskill_hq1 = [
    mean(
        (sum(wdtsim[s][t, :, x_q1, 1] .* wminsim[s][:, x_q1]) + sum(wdtsim[s][t, :, x_q1, 2] .* wmaxsim[s][:, x_q1]))/sum(wdtsim[s][t, :, x_q1, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q1, :]) != 0 && thigh[t] == 1) for s in 1:length(subs)
]

wagebyskill_hq2 = [
    mean(
        (sum(wdtsim[s][t, :, x_q2, 1] .* wminsim[s][:, x_q2]) + sum(wdtsim[s][t, :, x_q2, 2] .* wmaxsim[s][:, x_q2]))/sum(wdtsim[s][t, :, x_q2, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q2, :]) != 0 && thigh[t] == 1) for s in 1:length(subs)
]

wagebyskill_hq3 = [
    mean(
        (sum(wdtsim[s][t, :, x_q3, 1] .* wminsim[s][:, x_q3]) + sum(wdtsim[s][t, :, x_q3, 2] .* wmaxsim[s][:, x_q3]))/sum(wdtsim[s][t, :, x_q3, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q3, :]) != 0 && thigh[t] == 1) for s in 1:length(subs)
]

wagebyskill_hq4 = [
    mean(
        (sum(wdtsim[s][t, :, x_q4, 1] .* wminsim[s][:, x_q4]) + sum(wdtsim[s][t, :, x_q4, 2] .* wmaxsim[s][:, x_q4]))/sum(wdtsim[s][t, :, x_q4, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q4, :]) != 0 && thigh[t] == 1) for s in 1:length(subs)
]

wagebyskill_hq5 = [
    mean(
        (sum(wdtsim[s][t, :, x_q5, 1] .* wminsim[s][:, x_q5]) + sum(wdtsim[s][t, :, x_q5, 2] .* wmaxsim[s][:, x_q5]))/sum(wdtsim[s][t, :, x_q5, :]
    ) 
    for t in 1:T if sum(wdtsim[s][t, :, x_q5, :]) != 0 && thigh[t] == 1) for s in 1:length(subs)
]

pwagebyskill_h = plot(subs, wagebyskill_hq1./wagebyskill_hq1[1], label = "First Quintile")
plot!(subs, wagebyskill_hq2./wagebyskill_hq2[1], label = "Second Quintile")
plot!(subs, wagebyskill_hq3./wagebyskill_hq3[1], label = "Third Quintile")
plot!(subs, wagebyskill_hq4./wagebyskill_hq4[1], label = "Fourth Quintile")
plot!(subs, wagebyskill_hq5./wagebyskill_hq5[1], label = "Fifth Quintile")
plot!(title = "High State", legend = false)


pwagebyskill_all = plot(pwagebyskill, pwagebyskill_l, pwagebyskill_m, pwagebyskill_h, layout = 4)
plot!(xlabel = "Subsidy in Low State", plot_title = "Effect on Average Wages by Worker Skill", size = (700, 500), guidefontsize = 8)

savefig("$outdir/avwage_skill.svg")
savefig("$outdir/avwage_skill.png")

## Income by Skill
incbyskill_q1 = [
    mean(
        sum(wdtsim[s][t, :, x_q1, 1] .* wminsim[s][:, x_q1]) + sum(wdtsim[s][t, :, x_q1, 2] .* wmaxsim[s][:, x_q1])
    
    for t in 1:T) for s in 1:length(subs)
]

incbyskill_q2 = [
    mean(
        sum(wdtsim[s][t, :, x_q2, 1] .* wminsim[s][:, x_q2]) + sum(wdtsim[s][t, :, x_q2, 2] .* wmaxsim[s][:, x_q2])
    
    for t in 1:T) for s in 1:length(subs)
]

incbyskill_q3 = [
    mean(
        sum(wdtsim[s][t, :, x_q3, 1] .* wminsim[s][:, x_q3]) + sum(wdtsim[s][t, :, x_q3, 2] .* wmaxsim[s][:, x_q3])
   
    for t in 1:T) for s in 1:length(subs)
]

incbyskill_q4 = [
    mean(
        sum(wdtsim[s][t, :, x_q4, 1] .* wminsim[s][:, x_q4]) + sum(wdtsim[s][t, :, x_q4, 2] .* wmaxsim[s][:, x_q4])
    
    for t in 1:T) for s in 1:length(subs)
]

incbyskill_q5 = [
    mean(
        sum(wdtsim[s][t, :, x_q5, 1] .* wminsim[s][:, x_q5]) + sum(wdtsim[s][t, :, x_q5, 2] .* wmaxsim[s][:, x_q5])
    
    for t in 1:T) for s in 1:length(subs)
]

pincbyskill = plot(subs, incbyskill_q1./incbyskill_q1[1], label = "First Quintile")
plot!(subs, incbyskill_q2./incbyskill_q2[1], label = "Second Quintile")
plot!(subs, incbyskill_q3./incbyskill_q3[1], label = "Third Quintile")
plot!(subs, incbyskill_q4./incbyskill_q4[1], label = "Fourth Quintile")
plot!(subs, incbyskill_q5./incbyskill_q5[1], label = "Fifth Quintile")
plot!(title = "All States", legend = false)

incbyskill_lq1 = [
    mean(
        sum(wdtsim[s][t, :, x_q1, 1] .* wminsim[s][:, x_q1]) + sum(wdtsim[s][t, :, x_q1, 2] .* wmaxsim[s][:, x_q1])
    for t in 1:T if tlow[t] == 1) for s in 1:length(subs)
]

incbyskill_lq2 = [
    mean(
        sum(wdtsim[s][t, :, x_q2, 1] .* wminsim[s][:, x_q2]) + sum(wdtsim[s][t, :, x_q2, 2] .* wmaxsim[s][:, x_q2])
    
    for t in 1:T if tlow[t] == 1) for s in 1:length(subs)
]

incbyskill_lq3 = [
    mean(
        sum(wdtsim[s][t, :, x_q3, 1] .* wminsim[s][:, x_q3]) + sum(wdtsim[s][t, :, x_q3, 2] .* wmaxsim[s][:, x_q3])
    
    for t in 1:T if tlow[t] == 1) for s in 1:length(subs)
]

incbyskill_lq4 = [
    mean(
        sum(wdtsim[s][t, :, x_q4, 1] .* wminsim[s][:, x_q4]) + sum(wdtsim[s][t, :, x_q4, 2] .* wmaxsim[s][:, x_q4])
    
    for t in 1:T if tlow[t] == 1) for s in 1:length(subs)
]

incbyskill_lq5 = [
    mean(
        sum(wdtsim[s][t, :, x_q5, 1] .* wminsim[s][:, x_q5]) + sum(wdtsim[s][t, :, x_q5, 2] .* wmaxsim[s][:, x_q5])
     
    for t in 1:T if tlow[t] == 1) for s in 1:length(subs)
]

pincbyskill_l = plot(subs, incbyskill_lq1./incbyskill_lq1[1], label = "First Quintile")
plot!(subs, incbyskill_lq2./incbyskill_lq2[1], label = "Second Quintile")
plot!(subs, incbyskill_lq3./incbyskill_lq3[1], label = "Third Quintile")
plot!(subs, incbyskill_lq4./incbyskill_lq4[1], label = "Fourth Quintile")
plot!(subs, incbyskill_lq5./incbyskill_lq5[1], label = "Fifth Quintile")
plot!(title = "Low State", legend = :topleft)


incbyskill_mq1 = [
    mean(
        sum(wdtsim[s][t, :, x_q1, 1] .* wminsim[s][:, x_q1]) + sum(wdtsim[s][t, :, x_q1, 2] .* wmaxsim[s][:, x_q1])
     
    for t in 1:T if tmed[t] == 1) for s in 1:length(subs)
]

incbyskill_mq2 = [
    mean(
        sum(wdtsim[s][t, :, x_q2, 1] .* wminsim[s][:, x_q2]) + sum(wdtsim[s][t, :, x_q2, 2] .* wmaxsim[s][:, x_q2])
     
    for t in 1:T if tmed[t] == 1) for s in 1:length(subs)
]

incbyskill_mq3 = [
    mean(
        sum(wdtsim[s][t, :, x_q3, 1] .* wminsim[s][:, x_q3]) + sum(wdtsim[s][t, :, x_q3, 2] .* wmaxsim[s][:, x_q3])
    
    for t in 1:T if tmed[t] == 1) for s in 1:length(subs)
]

incbyskill_mq4 = [
    mean(
        sum(wdtsim[s][t, :, x_q4, 1] .* wminsim[s][:, x_q4]) + sum(wdtsim[s][t, :, x_q4, 2] .* wmaxsim[s][:, x_q4])
    
    for t in 1:T if tmed[t] == 1) for s in 1:length(subs)
]

incbyskill_mq5 = [
    mean(
        sum(wdtsim[s][t, :, x_q5, 1] .* wminsim[s][:, x_q5]) + sum(wdtsim[s][t, :, x_q5, 2] .* wmaxsim[s][:, x_q5])
     
    for t in 1:T if tmed[t] == 1) for s in 1:length(subs)
]

pincbyskill_m = plot(subs, incbyskill_mq1./incbyskill_mq1[1], label = "First Quintile")
plot!(subs, incbyskill_mq2./incbyskill_mq2[1], label = "Second Quintile")
plot!(subs, incbyskill_mq3./incbyskill_mq3[1], label = "Third Quintile")
plot!(subs, incbyskill_mq4./incbyskill_mq4[1], label = "Fourth Quintile")
plot!(subs, incbyskill_mq5./incbyskill_mq5[1], label = "Fifth Quintile")
plot!(title = "Medium State", legend = false)


incbyskill_hq1 = [
    mean(
        sum(wdtsim[s][t, :, x_q1, 1] .* wminsim[s][:, x_q1]) + sum(wdtsim[s][t, :, x_q1, 2] .* wmaxsim[s][:, x_q1])
     
    for t in 1:T if thigh[t] == 1) for s in 1:length(subs)
]

incbyskill_hq2 = [
    mean(
        sum(wdtsim[s][t, :, x_q2, 1] .* wminsim[s][:, x_q2]) + sum(wdtsim[s][t, :, x_q2, 2] .* wmaxsim[s][:, x_q2])
     
    for t in 1:T if thigh[t] == 1) for s in 1:length(subs)
]

incbyskill_hq3 = [
    mean(
        sum(wdtsim[s][t, :, x_q3, 1] .* wminsim[s][:, x_q3]) + sum(wdtsim[s][t, :, x_q3, 2] .* wmaxsim[s][:, x_q3])
     
    for t in 1:T if thigh[t] == 1) for s in 1:length(subs)
]

incbyskill_hq4 = [
    mean(
        sum(wdtsim[s][t, :, x_q4, 1] .* wminsim[s][:, x_q4]) + sum(wdtsim[s][t, :, x_q4, 2] .* wmaxsim[s][:, x_q4])
     
    for t in 1:T if thigh[t] == 1) for s in 1:length(subs)
]

incbyskill_hq5 = [
    mean(
        sum(wdtsim[s][t, :, x_q5, 1] .* wminsim[s][:, x_q5]) + sum(wdtsim[s][t, :, x_q5, 2] .* wmaxsim[s][:, x_q5])
     
    for t in 1:T if thigh[t] == 1) for s in 1:length(subs)
]

pincbyskill_h = plot(subs, incbyskill_hq1./incbyskill_hq1[1], label = "First Quintile")
plot!(subs, incbyskill_hq2./incbyskill_hq2[1], label = "Second Quintile")
plot!(subs, incbyskill_hq3./incbyskill_hq3[1], label = "Third Quintile")
plot!(subs, incbyskill_hq4./incbyskill_hq4[1], label = "Fourth Quintile")
plot!(subs, incbyskill_hq5./incbyskill_hq5[1], label = "Fifth Quintile")
plot!(title = "High State", legend = false)


pincbyskill_all = plot(pincbyskill, pincbyskill_l, pincbyskill_m, pincbyskill_h, layout = 4)
plot!(xlabel = "Subsidy in Low State", plot_title = "Effect on Wage Income by Worker Skill", size = (700, 500), guidefontsize = 8)

savefig("$outdir/wageinc_skill.svg")
savefig("$outdir/wageinc_skill.png")

## Zoom In

incbyskillzoom_low = [
    mean(
        sum(wdtsim[s][t, :, x_low, 1] .* wminsim[s][:, x_low]) + sum(wdtsim[s][t, :, x_low, 2] .* wmaxsim[s][:, x_low])
    
    for t in 1:T) for s in 1:length(subs)
]

incbyskillzoom_med = [
    mean(
        sum(wdtsim[s][t, :, x_med, 1] .* wminsim[s][:, x_med]) + sum(wdtsim[s][t, :, x_med, 2] .* wmaxsim[s][:, x_med])
    
    for t in 1:T) for s in 1:length(subs)
]

incbyskillzoom_high = [
    mean(
        sum(wdtsim[s][t, :, x_high, 1] .* wminsim[s][:, x_high]) + sum(wdtsim[s][t, :, x_high, 2] .* wmaxsim[s][:, x_high])
   
    for t in 1:T) for s in 1:length(subs)
]

pincbyskillzoom = plot(subs, incbyskillzoom_low./incbyskillzoom_low[1], label = "First Ventile")
plot!(subs, incbyskillzoom_med./incbyskillzoom_med[1], label = "Second Ventile")
plot!(subs, incbyskillzoom_high./incbyskillzoom_high[1], label = "Third Ventile")

plot!(title = "All States")

incbyskillzoom_llow = [
    mean(
        sum(wdtsim[s][t, :, x_low, 1] .* wminsim[s][:, x_low]) + sum(wdtsim[s][t, :, x_low, 2] .* wmaxsim[s][:, x_low])
    for t in 1:T if tlow[t] == 1) for s in 1:length(subs)
]

incbyskillzoom_lmed = [
    mean(
        sum(wdtsim[s][t, :, x_med, 1] .* wminsim[s][:, x_med]) + sum(wdtsim[s][t, :, x_med, 2] .* wmaxsim[s][:, x_med])
    
    for t in 1:T if tlow[t] == 1) for s in 1:length(subs)
]

incbyskillzoom_lhigh = [
    mean(
        sum(wdtsim[s][t, :, x_high, 1] .* wminsim[s][:, x_high]) + sum(wdtsim[s][t, :, x_high, 2] .* wmaxsim[s][:, x_high])
    
    for t in 1:T if tlow[t] == 1) for s in 1:length(subs)
]


pincbyskillzoom_l = plot(subs, incbyskillzoom_llow./incbyskillzoom_llow[1], label = "First Ventile")
plot!(subs, incbyskillzoom_lmed./incbyskillzoom_lmed[1], label = "Second Ventile")
plot!(subs, incbyskillzoom_lhigh./incbyskillzoom_lhigh[1], label = "Third Ventile")

plot!(title = "Low State", legend = false)


incbyskillzoom_mlow = [
    mean(
        sum(wdtsim[s][t, :, x_low, 1] .* wminsim[s][:, x_low]) + sum(wdtsim[s][t, :, x_low, 2] .* wmaxsim[s][:, x_low])
     
    for t in 1:T if tmed[t] == 1) for s in 1:length(subs)
]

incbyskillzoom_mmed = [
    mean(
        sum(wdtsim[s][t, :, x_med, 1] .* wminsim[s][:, x_med]) + sum(wdtsim[s][t, :, x_med, 2] .* wmaxsim[s][:, x_med])
     
    for t in 1:T if tmed[t] == 1) for s in 1:length(subs)
]

incbyskillzoom_mhigh = [
    mean(
        sum(wdtsim[s][t, :, x_high, 1] .* wminsim[s][:, x_high]) + sum(wdtsim[s][t, :, x_high, 2] .* wmaxsim[s][:, x_high])
    
    for t in 1:T if tmed[t] == 1) for s in 1:length(subs)
]


pincbyskillzoom_m = plot(subs, incbyskillzoom_mlow./incbyskillzoom_mlow[1], label = "First Ventile")
plot!(subs, incbyskillzoom_mmed./incbyskillzoom_mmed[1], label = "Second Ventile")
plot!(subs, incbyskillzoom_mhigh./incbyskillzoom_mhigh[1], label = "Third Ventile")
plot!(title = "Medium State", legend = false)


incbyskillzoom_hlow = [
    mean(
        sum(wdtsim[s][t, :, x_low, 1] .* wminsim[s][:, x_low]) + sum(wdtsim[s][t, :, x_low, 2] .* wmaxsim[s][:, x_low])
     
    for t in 1:T if thigh[t] == 1) for s in 1:length(subs)
]

incbyskillzoom_hmed = [
    mean(
        sum(wdtsim[s][t, :, x_med, 1] .* wminsim[s][:, x_med]) + sum(wdtsim[s][t, :, x_med, 2] .* wmaxsim[s][:, x_med])
     
    for t in 1:T if thigh[t] == 1) for s in 1:length(subs)
]

incbyskillzoom_hhigh = [
    mean(
        sum(wdtsim[s][t, :, x_high, 1] .* wminsim[s][:, x_high]) + sum(wdtsim[s][t, :, x_high, 2] .* wmaxsim[s][:, x_high])
     
    for t in 1:T if thigh[t] == 1) for s in 1:length(subs)
]


pincbyskillzoom_h = plot(subs, incbyskillzoom_hlow./incbyskillzoom_hlow[1], label = "First Ventile")
plot!(subs, incbyskillzoom_hmed./incbyskillzoom_hmed[1], label = "Second Ventile")
plot!(subs, incbyskillzoom_hhigh./incbyskillzoom_hhigh[1], label = "Third Ventile")

plot!(title = "High State", legend = false)


pincbyskillzoom_all = plot(pincbyskillzoom, pincbyskillzoom_l, pincbyskillzoom_m, pincbyskillzoom_h, layout = 4)
plot!(xlabel = "Subsidy in Low State", plot_title = "Effect on Wage Income by Worker Skill", size = (700, 500), guidefontsize = 8)

savefig("$outdir/wageinc_skillzoom.svg")
savefig("$outdir/wageinc_skillzoom.png")

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
    plot!(subs, y, linecolor = c[i], label = "χ = $w", legend = :bottomleft)
end

xlabel!("Subsidy in Low State")
ylabel!("Weighted Social Welfare (Normalised)")
title!("Social Welfare")

savefig("$outdir/swf.svg")
savefig("$outdir/swf.png")