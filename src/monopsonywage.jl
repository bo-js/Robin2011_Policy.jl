using Robin2011_RepPackage, Random, Plots, Base.Threads

subs = parse.(Float64, readlines("output/sub_grid.txt"))
sub_taxes = parse.(Float64, readlines("output/sub_taxes.txt"))

N = 50
M = 50

b = params_default()
g = grids(p; N = N, M = M)

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

function monopsony_sub(sub, subtax; wmin = 0)
    
    monop = zeros(length(sub))
    monop_low = zeros(length(sub))
    monop_med = zeros(length(sub))
    monop_high = zeros(length(sub))

    @threads for i in 1:lastindex(sub)
        sim = optCrit([y[statelow], sub[i]], zeros(N, M), wmin, subtax[i]; M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g, fullinfo = true)
        sel = burn+1:burn+T
        
        wdt = sim.wdt
        wdt = wdt[sel, :, :, :]
        statet = sim.statet
        statet = statet[sel]

        low = (statet .≤ statelow)
        med = (statelow .< statet .≤ statehigh)
        high = (statet .> statehigh)

        monop[i] = sum(wdt[:, :, :, 1])/sum(wdt[:, :, :, :])

        monop_low[i] = sum(wdt[low, :, :, 1])/sum(wdt[low, :, :, :])
        monop_med[i] = sum(wdt[med, :, :, 1])/sum(wdt[med, :, :, :])
        monop_high[i] = sum(wdt[high, :, :, 1])/sum(wdt[high, :, :, :])

    end

    return (monop = monop, low = monop_low, med = monop_med, high = monop_high)

end

monpsony_props = monopsony_sub(subs, sub_taxes)

plot(subs, monpsony_props.monop)
plot(subs, monpsony_props.low)
plot(subs, monpsony_props.med)
plot(subs, monpsony_props.high)

