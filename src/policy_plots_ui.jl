using Robin2011_RepPackage, Random, Plots, Base.Threads, Statistics

subs = parse.(Float64, readlines("output/sub_grid.txt"))
sub_taxes = parse.(Float64, readlines("output/sub_taxes.txt"))

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

# UI

ui = parse.(Float64, readlines("output/ui_grid.txt"))
ui_taxes = parse.(Float64, readlines("output/ui_taxes.txt"))


p = matchprod(x, y; B = b.B, C = b.C)

function policy_sims_ui(ui, uitax; wmin = 0)
    
    monop = zeros(length(ui))
    monop_low = zeros(length(ui))
    monop_med = zeros(length(ui))
    monop_high = zeros(length(ui))

    avwages = zeros((length(ui)))
    avwages_low = zeros((length(ui)))
    avwages_med = zeros((length(ui)))
    avwages_high = zeros((length(ui)))

    mean_u = zeros((length(ui)))
    mean_u_low = zeros((length(ui)))
    mean_u_med = zeros((length(ui)))
    mean_u_high = zeros((length(ui)))

    @threads for i in 1:lastindex(ui)
        sim = optCrit([0, 0], uiyx(ui[i], p), wmin, uitax[i]; M = M, N = N, T = T, burn = burn, draw = draw, b = b, grid = g, fullinfo = true)
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

        avwages[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, :]) for t in 1:T)/T
        avwages_low[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, :]) for t in 1:T if statet[t] ≤ statelow)/sum(low)
        avwages_med[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, :]) for t in 1:T if statelow < statet[t] ≤ statehigh)/sum(med)
        avwages_high[i] = sum(sum(wdt[t, :, :, 1] .* wd[:wmin] + wdt[t, :, :, 2] .* wd[:wmax])/sum(wdt[t, :, :, :]) for t in 1:T if statehigh < statet[t])/sum(high)
        
        
        mean_u[i] = mean(ut)
        mean_u_low[i] = mean(ut[low])
        mean_u_med[i] = mean(ut[med])
        mean_u_high[i] = mean(ut[high])


    end

    return (
        monop = monop, 
        monop_low = monop_low, 
        monop_med = monop_med, 
        monop_high = monop_high, 
        avwages = avwages,
        avwages_low = avwages_low,
        avwages_med = avwages_med,
        avwages_high = avwages_high,
        mean_u = mean_u,
        mean_u_low = mean_u_low,
        mean_u_med = mean_u_med,
        mean_u_high = mean_u_high
    )

end

policy_effects_ui = policy_sims_ui(ui, ui_taxes)

pmonopui = plot(ui, policy_effects_ui.monop; title = "Proportion on Monopsony Wage", xlabel ="UI Replacement Rate" )
pmonopui_low = plot(ui, policy_effects_ui.monop_low; title = "Proportion on Monopsony Wage (Low State)", xlabel ="UI Replacement Rate")
pmonopui_med = plot(ui, policy_effects_ui.monop_med; title = "Proportion on Monopsony Wage (Medium State)", xlabel ="UI Replacement Rate")
pmonopui_high = plot(ui, policy_effects_ui.monop_high; title = "Proportion on Monopsony Wage (High State)", xlabel ="UI Replacement Rate")

pavwageui = plot(ui, policy_effects_ui.avwages; title = "Mean Wages", xlabel ="UI Replacement Rate" )
pavwageui_low = plot(ui, policy_effects_ui.avwages_low; title = "Mean Wage (Low State)", xlabel ="UI Replacement Rate")
pavwageui_med = plot(ui, policy_effects_ui.avwages_med; title = "Mean Wage (Medium State)", xlabel ="UI Replacement Rate")
pavwageui_high = plot(ui, policy_effects_ui.avwages_high; title = "Mean (High State)", xlabel ="UI Replacement Rate")

