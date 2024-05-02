using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Robin2011_RepPackage

M = 500
N = 100

r::Number = 0.05/4

## Parameter Values - Change for when we estimate
δ::Number = 0.041563759920623
λ0::Number = 0.994544861919718
λ1::Number = 0.119345383430366
ρ = 0.913702234286476
σ = 0.0257
ν::Number = 2.019365636076711
μ::Number = 5.786082109731152
τ::Number = 0.5


dta = CSV.read("/Users/bojs/Desktop/Robin 2011 Rep Files/matlab/USquarterly.csv", DataFrame, header = false)

filter!(row -> !isnan(row.Column1) && !isnan(row.Column4), dta);

prod = collect(dta.Column1);

prod = prod/mean(prod);

T = length(prod);

wage = collect(dta.Column2);
unr = collect(dta.Column3);
vac = collect(dta.Column4);
tightness = collect(dta.Column5);
f_q = collect(dta.Column6);
s_q = collect(dta.Column7);
year = collect(dta.Column8);
qtr = collect(dta.Column9);

time = year + qtr./4;

### Define Grids
grid = grids(; M = M, N = N, ν = ν, μ = μ, ρ = ρ, σ = σ)
x = grid[:x]
y = grid[:y]
Π = grid[:Π]
F = grid[:F]
l = grid[:l]

##Production and Surplus
p = matchprod(x, y)
z = homeprod(x, y)
Sx = SurplusVFI(p, z, Π)
Ux = (I(N) - Π./(1 + r))\z

## Wages
wd = WageVFI(Sx, Π, z)

##Steady State
ux = (δ/(δ + λ0)) .* (Sx .> 0) + (Sx .<= 0)

## Initial Conditions
burn = 0;
T1 = burn+T;
Ft = ones(T1)
uxt = ones(T1+1, M)
yt = ones(T1)
statet = zeros(Int, T1)
Wqua=zeros(T1,9);
wqua=zeros(T1,9);
# wminqua=zeros(T1,9);
# wmaxqua=zeros(T1,9);
wagemean=zeros(T1);
wagevar=zeros(T1);
meanwvar=zeros(T1);
varwmean=zeros(T1);

i = min(sum(y.-prod[1] .<= 0), N);
Ft[1] = F[i];
yt[1] = y[i];
statet[1] = i;
uxt[1, :] = ux[i,:]

## Productivity and Unemployment Dynamics
for t in 1:(burn+T)

    uxt_1  = 1 .- (Sx .> 0) .* repeat((1 - δ) * (1 .- uxt[t, :]) + λ0 * uxt[t, :], 1, N)'
    ext = repeat(l', N, 1) .* (1 .- uxt_1)
    et = sum(ext, dims = 2)
    xt = vec(sum(p .* ext, dims = 2)./et)
        
    global i = argmin(abs.(xt .- prod[t]))
    statet[t] = i
    yt[t] = y[i]

    uxt[t+1, :] = [1 - (Sx[i, m] > 0) * ((1 - δ) * (1 - uxt[t, m]) + λ0 * uxt[t, m]) for m in 1:M]
end

ut = uxt * l
u = mean(ut)
## Wage Dynamics
wdt = wage_dens_path(Sx, uxt, wd, l, Ux, statet, T1)

paidwagest = [
    sum(wdt[t, :, :, 1] .* wd[:wmin]) + sum(wdt[t, :, :, 2] .* wd[:wmin]) for t in 1:length(prod)
]

avwagest = paidwagest./(1 .- ut[1:length(prod)])

avwages = mean(avwagest)

propmonopsony = sum(wdt[:, :, :, 1])/sum(wdt)

## Turnover Dynamics
ft = [λ0 * sum((Sx[statet[t], m] > 0) * uxt[t, m] * l[m] for m in 1:M)/ut[t] for t in 1:T1]
qt = [τ * λ1 * (1 - δ) * sum((Sx[statet[t], m] > 0) * (1 - uxt[t, m]) * l[m] for m in 1:M)/(1 - ut[t]) for t in 1:T1]
st = [δ + (1 - δ) * sum((Sx[statet[t], m] ≤ 0) * (1 - uxt[t, m]) * l[m] for m in 1:M)/(1 - ut[t]) for t in 1:T1]
