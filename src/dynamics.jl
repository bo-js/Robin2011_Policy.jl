function unemp_path(S::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, δ::Number = 0.041563759920623)
    M = length(S[1, :])

    uxt = zeros(T, M)

    ## First Period Starts at Steady State
    uxt[1, :] = [(δ/(δ + λ0))*(S[statet[1], m] > 0) + 1 - (S[statet[1], m] > 0) for m in 1:M]

    ## Define Further Periods Recursively
    uxt[2:T, :] = [1 - ((1 - δ)*(1 - uxt[t-1, m]) + λ0 * uxt[t-1, m])*(S[statet[t-1], m] > 0) for t in 2:T, m in 1:M]

    return uxt

end

function unemp_path(S::Matrix, Smin::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, δ::Number = 0.041563759920623)
    M = length(S[1, :])

    uxt = zeros(T, M)

    ## First Period Starts at Steady State
    uxt[1, :] = [(δ/(δ + λ0))*(S[statet[1], m] > max(0, Smin[statet[1], m])) + 1 - (S[statet[1], m] > max(0, Smin[statet[1], m])) for m in 1:M]

    ## Define Further Periods Recursively
    uxt[2:T, :] = [1 - ((1 - δ)*(1 - uxt[t-1, m]) + λ0 * uxt[t-1, m])*(S[statet[t-1], m] > Smin[statet[t-1], m]) for t in 2:T, m in 1:M]

    return uxt

end

function wage_dens_path(S::Matrix, uxt::Matrix, wd::Dict, l::Vector, U::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, λ1::Number = 0.119345383430366, δ::Number = 0.041563759920623)
    M = length(S[1, :])

    N = length(S[:, 1])

    gt = zeros(T, N, M, 2)

    Wmin = wd[:Wmin]
    Wmax = wd[:Wmax]

    # Define First Period from Steady State
    gt[1, statet[1], :, 1] = [(λ0 * (S[statet[1], m] > 0) * uxt[1, m] * l[m])/(1 - (1 - δ) * (S[statet[1], m] > 0) * (1 - λ1)) for m in 1:M]
    
    gt[1, statet[1], :, 2] = [((1 - δ) * λ1 * (1 - uxt[1, m]) * (S[statet[1], m] > 0) * l[m])/(1 - (1 - δ) * (S[statet[1], m] > 0) * (1 - λ1)) for m in 1:M]

    # Define subsequent Periods Recurively
    for t in 2:T
        for i in 1:N
            if statet[t] == i # For wages being assigned in current period's state
                
                gt[t, i, :, 1] = [(S[statet[t], m] > 0) * (λ0 * uxt[t-1, m] * l[m] + (1 - δ) * (1 - λ1) * (
                    gt[t-1, i, m, 1] + sum(
                        
                        (Wmin[statet[t], k, m] - U[statet[t], m] < 0) * gt[t-1, k, m, 1] + 
                        
                        (Wmax[statet[t], k, m] - U[statet[t], m] < 0) * gt[t-1, k, m, 2] for k in 1:N
                    )
                )) for m in 1:M]

                gt[t, i, :, 2] = [
                    (S[statet[t], m] > 0) * (1 - δ) * (
                        λ1 * (1 - uxt[t - 1, m]) * l[m] + (1 - λ1) * (
                            gt[t - 1, i, m, 2] + sum(
                                
                                (Wmin[statet[t], k, m] - U[statet[t], m] > S[statet[t], m]) * gt[t-1, k, m, 1] + 
                        
                                (Wmax[statet[t], k, m] - U[statet[t], m] > S[statet[t], m]) * gt[t-1, k, m, 2] for k in 1:N
                            )
                        )
                    )
                for m in 1:M]

            else # For wages still viable in this state, but not being actively assigned.

                gt[t, i, :, 1] = [
                    (S[statet[t], m] > 0) * (1 - δ) * (1 - λ1) * (0 ≤ Wmin[statet[t], i, m] - U[statet[t], m] ≤ S[statet[t], m])
                    for m in 1:M
                ]

                gt[t, i, :, 2] = [
                    (S[statet[t], m] > 0) * (1 - δ) * (1 - λ1) * (0 ≤ Wmax[statet[t], i, m] - U[statet[t], m] ≤ S[statet[t], m])
                    for m in 1:M
                ]

            end
        end
    end

    return gt

end

function wage_dens_path(S::Matrix,Smin::Matrix, uxt::Matrix, wd::Dict, l::Vector, U::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, λ1::Number = 0.119345383430366, δ::Number = 0.041563759920623)
    M = length(S[1, :])

    N = length(S[:, 1])

    gt = zeros(T, N, M, 2)

    Wmin = wd[:Wmin]
    Wmax = wd[:Wmax]

    # Define First Period from Steady State
    gt[1, statet[1], :, 1] = [(λ0 * (S[statet[1], m] > max(Smin[statet[1], m], 0)) * uxt[1, m] * l[m])/(1 - (1 - δ) * (S[statet[1], m] > max(Smin[statet[1], m], 0)) * (1 - λ1)) for m in 1:M]
    
    gt[1, statet[1], :, 2] = [((1 - δ) * λ1 * (1 - uxt[1, m]) * (S[statet[1], m] > max(Smin[statet[1], m], 0)) * l[m])/(1 - (1 - δ) * (S[statet[1], m] > max(Smin[statet[1], m], 0)) * (1 - λ1)) for m in 1:M]

    # Define subsequent Periods Recurively
    for t in 2:T
        @threads for i in 1:N
            if statet[t] == i # For wages being assigned in current period's state
                
                gt[t, i, :, 1] = [(S[statet[t], m] > max(Smin[statet[t], m], 0)) * (λ0 * uxt[t-1, m] * l[m] + (1 - δ) * (1 - λ1) * (
                    gt[t-1, i, m, 1] + sum(
                        
                        (Wmin[statet[t], k, m] - U[statet[t], m] < max(Smin[statet[t], m], 0)) * gt[t-1, k, m, 1] + 
                        
                        (Wmax[statet[t], k, m] - U[statet[t], m] < max(Smin[statet[t], m], 0)) * gt[t-1, k, m, 2] for k in 1:N
                    )
                )) for m in 1:M]

                gt[t, i, :, 2] = [
                    (S[statet[t], m] > max(Smin[statet[t], m], 0)) * (1 - δ) * (
                        λ1 * (1 - uxt[t - 1, m]) * l[m] + (1 - λ1) * (
                            gt[t - 1, i, m, 2] + sum(
                                
                                (Wmin[statet[t], k, m] - U[statet[t], m] > S[statet[t], m]) * gt[t-1, k, m, 1] + 
                        
                                (Wmax[statet[t], k, m] - U[statet[t], m] > S[statet[t], m]) * gt[t-1, k, m, 2] for k in 1:N
                            )
                        )
                    )
                for m in 1:M]

            else # For wages still viable in this state, but not being actively assigned.

                gt[t, i, :, 1] = [
                    (S[statet[t], m] > 0) * (1 - δ) * (1 - λ1) * (max(Smin[statet[t], m], 0) ≤ Wmin[statet[t], i, m] - U[statet[t], m] ≤ S[statet[t], m])
                    for m in 1:M
                ]

                gt[t, i, :, 2] = [
                    (S[statet[t], m] > 0) * (1 - δ) * (1 - λ1) * (max(Smin[statet[t], m], 0) ≤ Wmax[statet[t], i, m] - U[statet[t], m] ≤ S[statet[t], m])
                    for m in 1:M
                ]

            end
        end
    end

    return gt

end