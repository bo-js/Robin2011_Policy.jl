function WageVFI(S::Matrix, Π::Matrix, z::Matrix; λ1 = 0.119345383430366, β::Number = 0.946603693905558)
    
    N = length(Π[1, :])
    M = length(S[1, :])

    Wmin = zeros(N, N, M)
    Wmax = zeros(N, N, M)
    wmin = zeros(N, M)
    wmax = zeros(N, M)
    kf = zeros(N)
    kf .= λ1
    for m in 1:M
        for i in 1:N

            # Minimum Wages - W is a vector giving the worker surplus of being in each state with the minimum wage set for state i
            A = (Π - repeat(Π[i, :]', N, 1)) * repeat((S[:, m] .> 0)', N, 1)
            W = (I(length(Π[1, :])) - β * A * repeat((1 .- kf)', N, 1)) \ (z[i, m] .- z[:, m] + β * A *( kf .* S[:, m]))
            W0 = min.(max.(W, 0), S[:, m])
            e = norm(W - W0, 2)

            while e > 0.001
                W1 = W
                W = z[i, m] .- z[:, m] + β * A * (kf .* S[:, m] + (1 .- kf) .* W0)
                W0 = min.(max.(W, 0), S[:, m])
                e = norm(W - W1, 2)
            end

            Wmin[:, i, m] = W
            wmin[i, m] = z[i, m] - first(β * Π[i, :]' .* (S[:, m] .> 0)' * (kf .* S[:, m] + (1 .- kf) .* W0))


            ## Maximum Wages - Likewise as above

            W = (I(length(Π[1, :])) - β * A * repeat((1 .- kf)', N, 1)) \ (S[i, m] + z[i, m] .- z[:, m] + β * A * (kf .* S[:, m]))
            W0 = min.(max.(W, 0), S[:, m])
            e = norm(W - W0, 2)

            while e > 0.001
                W1 = W
                W = S[i, m] + z[i, m] .- z[:, m] + β * A * (kf .* S[:, m] + (1 .- kf) .* W0)
                W0 = min.(max.(W, 0), S[:, m])
                e = norm(W - W1, 2)
            end

            Wmax[:, i, m] = W
            wmax[i, m] = S[i, m] + z[i, m] - first(β * Π[i, :]' .* (S[:, m] .> 0)' * (kf .* S[:, m] + (1 .- kf) .* W0))


        end
    end
 
    return Dict(:Wmin => Wmin, :wmin => wmin, :Wmax => Wmax, :wmax => wmax)
end

function WageVFI(S::Matrix, Smin::Matrix, Π::Matrix, z::Matrix; λ1 = 0.119345383430366, β::Number = 0.946603693905558)
    
    N = length(Π[1, :])
    M = length(S[1, :])

    Wmin = zeros(N, N, M)
    Wmax = zeros(N, N, M)
    wmin = zeros(N, M)
    wmax = zeros(N, M)
    kf = zeros(N)
    kf .= λ1
    for m in 1:M
        for i in 1:N

            # Minimum Wages - W is a vector giving the worker surplus of being in each state with the minimum wage set for state i
            A = (Π - repeat(Π[i, :]', N, 1)) * repeat((S[:, m] .> max.(Smin[:, m], 0))', N, 1)
            W = (I(length(Π[1, :])) - β * A * repeat((1 .- kf)', N, 1)) \ (max(Smin[i, m], 0) + z[i, m] .- z[:, m] + β * A *( kf .* max.(Smin[:, m], S[:, m])))
            W0 = min.(max.(W, 0, Smin[:, m]), max.(Smin[:, m], S[:, m]))
            e = norm(W - W0, 2)

            while e > 0.001
                W1 = W
                W = max(Smin[i,m], 0) + z[i, m] .- z[:, m] + β * A * (kf .* max.(Smin[:, m], S[:, m]) + (1 .- kf) .* W0)
                W0 = min.(max.(W, 0, Smin[:, m]), max.(Smin[:, m], S[:, m]))
                e = norm(W - W1, 2)
            end

            Wmin[:, i, m] = W
            wmin[i, m] = max(0, Smin[i, m]) + z[i, m] - first(β * Π[i, :]' .* (S[:, m] .> max.(Smin[:, m], 0))' * (kf .* max.(Smin[:, m], S[:, m]) + (1 .- kf) .* W0))


            ## Maximum Wages - Likewise as above

            W = (I(length(Π[1, :])) - β * A * repeat((1 .- kf)', N, 1)) \ (max.(Smin[i, m], S[i, m]) + z[i, m] .- z[:, m] + β * A * (kf .* S[:, m]))
            W0 = min.(max.(W, Smin[:,m], 0), S[:, m])
            e = norm(W - W0, 2)

            while e > 0.001
                W1 = W
                W = max.(Smin[i, m], S[i, m]) + z[i, m] .- z[:, m] + β * A * (kf .* S[:, m] + (1 .- kf) .* W0)
                W0 = min.(max.(W, Smin[:,m], 0), max.(Smin[:, m], S[:, m]))
                e = norm(W - W1, 2)
            end

            Wmax[:, i, m] = W
            wmax[i, m] = max.(Smin[i, m], S[i, m])+ z[i, m] - first(β * Π[i, :]' .* (S[:, m] .> max.(Smin[:, m], 0))' * (kf .* max.(Smin[:, m], S[:, m]) + (1 .- kf) .* W0))


        end
    end
 
    return Dict(:Wmin => Wmin, :wmin => wmin, :Wmax => Wmax, :wmax => wmax)
end