function SminVFI(wmin, p::Matrix, z::Matrix, Π::Matrix;tol = 0.001, MaxIter = 10000, β::Number = 0.946603693905558, λ1::Number = 0.119345383430366, λ0::Number = 0.994544861919718, r = 0.05/4)
    S = (I(length(Π[1,:])) - β * Π )\(p - z)

    Smin = wmin .-z + β * λ1 .* Π * ((S .> 0) .* S)

    e1 = norm(S - max.(S, Smin, 0), 2)

    e2 = norm(Smin, 2)

    for iter in 1:MaxIter
        S1 = S
        Smin1 = Smin
        S = p - z + β * Π * ((S1 .> max.(Smin1, 0)) .* S1) - (1/(1+r)) .* Π * (λ0 .* (S1 .> max.(Smin1, 0)) .* max.(Smin1, 0))
        e1 = norm(S - S1, 2)
        Smin = wmin .- z + β * Π * ((S .> max.(Smin, 0)) .* (λ1 .* max.(S, Smin) + (1 - λ1) .* max.(Smin, 0))) - (1/(1+r)) .* Π * (λ1 .* (S .> max.(Smin, 0)) .* max.(Smin, 0))
        e2 = norm(Smin - Smin1, 2)
        if e1 ≤ tol && e2 ≤ tol
            @info "The Surplus Functions Converged in $iter iterations."
            return (S = S, Smin = Smin, success = true)
        elseif iter == MaxIter
            @warn "The Surplus Functions failed to converge after $iter iterations."
            return (S = S, Smin = Smin, success = false)
        end
    end

end