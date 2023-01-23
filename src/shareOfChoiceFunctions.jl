function shareOfChoise(lambdas::Vector{Float64}, betas::Array{Float64, 3},
        Us::Vector{Float64}, C::Float64;
        points::Int64=9, method=:pwl, times::Int64=1,
        pwl_method=:Logarithmic, timeLimit::Float64=60.0)
    if (length(lambdas) != length(Us))
        error("The size of lambda is different from the size of U.")
    end
    v, S, eta = size(betas)
    m = direct_model(Gurobi.Optimizer())
    set_optimizer_attribute(m, "OutputFlag", 1)
    # set_optimizer_attribute(m, "PreCrush", 1)
    set_optimizer_attribute(m, "Threads", 4)
    set_optimizer_attribute(m, "TimeLimit", timeLimit)

    x = @variable(m, [1:eta], lower_bound=0.0, upper_bound=1.0)
    muBar = @variable(m, [1:v])
    mu = @variable(m, [1:v, 1:S])
    pBar = @variable(m, [1:v])
    p = @variable(m, [1:v, 1:S])
    xLowB = [0.0 for i in 1:eta]
    xUppB = [1.0 for i in 1:eta]
    for i in 1:v
        act = x -> 1 / (1 + exp(Us[i] - x))
        dAct = x -> exp(Us[i] - x) / (1 + exp(Us[i] - x))^2
        lowR, uppR = calculateRanges(
                    sum((1/S) * betas[i, s, :] for i in 1:v for s in 1:S),
                    xLowB, xUppB)
        nonlinearRelaxation!(m, muBar[i], pBar[i],
                        lowR, uppR, act, dAct,
                        points, [Us[i]], method=method,
                        times=times, pwl_method=pwl_method)
        for s in 1:S
            lowR, uppR = calculateRanges(betas[i, s, :], xLowB, xUppB)
            nonlinearRelaxation!(m, mu[i, s], p[i, s],
                            lowR, uppR, act, dAct,
                            points, [Us[i]], method=method,
                            times=times, pwl_method=pwl_method)
        end
    end
    @constraint(m, [i=1:v], muBar[i] == (1 / S) * sum(betas[i, s, :]' * x for s in 1:S))
    @constraint(m, [i=1:v, s=1:S], mu[i, s] == betas[i, s, :]' * x)
    @constraint(m, [s=1:S], sum(lambdas[i] * p[i, s] for i in 1:v) >=
                            C * sum(lambdas[i] * pBar[i] for i in 1:v))
    @objective(m, Max, sum(lambdas .* pBar))
    optimize!(m)
    return m
end

function removeClosePoints(xSort::Vector{Float64}; threshold=1e-8)
    if (length(xSort) == 0)
        return xSort
    end
    res = Vector{Float64}([xSort[1]])
    for i in 2:length(xSort)
        if abs(xSort[i] - res[length(res)]) > threshold
            append!(res, xSort[i])
        end
    end
    return res
end

function calculateRanges(coef::Vector{Float64}, lower::Vector{Float64}, upper::Vector{Float64})
    low, upp = 0.0, 0.0
    cLen = length(coef)
    for i in 1:cLen
        if (coef[i] >= 0.0)
            low += coef[i] * lower[i]
            upp += coef[i] * upper[i]
        else
            low += coef[i] * upper[i]
            upp += coef[i] * lower[i]
        end
    end
    return low, upp
end
