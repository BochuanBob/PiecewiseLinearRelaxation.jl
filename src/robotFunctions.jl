function calculateRanges(angleRanges::Vector{Tuple{Float64, Float64}})
    lower, upper = angleRanges[1]
    res = Vector{Tuple{Float64, Float64}}([(lower, upper)])
    for i in 2:length(angleRanges)
        lower += angleRanges[i][1]
        upper += angleRanges[i][2]
        append!(res, [(lower, upper)])
    end
    return res
end

function forward2D(thetas::Vector{Float64}, arms::Vector{Float64}, initAngle::Float64)
    locations = [[0.0, 0.0]]
    angle = 0.0
    unitVec = [cos(initAngle) -sin(initAngle);
                 sin(initAngle) cos(initAngle)] * [1.0, 0.0]
    for i in 1:length(thetas)
        angle += thetas[i]
        loc = locations[length(locations)]
        newLoc = [cos(angle) -sin(angle);
                     sin(angle) cos(angle)] * arms[i] * unitVec + loc
        append!(locations, [newLoc])
    end
    return locations, angle
end

function inverseKinematics2DNLP(arms::Vector{Float64},
        angleRanges::Vector{Tuple{Float64, Float64}},
        targetPosition::Tuple{Float64, Float64}, targetAngle::Float64;
        beta::Float64=0.0, timeLimit::Float64=60.0, initAngle::Float64=pi/2)
    if (length(angleRanges) == 0)
        error("The number of angles should be at least one.")
    end
    if (length(arms) != length(angleRanges))
        error("The number of arms is different from the number of angles.")
    end
    aLen = length(arms)
    unitVec = [cos(initAngle) -sin(initAngle);
                 sin(initAngle) cos(initAngle)] * [1.0, 0.0]
    cumAngleRanges = calculateRanges(angleRanges)
    m = direct_model(SCIP.Optimizer())
    theta = @variable(m, [1:aLen])

    @constraint(m, [i=1:aLen], theta[i] >= angleRanges[i][1])
    @constraint(m, [i=1:aLen], theta[i] <= angleRanges[i][2])

    cumTheta = [sum(theta[i] for i in 1:j) for j in 1:aLen]
    xPos, yPos = 0.0, 0.0
    sinVars = Vector{VariableRef}(undef, aLen)
    cosVars = Vector{VariableRef}(undef, aLen)
    for i in 1:aLen
        sinVars[i] = @variable(m)
        @NLconstraint(m, sinVars[i] == sin(cumTheta[i]))

        cosVars[i] = @variable(m)
        @NLconstraint(m, cosVars[i] == cos(cumTheta[i]))

        armVec = arms[i] * unitVec
        xPos += cosVars[i] * armVec[1] - sinVars[i] * armVec[2]
        yPos += sinVars[i] * armVec[1] + cosVars[i] * armVec[2]
        # xPos += cosVar * arms[i]
        # yPos += sinVar * arms[i]
    end
    absDiff = @variable(m, [1:3], lower_bound=0.0)
    @constraint(m, absDiff .>= [targetPosition[1], targetPosition[2], targetAngle]
                            .- [xPos, yPos, cumTheta[aLen] + initAngle])
    @constraint(m, absDiff .>= [xPos, yPos, cumTheta[aLen] + initAngle]
                    .- [targetPosition[1], targetPosition[2], targetAngle])
    @objective(m, Min, absDiff[1] + absDiff[2] + beta * absDiff[3])
    optimize!(m)
    println("Position x: ", value(xPos), "Position y: ", value(yPos))
    println("Final Angle: ", value(cumTheta[aLen]))
    for i in 1:aLen
        println("Actual Sin: ", sin(value(cumTheta[i])) , " MIP Sin: ", value(sinVars[i]))
        println("Actual Cos: ", cos(value(cumTheta[i])) , " MIP Cos: ", value(cosVars[i]))
    end
    return m, value.(theta)
end

function inverseKinematics2D(arms::Vector{Float64},
        angleRanges::Vector{Tuple{Float64, Float64}},
        targetPosition::Tuple{Float64, Float64}, targetAngle::Float64;
        beta::Float64=0.0, points::Int64=9, method=:pwl, times::Int64=1,
        pwl_method=:Logarithmic, timeLimit::Float64=60.0, initAngle::Float64=pi/2)
    if (length(angleRanges) == 0)
        error("The number of angles should be at least one.")
    end
    if (length(arms) != length(angleRanges))
        error("The number of arms is different from the number of angles.")
    end
    aLen = length(arms)
    unitVec = [cos(initAngle) -sin(initAngle);
                 sin(initAngle) cos(initAngle)] * [1.0, 0.0]
    cumAngleRanges = calculateRanges(angleRanges)
    m = direct_model(Gurobi.Optimizer())
    set_optimizer_attribute(m, "OutputFlag", 1)
    # set_optimizer_attribute(m, "PreCrush", 1)
    set_optimizer_attribute(m, "Threads", 4)
    set_optimizer_attribute(m, "TimeLimit", timeLimit)
    theta = @variable(m, [1:aLen])

    @constraint(m, [i=1:aLen], theta[i] >= angleRanges[i][1])
    @constraint(m, [i=1:aLen], theta[i] <= angleRanges[i][2])

    cumTheta = [sum(theta[i] for i in 1:j) for j in 1:aLen]
    xPos, yPos = 0.0, 0.0
    sinVars = Vector{VariableRef}(undef, aLen)
    cosVars = Vector{VariableRef}(undef, aLen)
    for i in 1:aLen
        sinVars[i] = @variable(m)
        sinVar = sinVars[i]
        breaks = getSinBreaks(cumAngleRanges[i][1], cumAngleRanges[i][2])
        # println(breaks)
        nonlinearRelaxation!(m, cumTheta[i], sinVar,
                        cumAngleRanges[i][1], cumAngleRanges[i][2], sin, dsin,
                        points, breaks, method=method,
                        times=times, pwl_method=pwl_method)
        cosVars[i] = @variable(m)
        cosVar = cosVars[i]
        breaks = getCosBreaks(cumAngleRanges[i][1], cumAngleRanges[i][2])
        # println(breaks)
        nonlinearRelaxation!(m, cumTheta[i], cosVar,
                        cumAngleRanges[i][1], cumAngleRanges[i][2], cos, dcos,
                        points, breaks, method=method,
                        times=times, pwl_method=pwl_method)
        armVec = arms[i] * unitVec
        xPos += cosVar * armVec[1] - sinVar * armVec[2]
        yPos += sinVar * armVec[1] + cosVar * armVec[2]
        # xPos += cosVar * arms[i]
        # yPos += sinVar * arms[i]
    end
    absDiff = @variable(m, [1:3], lower_bound=0.0)
    @constraint(m, absDiff .>= [targetPosition[1], targetPosition[2], targetAngle]
                            .- [xPos, yPos, cumTheta[aLen] + initAngle])
    @constraint(m, absDiff .>= [xPos, yPos, cumTheta[aLen] + initAngle]
                    .- [targetPosition[1], targetPosition[2], targetAngle])
    @objective(m, Min, absDiff[1] + absDiff[2] + beta * absDiff[3])
    optimize!(m)
    println("Position x: ", value(xPos), "Position y: ", value(yPos))
    println("Final Angle: ", value(cumTheta[aLen]))
    for i in 1:aLen
        println("Actual Sin: ", sin(value(cumTheta[i])) , " MIP Sin: ", value(sinVars[i]))
        println("Actual Cos: ", cos(value(cumTheta[i])) , " MIP Cos: ", value(cosVars[i]))
    end
    return m, value.(theta)
end

function getSinBreaks(lower::Float64, upper::Float64)
    res = Vector{Float64}([])
    l = ceil(Int64, lower / pi)
    u = floor(Int64, upper / pi)
    for i in l:u
        append!(res, i * pi)
    end
    return res
end

function getCosBreaks(lower::Float64, upper::Float64)
    res = Vector{Float64}([])
    l = ceil(Int64, lower / pi - 0.5)
    u = floor(Int64, upper / pi - 0.5)
    for i in l:u
        append!(res, (i + 0.5) * pi)
    end
    return res
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
