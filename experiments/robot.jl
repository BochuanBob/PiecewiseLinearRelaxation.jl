using Pkg
using Random
using CSV
using StatsBase
using DataFrames
using SparseArrays
using JuMP, Gurobi
using Graphs

const VarOrAff = Union{JuMP.VariableRef, JuMP.AffExpr}
const VarOrAffArr = Union{Array{JuMP.VariableRef}, Array{JuMP.AffExpr}}
const VarOrAff1DArr = Union{Vector{JuMP.VariableRef}, Vector{JuMP.AffExpr}}

include("../src/PiecewiseLinearOpt.jl")
using Main.PiecewiseLinearOpt

include("../src/relaxation.jl")
include("../src/robotFunctions.jl")


function main()
    # Make results folder.

    arms = [5.0, 5.0, 3.0, 1.0]
    angleRanges = [(-0.5 * pi, 0.5 * pi), (-0.25 * pi, 0.25 * pi), (-0.25 * pi, 0.25 * pi),
                    (-0.25 * pi, 0.25 * pi)]
    initAngle = pi/2
    beta = 0.1
    points = parse(Int64, ARGS[1])
    times = parse(Int64, ARGS[2])
    i = parse(Int64, ARGS[3])
    timeLimit = parse(Float64, ARGS[4])
    
    j = 1
    isdir("../results") || mkdir("../results")
    folderName = String(ARGS[5])
    isdir(folderName) || mkdir(folderName)
    outputPath = string(folderName, String("result_$(points)_$(times)/"))
    isdir(outputPath) || mkdir(outputPath)


    # :pwl, :pwr, :pwlMerge
    methodList = [:pwl, :pwlMerge, :pwr]
    pwlMethods = Dict(:pwl=>[:Incremental, :CC, :MC, :DLog, :Logarithmic, :ZigZag, :ZigZagInteger],
                      :pwlMerge=>[:Incremental, :DLog, :Logarithmic, :SOS2, :ZigZag, :ZigZagInteger],
                      :pwr=>[:Incremental, :DLog, :BRGC, :balancedCode, :biclique])


    seed = 2022 + i
    Methods, Approachs,
            PointList, TimesList,
            Indices, ThetaList,
            Targets, Actuals,
            Bounds, GurobiTimes = [], [], [], [], [], [], [], [], [], []
    for method in methodList
        for pwl_method in pwlMethods[method]
            Random.seed!(seed)
            targetPosition = (Random.rand(Float64)*10, Random.rand(Float64)*10)
            targetAngle = Random.rand(Float64) * 2 * pi + initAngle - pi
            println(" The method: ", method, " the PWL/PWR Approach: ", pwl_method)
            println("The number of points: ", points,
                    "The number of extra points in each relaxation:", times)
            println("The samples: ", i)
            m, thetas = inverseKinematics2D(arms, angleRanges, targetPosition, targetAngle, initAngle=initAngle,
                    beta=beta, points=points, method=method, pwl_method=pwl_method, times=times, timeLimit=timeLimit)
            locations, actualAng = forward2D(thetas, arms, initAngle)
            actualAng += initAngle
            lLen = length(locations)
            actualPos = locations[lLen]

            push!(Methods, method)
            push!(Approachs, pwl_method)
            append!(Indices, i)
            append!(PointList, points)
            append!(TimesList, times)
            push!(ThetaList, thetas)
            push!(Targets, [targetPosition[1], targetPosition[2], targetAngle])
            push!(Actuals, [actualPos[1], actualPos[2], actualAng])
            push!(Bounds, MOI.get(m, Gurobi.ModelAttribute("ObjBound")))
            append!(GurobiTimes, MOI.get(m, Gurobi.ModelAttribute("Runtime")))
        end
    end
    df = DataFrame(Indices=Indices, Pieces=PointList, ExtraPoints=TimesList,
                    Methods=Methods, Approaches=Approachs,
                    Thetas=ThetaList, Targets=Targets, Actuals=Actuals,
                    Bounds=Bounds, GurobiTimes=GurobiTimes)
    CSV.write(string(outputPath, "results_$(points)_$(times)_$(i)", ".csv"), df)

end

main()
