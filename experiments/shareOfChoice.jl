using Pkg
using Random
using CSV
using StatsBase
using DataFrames
using SparseArrays
using Dates
using JuMP, Gurobi
using Graphs

const VarOrAff = Union{JuMP.VariableRef, JuMP.AffExpr}
const VarOrAffArr = Union{Array{JuMP.VariableRef}, Array{JuMP.AffExpr}}
const VarOrAff1DArr = Union{Vector{JuMP.VariableRef}, Vector{JuMP.AffExpr}}

include("../src/PiecewiseLinearOpt.jl")
using Main.PiecewiseLinearOpt

include("../src/relaxation.jl")
include("../src/shareOfChoiceFunctions.jl")


function main()
    # Make results folder.

    v, S, eta, C = 10, 6, 15, 0.2

    points = parse(Int64, ARGS[1])
    times = parse(Int64, ARGS[2])
    samples = parse(Int64, ARGS[3])
    timeLimit = parse(Float64, ARGS[4])

    todayStr = Dates.format(Dates.today(), "yyyy_mm_dd")
    j = 1
    isdir("../results") || mkdir("../results")
    folderName = String(ARGS[5])
    isdir(folderName) || mkdir(folderName)
    outputPath = string(folderName, String("result_$(points)_$(times)/"))
    isdir(outputPath) || mkdir(outputPath)


    # :pwl, :pwr, :pwlMerge
    methodList = [:pwlMerge, :pwr]
    pwlMethods = Dict(:pwlMerge=>[:Incremental, :Logarithmic, :SOS2, :ZigZag, :ZigZagInteger],
                      :pwr=>[:BRGC, :balancedCode, :biclique])

    for i in 1:samples
        seed = 2022 + i
        Methods, Approachs,
                PointList, TimesList,
                Indices, ParaList, Objs,
                Bounds, GurobiTimes = [], [], [], [], [], [], [], [], []

        Random.seed!(seed)
        betas = rand(Float64, (v, S, eta)) * 10 .- 5.0
        Us = rand(Float64, v) * 10
        lambdas = rand(Float64, v)
        for method in methodList
            for pwl_method in pwlMethods[method]
                println(" The method: ", method, " the PWL/PWR Approach: ", pwl_method)
                println("The number of points: ", points,
                        "The number of extra points in each relaxation:", times)
                println("The samples: ", i)

                m = shareOfChoise(lambdas, betas, Us, C, points=points, method=method, times=times,
                        pwl_method=pwl_method, timeLimit=timeLimit)

                push!(Methods, method)
                push!(Approachs, pwl_method)
                append!(Indices, i)
                append!(PointList, points)
                append!(TimesList, times)
                push!(ParaList, [v, S, eta, C])
                if (MOI.get(m, Gurobi.ModelAttribute("SolCount")) > 0)
                    append!(Objs, MOI.get(m, Gurobi.ModelAttribute("ObjVal")))
                else
                    append!(Objs, Inf)
                end
                push!(Bounds, MOI.get(m, Gurobi.ModelAttribute("ObjBound")))
                append!(GurobiTimes, MOI.get(m, Gurobi.ModelAttribute("Runtime")))
            end
        end
        df = DataFrame(Indices=Indices, Pieces=PointList, ExtraPoints=TimesList,
                        Methods=Methods, Approaches=Approachs,
                        Parameters=ParaList, Objectives=Objs,
                        Bounds=Bounds, GurobiTimes=GurobiTimes)
        CSV.write(string(outputPath, "results_$(points)_$(times)_$(i)", ".csv"), df)
    end
end

main()
