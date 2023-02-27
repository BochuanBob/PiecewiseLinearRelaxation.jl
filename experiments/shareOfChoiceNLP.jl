using Pkg
using Random
using CSV
using StatsBase
using DataFrames
using SparseArrays
using JuMP, Gurobi, SCIP
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

    samples = parse(Int64, ARGS[1])
    timeLimit = parse(Float64, ARGS[2])

    j = 1
    isdir("../results") || mkdir("../results")
    folderName = String(ARGS[3])
    isdir(folderName) || mkdir(folderName)
    outputPath = string(folderName, String("result/"))
    isdir(outputPath) || mkdir(outputPath)

    Formulations = ["NLP", "Gurobi_tiny", "Gurobi_small", "Gurobi_large"]
    Params = Dict("SCIP_small"=>["scip", 10, 2, :pwlMerge, :Incremental],
                "Gurobi_tiny"=>["gurobi", 10, 1, :pwlMerge, :Incremental],
                "Gurobi_small"=>["gurobi", 10, 2, :pwlMerge, :Incremental],
                "Gurobi_large"=>["gurobi", 50, 2, :pwr, :balancedCode])
    for i in 1:samples
        seed = 2022 + i
        Indices, MethodList, ParaList, Objs,
                Bounds, SolverTimes = [], [], [], [], [], []

        Random.seed!(seed)
        betas = rand(Float64, (v, S, eta)) * 10 .- 5.0
        Us = rand(Float64, v) * 10
        lambdas = rand(Float64, v)
        lambdas = lambdas / sum(lambdas)
        println("The samples: ", i)
        for formulation in Formulations
            if (formulation == "NLP")
                m = shareOfChoiceNLP(lambdas, betas, Us, C, timeLimit=timeLimit, threads=1)
            else
                solver, points, times, method, pwl_method = Params[formulation]
                m = shareOfChoice(lambdas, betas, Us, C, points=points, method=method, times=times,
                        pwl_method=pwl_method, timeLimit=timeLimit, solver=solver, threads=1)
            end
            optimizer = backend(m)
            append!(Indices, i)
            push!(ParaList, [v, S, eta, C])
            if MOI.get(optimizer, MOI.ResultCount()) > 0
                append!(Objs, MOI.get(optimizer, MOI.ObjectiveValue()))
            else
                append!(Objs, -Inf)
            end
            push!(Bounds, MOI.get(optimizer, MOI.ObjectiveBound()))
            push!(MethodList, formulation)
            append!(SolverTimes, MOI.get(optimizer, MOI.SolveTimeSec()))
        end
        df = DataFrame(Indices=Indices, Methods=MethodList,
                        Parameters=ParaList, Objectives=Objs,
                        Bounds=Bounds, SolverTimes=SolverTimes)
        CSV.write(string(outputPath, "results_$(i)", ".csv"), df)
    end
end

main()
