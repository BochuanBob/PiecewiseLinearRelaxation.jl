using Pkg
using Random
using CSV
using MAT
using DataFrames
using SparseArrays

using JuMP, Gurobi, SCIP
using Graphs

const VarOrAff = Union{JuMP.VariableRef, JuMP.AffExpr}
const VarOrAffArr = Union{Array{JuMP.VariableRef}, Array{JuMP.AffExpr}}
const VarOrAff1DArr = Union{Vector{JuMP.VariableRef}, Vector{JuMP.AffExpr}}

include("PiecewiseLinearOpt.jl")
using Main.PiecewiseLinearOpt

include("relaxation.jl")
include("shareOfChoiceFunctions.jl")

v, S, eta, C = 10, 6, 15, 0.2
Random.seed!(2022)
betas = rand(Float64, (v, S, eta)) * 10 .- 5.0
Us = rand(Float64, v) * 10
lambdas = rand(Float64, v)

# :pwl, :pwr, :pwlMerge
method = :pwr
points = 50
times=4
# PWL: :LogarithmicIB, :ZigZag, :ZigZagInteger, :GeneralizedCelaya,
# :SOS2, :Incremental, :CC, :DisaggLogarithmic, :MC, :DLog
# PWR: :general, :BRGC, :balancedCode, :biclique
pwl_method = :biclique
timeLimit = 1200.0

println("The number of points:", points, " The method: ", method, " the PWL/PWR Approach: ", pwl_method)

shareOfChoise(lambdas, betas, Us, C, points=points, method=method, times=times,
        pwl_method=pwl_method, timeLimit=timeLimit)
