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
include("robotFunctions.jl")



# arms = [10.0, 5.0, 4.0, 3.0, 2.0, 1.0]
# angleRanges = [(-0.5 * pi, 0.5 * pi), (-0.1 * pi, 0.1 * pi), (-0.1 * pi, 0.1 * pi),
#                (-0.1 * pi, 0.1 * pi), (-0.1 * pi, 0.1 * pi), (-0.1 * pi, 0.1 * pi)]
arms = [5.0, 5.0, 3.0, 3.0]
angleRanges = [(-0.5 * pi, 0.5 * pi), (-0.25 * pi, 0.25 * pi), (-0.25 * pi, 0.25 * pi),
                (-0.25 * pi, 0.25 * pi)]



# arms = [10.0, 5.0]
# angleRanges = [(-0.5 * pi, 0.5 * pi), (-0.1 * pi, 0.1 * pi)]
targetPosition = (3.4, 5.2)# (3.5, 6.0)
targetAngle = pi/2
beta = 0.1
points = 500
# :pwl, :pwr, :pwlMerge
method = :pwlMerge
initAngle = pi/2
times=2
# PWL: :LogarithmicIB, :ZigZag, :ZigZagInteger, :GeneralizedCelaya,
# :SOS2, :Incremental, :CC, :DisaggLogarithmic, :MC, :DLog
# PWR: :general, :BRGC, :balancedCode, :biclique
pwl_method = :ZigZagInteger
timeLimit = 60.0

println("The number of points:", points, " The method: ", method, " the PWL/PWR Approach: ", pwl_method)

# thetas = inverseKinematics2DNLP(arms, angleRanges, targetPosition, targetAngle, initAngle=initAngle,
#         beta=beta, timeLimit=timeLimit)

m, thetas = inverseKinematics2D(arms, angleRanges, targetPosition, targetAngle, initAngle=initAngle,
        beta=beta, points=points, method=method, pwl_method=pwl_method, times=times, timeLimit=timeLimit)

locations = forward2D(thetas, arms, initAngle)
