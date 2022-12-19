__precompile__()

module PiecewiseLinearOpt

using JuMP
import MathOptInterface
const MOI = MathOptInterface
using LinearAlgebra
using Random

export PWLFunction, UnivariatePWLFunction, BivariatePWLFunction, piecewiselinear, piecewiselinear2

include("types.jl")
include("jump.jl")

end # module
