Copyright 2023 Bochuan Lyu, Illya V. Hicks, and Joey Huchette.
# PiecewiseLinearRelaxation.jl

Code for the paper "Building Formulations for Piecewise Linear Relaxations of Nonlinear Functions" by Bochuan Lyu, Illya V. Hicks, and Joey Huchette.

This code uses [Julia](https://julialang.org/), [JuMP](https://jump.dev/), and requires [Gurobi](https://www.gurobi.com/) and [SCIP](https://scipopt.org/) solvers.

***

## Required Packages

```julia
using Pkg
using Random
using CSV
using MAT
using LinearAlgebra
using StatsBase
using DataFrames
using SparseArrays
using Dates
using JuMP, Gurobi, SCIP
using Graphs
```

An example to install StatsBase package:

```julia
using Pkg
Pkg.add(StatsBase)
```

***

## Summary of Repository
- `src/` contains all the formulations of piecewise linear relaxation approaches that are listed in the Computational Results section of the paper.
  - `types.jl`, `jump.jl` and `PiecewiseLinearOpt.jl` are the modified code from [PiecewiseLinearOpt.jl](https://github.com/joehuchette/PiecewiseLinearOpt.jl). We implement Merged method in the paper.
  - `robotFunctions.jl` contains the model of 2D inverse kinematics problem.
  - `robotMain.jl` contains the code to solve a sample instance of 2D inverse kinematics problem.
  - `shareOfChoiceFunctions.jl` contains the model of share-of-choice problem.
  - `shareOfChoiceMain.jl` contains the code to solve a sample instance of share-of-choice problem.
- `experiments/` contains the code for computational experiments in the paper.
  - `robot` and `robot.jl`: computational results for 2D inverse kinematics problems.
  - `shareOfChoice` and `shareOfChoice.jl`: computational results for share-of-choice problems.
  - `shareOfChoiceNLP` and `shareOfChoiceNLP.jl`: computational comparison between the original share-of-choice problems solved by SCIP MINLP solver and the relaxed problems solved by Gurobi MILP solver.
- `analysis/` contains the code to summarize the computational results.

***

## Running Computational Experiments

```bash
cd experiments/
# Computational Results in Table 1
./robot
# Computational Results in Table 2
./shareOfChoice
# Computational Results in Table 3
./shareOfChoise
```

***

## Acknowledgements

`types.jl`, `jump.jl` and `PiecewiseLinearOpt.jl` are the modified code from [PiecewiseLinearOpt.jl](https://github.com/joehuchette/PiecewiseLinearOpt.jl), which is the code for the paper ["Nonconvex Piecewise Linear Functions: Advanced Formulations and Simple Modeling Tools"]("https://pubsonline.informs.org/doi/abs/10.1287/opre.2019.1973").
