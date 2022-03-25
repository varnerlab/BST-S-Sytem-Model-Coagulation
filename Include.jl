# what is the path?
_BASE_PATH = pwd()
_PATH_TO_SRC = joinpath(_BASE_PATH, "src")
_PATH_TO_DATA = joinpath(_BASE_PATH, "data")
_PATH_TO_MODEL = joinpath(_BASE_PATH, "model")
_PATH_TO_TMP = joinpath(_BASE_PATH, "tmp")

# load external packages -
using CSV
using DataFrames
using LinearAlgebra
using Plots
using DataFrames
using Distributions
using DifferentialEquations
using Tables
using NumericalIntegration
using Optim

# load my codes -
include(joinpath(_PATH_TO_SRC, "Balances.jl"))
include(joinpath(_PATH_TO_SRC, "Kinetics.jl"))
include(joinpath(_PATH_TO_SRC, "Factory.jl"))
include(joinpath(_PATH_TO_SRC, "Utility.jl"))
include(joinpath(_PATH_TO_SRC, "Simulation.jl"))
include(joinpath(_PATH_TO_SRC, "Compute.jl"))
include(joinpath(_PATH_TO_SRC, "Learn.jl"))