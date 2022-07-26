# what is the path?
const _BASE_PATH = pwd()
const _PATH_TO_SRC = joinpath(_BASE_PATH, "src")
const _PATH_TO_DATA = joinpath(_BASE_PATH, "data")
const _PATH_TO_MODEL = joinpath(_BASE_PATH, "model")
const _PATH_TO_TMP = joinpath(_BASE_PATH, "tmp")
const _PATH_TO_SYNTHETIC_ENSEMBLE = joinpath(_BASE_PATH, "notebooks", "synthetic_10K_ensemble_s_system")
const _PATH_TO_ACTUAL_ENSEMBLE = joinpath(_BASE_PATH, "actual_ensemble_s_system")
const _PATH_TO_BASE_FIGS = joinpath(_BASE_PATH, "figs")

# load external packages -
using CSV
using DataFrames
using LinearAlgebra
using Plots
using Distributions
using DifferentialEquations
using Tables
using NumericalIntegration
using Optim
using Colors

# load my codes -
include(joinpath(_PATH_TO_SRC, "Balances.jl"))
include(joinpath(_PATH_TO_SRC, "Kinetics.jl"))
include(joinpath(_PATH_TO_SRC, "Factory.jl"))
include(joinpath(_PATH_TO_SRC, "Utility.jl"))
include(joinpath(_PATH_TO_SRC, "Simulation.jl"))
include(joinpath(_PATH_TO_SRC, "Compute.jl"))
include(joinpath(_PATH_TO_SRC, "Learn.jl"))