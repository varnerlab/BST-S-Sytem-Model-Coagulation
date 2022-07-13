include("Include.jl")

# This script requires Flux.jl (we manually include here, its takes a while to load)
using Flux
using Flux: @epochs
using BSON: @load

# load my model -
@load "FIIa-trajectory.bson" FIIa_trajectory_model

# how many training data sets do we have?
synthetic_patient_index = 10000

trajectory_range = range(1,stop=2001,step=10) |> collect
N = length(trajectory_range)

# load the FIIa -
simulation_df = CSV.read(joinpath(_PATH_TO_SYNTHETIC_ENSEMBLE,"SIM-Synthetic-1k-P$(synthetic_patient_index).csv"), DataFrame)
FIIa_vector = convert.(Float32, simulation_df[!,:FIIa])
T = convert.(Float32, simulation_df[!,:T])

# load the output -> (input to the FBB)
output_df = CSV.read(joinpath(_PATH_TO_SYNTHETIC_ENSEMBLE,"OUT-Synthetic-1k-P$(synthetic_patient_index).csv"), DataFrame)
output_vector = convert.(Float32, output_df[!,:simulated])

# compute the trajectory -
FIIa_estimated = FIIa_trajectory_model(output_vector)
T_estimated = convert.(Float32, simulation_df[trajectory_range,:T])


plot(T,FIIa_vector,c=:red,lw=2, label="BST model (P(i) = $(synthetic_patient_index))")
plot!(T_estimated, FIIa_estimated,c=:black,lw=2, label="ANN model (P(i) = $(synthetic_patient_index))")
xlabel!("Time [min]", fontsize=18)
ylabel!("FIIa (nM)", fontsize=18)
current()

# dump figures to disk -
fig_name = "FIIa-BST-ANN-MAE-E500-P$(synthetic_patient_index).pdf"
savefig(joinpath(_PATH_TO_TMP, fig_name))
