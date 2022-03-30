# include -
include("Include.jl")

# This script requires Flux.jl (we manually include here, its takes a while to load)
using Flux
using Flux: @epochs

# how many training data sets do we have?
number_of_synthetic_patients = 500

trajectory_range = range(1,stop=2001,step=100) |> collect
N = length(trajectory_range)

# initialize storage for the training data -
training_data = Vector{Tuple{Vector{Float32},Vector{Float32}}}()
for i ∈ 1:number_of_synthetic_patients

    # load the output -
    output_df = CSV.read(joinpath(_PATH_TO_TMP,"OUT-Synthetic-1k-P$(i).csv"), DataFrame)
    output_vector = convert.(Float32, output_df[!,:simulated])

    # load the FIIa -
    simulation_df = CSV.read(joinpath(_PATH_TO_TMP,"SIM-Synthetic-1k-P$(i).csv"), DataFrame)
    FIIa_vector = convert.(Float32, simulation_df[trajectory_range,:FIIa])

    # make a training tuple -
    training_tuple = (
        output_vector,FIIa_vector
    );

    # insert -
    push!(training_data, training_tuple)
end

# build the model -
FIIa_trajectory_model = Chain(Dense(5, 10, σ), Dense(10, N));

# setup a loss function -
loss(x, y) = Flux.Losses.mae(FIIa_trajectory_model(x), y; agg = mean)

# pointer to params -
ps = Flux.params(FIIa_trajectory_model)

# # use old school gradient descent -
opt = Momentum(0.1, 0.95)

@epochs 1000 Flux.train!(loss, ps, training_data, opt)
