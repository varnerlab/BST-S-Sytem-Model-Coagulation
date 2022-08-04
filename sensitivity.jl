# Include the include -
include("Include.jl")

function performance(κ, model::Dict{String,Any}, training_df::DataFrame, index::Int64)

    # main simulation logic -
    SF = 1e9

    # setup static -
    sfa = model["static_factors_array"]
    sfa[1] = training_df[index,:TFPI]       # 1 TFPI
    sfa[2] = training_df[index,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF                   # 3 TF
    sfa[6] = 0.005                          # 6 TRAUMA

    # grab the multiplier from the data -
    ℳ = model["number_of_dynamic_states"]
    xₒ = zeros(ℳ)
    xₒ[1] = training_df[index, :II]         # 1 FII 
    xₒ[2] = training_df[index, :VII]        # 2 FVII 
    xₒ[3] = training_df[index, :V]          # 3 FV
    xₒ[4] = training_df[index, :X]          # 4 FX
    xₒ[5] = training_df[index, :VIII]       # 5 FVIII
    xₒ[6] = training_df[index, :IX]         # 6 FIX
    xₒ[7] = training_df[index, :XI]         # 7 FXI
    xₒ[8] = training_df[index, :XII]        # 8 FXII 
    xₒ[9] = (1e-14)*SF                      # 9 FIIa
    model["initial_condition_vector"] = xₒ

    # get the G matrix -
    G = model["G"]

    # κ map
    # 1 - 9 : α vector
    model["α"] = κ[1:9]
    idx = indexin(model, "FVIIa")   # 10
    G[idx, 4] = κ[10]
    idx = indexin(model, "AT")      # 11
    G[idx, 9] = κ[11]
    idx = indexin(model, "TFPI")    # 12
    G[idx, 1] = -1*κ[12]

    # run the model -
    (T,U) = evaluate(model)
    Xₘ = hcat(U...)

    # compute the model ouput -
    Yₘ = model_output_vector(T,Xₘ[9,:]) # properties of the Thrombin curve 

    # return -
    return Yₘ 
end

# load the training data -
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Thrombin-TF-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)[!,3:end]

# build the model structure -
path_to_model_file = joinpath(_PATH_TO_MODEL, "Coagulation.net")
model_buffer = read_model_file(path_to_model_file)

# build the default model structure -
model = build_default_model_dictionary(model_buffer)

# setup sensitivity analysis -
# load a pset -
pset_filename = "PSET-Actual-P17.csv"
pset_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, pset_filename), DataFrame)

# create lower and upper bound array -
NP = length(pset_df[!,:parameters])
L = zeros(NP-1)
U = zeros(NP-1)
for pᵢ ∈ 1:(NP - 1)
    L[pᵢ] = 0.1*pset_df[pᵢ + 1,:parameters]
    U[pᵢ] = 10.0*pset_df[pᵢ + 1,:parameters]
end

# setup call to Morris method -
F(κ) =  performance(κ, model, training_df, 17)
m = gsa(F, Morris(num_trajectory=10000), [[L[i],U[i]] for i in 1:(NP-1)]);

# dump sensitivity data to disk -
jldsave("Sensitivity-Morris-P36-V2-N10K.jld2"; mean=m.means, variance=m.variances)