# load the include -
include("Include.jl")

# how many samples?
number_of_samples = 36
number_of_parameters = 12
ensemble_archive = zeros(number_of_parameters+1,1); # first row is the fitness 

# load the training data -
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Thrombin-TF-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)[!,3:end]

# build the model structure -
path_to_model_file = joinpath(_PATH_TO_MODEL, "Coagulation.net")
model_buffer = read_model_file(path_to_model_file)

# build the default model structure -
model = build_default_model_dictionary(model_buffer)

# main loop -
p_previous = nothing
for i ∈ 1:number_of_samples

    # run the learn routine -
    (p, T, Xₘ, Yₘ, Y) = learn_optim(i, model, training_df; pₒ = nothing)

    # compute the fitness -
    fitness = norm((Yₘ .- Y).^2)
    global ensemble_archive[1] = fitness
    
    # cache the parameters -
    for k ∈ 1:number_of_parameters
        global ensemble_archive[k+1] = p[k]
    end

    # dump parameters to disk (just in case) -
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "PSET-Actual-P$(i).csv"), 
        Tables.table(ensemble_archive); header=["parameters"])
    
    # dump output to disk -
    data_output = [Y Yₘ]
    data_output_header = ["actual", "simulated"]
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "OUT-Actual-P$(i).csv"), 
        Tables.table(data_output); header = data_output_header)

    # dump state to disk -
    data_state = transpose(vcat(reshape(T,1,length(T)), Xₘ))
    data_state_header = ["T", "FII", "FVII", "FV", "FX", "FVIII", "FIX", "FXI", "FXII", "FIIa", "FVIIa", "FVa", "FXa", "FVIIIa", "FIXa", "FXIa", "FXIIa", "FIIa_inactive"]
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "SIM-Actual-P$(i).csv"),  
        Tables.table(data_state); header=data_state_header)

    # update/clean up for the next patient -
    global p_previous = p
end