# load the include -
include("Include.jl")

# how many samples?
number_of_samples = 25
number_of_parameters = 12
ensemble_archive = zeros(number_of_parameters+1,1); # first row is the fitness 

# load the training data -
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Synthetic-Thrombin-TF-1K.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

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
    path_to_synthetic_ensemble = joinpath(_PATH_TO_TMP, "PSET-Synthetic-1k-P$(i).csv")
    CSV.write(path_to_synthetic_ensemble, Tables.table(ensemble_archive))
    
    # dump output to disk -
    data_output = [Y Yₘ]
    path_to_synthetic_patient_otuput = joinpath(_PATH_TO_TMP, "OUT-Synthetic-1k-P$(i).csv")
    CSV.write(path_to_synthetic_patient_otuput, Tables.table(data_output))

    # dump state to disk -
    data_state = vcat(reshape(T,1,length(T)), Xₘ)
    path_to_synthetic_patient_simulation = joinpath(_PATH_TO_TMP, "SIM-Synthetic-1k-P$(i).csv")
    CSV.write(path_to_synthetic_patient_simulation, Tables.table(data_state))

    # update/clean up for the next patient -
    global p_previous = p
end