# load the include -
include("Include.jl")

# how many samples?
number_of_samples = 25
number_of_parameters = 12
ensemble_archive = zeros(number_of_parameters+1,number_of_samples); # first row is the fitness 

# main loop -
p_previous = nothing
for i ∈ 1:number_of_samples

    # run the learn routine -
    (p, Yₘ, Y) = learn_optim(i; pₒ = p_previous)

    # compute the fitness -
    fitness = norm((Yₘ .- Y).^2)
    global ensemble_archive[1,i] = fitness
    
    # cache the parameters -
    for k ∈ 1:number_of_parameters
        global ensemble_archive[k+1,i] = p[k]
    end

    global p_previous = p

    
    # dump parameters to disk (just in case) -
    path_to_synthetic_ensemble = joinpath(_PATH_TO_TMP, "Ensemble-Synthetic-Thrombin-TF-1K-Optim.csv")
    CSV.write(path_to_synthetic_ensemble, Tables.table(ensemble_archive))
    
    # dump output to disk -
    data_output = [Y Yₘ]
    path_to_synthetic_patient_otuput = joinpath(_PATH_TO_TMP, "SIM-Synthetic-P$(i).csv")
    CSV.write(path_to_synthetic_patient_otuput, Tables.table(ensemble_archive))
end