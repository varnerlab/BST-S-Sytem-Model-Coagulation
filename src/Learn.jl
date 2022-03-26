function learn_optim(index::Int, model::Dict{String,Any}, training_df::DataFrame; 
    pₒ::Union{Nothing,Array{Float64,1}} = nothing)

    # main simulation loop -
    SF = 1e9

    # setup static -
    sfa = model["static_factors_array"]
    sfa[1] = training_df[index,:TFPI]       # 1 TFPI
    sfa[2] = training_df[index,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF               # 3 TF
    sfa[6] = 0.005                      # 6 TRAUMA

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

    # what is the output array?
    Y = Array{Float64,1}(undef,5)
    Y[1] = training_df[index, :Lagtime]
    Y[2] = training_df[index, :Peak]
    Y[3] = training_df[index, Symbol("T.Peak")]
    Y[4] = training_df[index, :Max]
    Y[5] = training_df[index, :AUC]

    # setup initial parameter values and bounds array -
    κ = [
            
            # default: hand fit set -
            0.061   0.01 10.0   ; # 1
            1.0     0.01 10.0   ; # 2
            1.0     0.01 10.0   ; # 3
            1.0     0.01 10.0   ; # 4
            1.0     0.01 10.0   ; # 5
            1.0     0.01 10.0   ; # 6
            1.0     0.01 10.0   ; # 7
            1.0     0.01 10.0   ; # 8
            0.70    0.01 10.0   ; # 9
            0.1     0.01 10.0   ; # 10
            0.045   0.01 10.0   ; # 11
            0.065   0.01 10.0   ; # 12
        ];

    # set default set as the start -
    if (isnothing(pₒ) == true)
        pₒ = κ[:,1]
    end

    # setup the objective function -
    inner_optimizer = LBFGS()
    obj_function(p) =  loss_scalar(p, Y, model)
    results = optimize(obj_function, κ[:,2], κ[:,3], pₒ, Fminbox(inner_optimizer), 
        Optim.Options(time_limit = 600, show_trace = true, trace_simplex = true, show_every = 10))
    
    # grab the best parameters -
    p_best = Optim.minimizer(results)
    
    # run the sim w/the best parameters -
    # 1 - 9 : α vector
    model["α"] = p_best[1:9]
    G = model["G"]
    idx = indexin(model, "FVIIa")   # 10
    G[idx, 4] = p_best[10]
    idx = indexin(model, "AT")      # 11
    G[idx, 9] = p_best[11]
    idx = indexin(model, "TFPI")    # 12
    G[idx, 1] = -1*p_best[12]

    # run the model -
    (T,U) = evaluate(model)
    Xₘ = hcat(U...)
    Yₘ = model_output_vector(T, Xₘ[9,:]) # properties of the Thrombin curve 
    
    return (p_best, Yₘ, Y)
end

function learn(index::Int; pₒ::Union{Nothing,Array{Float64,1}} = nothing)

    # let'e experiment with a simple evolutionary algorithm -
    NTRIALS = 12500
    σ = 0.05

    # load the training data -
    path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Synthetic-Thrombin-TF-1K.csv")
    training_df = CSV.read(path_to_training_data, DataFrame)

    # build the model structure -
    path_to_model_file = joinpath(_PATH_TO_MODEL, "Coagulation.net")
    model_buffer = read_model_file(path_to_model_file)

    # build the default model structure -
    model = build_default_model_dictionary(model_buffer)

    # main simulation loop -
    SF = 1e9

    # setup static -
    sfa = model["static_factors_array"]
    sfa[1] = training_df[index,:TFPI]       # 1 TFPI
    sfa[2] = training_df[index,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF               # 3 TF
    sfa[6] = 0.005                      # 6 TRAUMA

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

    # what is the output array?
    Y = Array{Float64,1}(undef,5)
    Y[1] = training_df[index, :Lagtime]
    Y[2] = training_df[index, :Peak]
    Y[3] = training_df[index, Symbol("T.Peak")]
    Y[4] = training_df[index, :Max]
    Y[5] = training_df[index, :AUC]
    
    # setup parameters, and parameter bounds -
    if (isnothing(pₒ) == true)
        κ = [
            
            # default: hand fit set -
            0.061   0.01 10.0   ; # 1
            1.0     0.01 10.0   ; # 2
            1.0     0.01 10.0   ; # 3
            1.0     0.01 10.0   ; # 4
            1.0     0.01 10.0   ; # 5
            1.0     0.01 10.0   ; # 6
            1.0     0.01 10.0   ; # 7
            1.0     0.01 10.0   ; # 8
            0.70    0.01 10.0   ; # 9
            0.1     0.01 10.0   ; # 10
            0.045   0.01 10.0   ; # 11
            0.065   0.01 10.0   ; # 12
        ];

        # set default set as the start -
        pₒ = κ[:,1]
    end

    # what is my initial loss?
    lₒ = loss(pₒ, Y, model)
    
    for i ∈ 1:NTRIALS
        
        # generate a new candidate solution -
        Δᵢ = abs.(pₒ .* (1 .+ σ*rand(-1:1,12)))
        
        # compute the loss with this new solution -
        lᵢ = loss(Δᵢ, Y, model)

        # ok, is norm(lᵢ) < norm(lₒ)?
        if (norm(lᵢ) < norm(lₒ))
            
            # replace the "mom" -
            pₒ = Δᵢ

            # update the loss -
            lₒ = lᵢ

            @info "New soln for training set $(index) found at iteration $(i) with || lᵢ || = $(norm(lᵢ))"
        end
    end

    # run the sim w/the best parameters -
    # 1 - 9 : α vector
    model["α"] = pₒ[1:9]
    G = model["G"]
    idx = indexin(model, "FVIIa")   # 10
    G[idx, 4] = pₒ[10]
    idx = indexin(model, "AT")      # 11
    G[idx, 9] = pₒ[11]
    idx = indexin(model, "TFPI")    # 12
    G[idx, 1] = -1*pₒ[12]

    # run the model -
    (T,U) = evaluate(model)
    Xₘ = hcat(U...)
    Yₘ = model_output_vector(T,Xₘ[9,:]) # properties of the Thrombin curve 

    return (pₒ, Yₘ, Y)
end