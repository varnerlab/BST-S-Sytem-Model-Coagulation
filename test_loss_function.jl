# load the include -
include("Include.jl")

# load the training data -
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Synthetic-Thrombin-TF-1K.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

# size of training set -
(R,C) = size(training_df)

# build the model structure -
path_to_model_file = joinpath(_PATH_TO_MODEL, "Coagulation.net")
model_buffer = read_model_file(path_to_model_file)

# build the default model structure -
model = build_default_model_dictionary(model_buffer)

# main simulation loop -
SF = 1e9
for i ∈ 1:10

    # build new model -
    dd = deepcopy(model)

    # setup static -
    sfa = dd["static_factors_array"]
    sfa[1] = training_df[i,:TFPI]       # 1 TFPI
    sfa[2] = training_df[i,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF               # 3 TF
    sfa[6] = 0.005                      # 6 TRAUMA

    # grab the multiplier from the data -
    ℳ = dd["number_of_dynamic_states"]
    xₒ = zeros(ℳ)
    xₒ[1] = training_df[i, :II]         # 1 FII 
    xₒ[2] = training_df[i, :VII]        # 2 FVII 
    xₒ[3] = training_df[i, :V]          # 3 FV
    xₒ[4] = training_df[i, :X]          # 4 FX
    xₒ[5] = training_df[i, :VIII]       # 5 FVIII
    xₒ[6] = training_df[i, :IX]         # 6 FIX
    xₒ[7] = training_df[i, :XI]         # 7 FXI
    xₒ[8] = training_df[i, :XII]        # 8 FXII 
    xₒ[9] = (1e-14)*SF                  # 9 FIIa
    dd["initial_condition_vector"] = xₒ

    # generate a κ vector -
    # κ = [
    #     0.061   ; # 1
    #     1.0     ; # 2
    #     1.0     ; # 3
    #     1.0     ; # 4
    #     1.0     ; # 5
    #     1.0     ; # 6
    #     1.0     ; # 7
    #     1.0     ; # 8
    #     0.70    ; # 9
    #     0.1     ; # 10
    #     0.045   ; # 11
    #     -0.065  ; # 12
    # ];

    κ = [
        0.7478700663653719
        1.7676181108566378
        1.6647348616139388
        1.777505405930172
        1.7730138521648193
        1.7272541526011016
        1.691479136564128
        1.7180148264997843
        1.4655502130381282
        0.3795056851485291
        0.6631825691926526
        0.7002036750754579
    ];

    # what is the output array?
    Y = Array{Float64,1}(undef,5)
    Y[1] = training_df[i,:Lagtime]
    Y[2] = training_df[i,:Peak]
    Y[3] = training_df[i,Symbol("T.Peak")]
    Y[4] = training_df[i,:Max]
    Y[5] = training_df[i,:AUC]
    
    # compute the loss -
    l = loss(κ, Y, dd)

    @show l
end

