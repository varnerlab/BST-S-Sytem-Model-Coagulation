# Where are we?
_BASE_PATH = pwd()
_PATH_TO_DATA = joinpath(_BASE_PATH, "data")

# Package that we need -
using DataFrames
using CSV
using Distributions
using Statistics
using LinearAlgebra

# load the actual data -
path_to_experimental_data = joinpath(_PATH_TO_DATA, "Training-Thrombin-TF.csv")
full_df = CSV.read(path_to_experimental_data, DataFrame)

# what are the col names? (exclude subject id and visit id)
data_col_symbols = Symbol.(names(full_df)[3:end])
data_col_names_tuple = Tuple(Symbol.(names(full_df)[3:end]))
number_of_fields = length(data_col_symbols)

# what visit are we looking at?
# visit_df = filter(:visitid => x->(x==2 || x==3), full_df)

# sample -
number_of_synthetic_samples = 1000
data_array = Matrix(full_df[!, 3:end])
μ = mean(data_array,dims=1)
Σ = cov(data_array)
D = MvNormal(reshape(μ,16), Σ)
synthetic_sampled_data = rand(D, number_of_synthetic_samples)

# create a new data frame -
synthetic_data_frame = DataFrame()
for i ∈ 1:number_of_synthetic_samples

    # grab a col of data (a random sample)
    data_col = max.(0.0, synthetic_sampled_data[:,i])
    
    # build a data tuple to add to df -
    data_tuple = NamedTuple{data_col_names_tuple}(data_col)
    push!(synthetic_data_frame, data_tuple)
end

# setup scale factor dictionary -
SF = 1e9
scale_factor_dict = Dict()
scale_factor_dict[:II] = (1.4e-6)*SF
scale_factor_dict[:V] = (2e-8)*SF
scale_factor_dict[:VII] = (1e-8)*SF
scale_factor_dict[:VIII] = (7e-10)*SF
scale_factor_dict[:IX] = (9e-8)*SF
scale_factor_dict[:X] = (1.6e-7)*SF
scale_factor_dict[:XI] = (1e-8)*SF
scale_factor_dict[:XII] = (1e-8)*SF
scale_factor_dict[:AT] = (3.4e-6)*SF
scale_factor_dict[:PC] = (6.3e-8)*SF
scale_factor_dict[:TFPI] = 1.0


# ok, finally - let's convert these data into the correct units - TFPI is already in nM, the rest 
# are in percentages of nominal -
for i ∈ 1:number_of_synthetic_samples

    for j ∈ 1:number_of_fields
        
        # look up the scale factor -
        fld_symbol = data_col_symbols[j]
        if (fld_symbol != :TFPI && haskey(scale_factor_dict,fld_symbol) == true) # TFPI is already in the correct units -
            SF_val = scale_factor_dict[fld_symbol]
            old_value = synthetic_data_frame[i,fld_symbol]
            new_value = SF_val*(old_value/100)
            synthetic_data_frame[i,fld_symbol] = new_value
        end 
    end
end



# dump sample data to disk -
path_to_synthetic_data = joinpath(_PATH_TO_DATA, "Training-Synthetic-Thrombin-TF.csv")
CSV.write(path_to_synthetic_data, synthetic_data_frame)