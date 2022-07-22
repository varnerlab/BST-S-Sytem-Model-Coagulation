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
(number_of_records, _) = size(full_df) 

# what are the col names? (exclude subject id and visit id)
data_col_symbols = Symbol.(names(full_df)[1:end])
number_of_fields = length(data_col_symbols)

# setup scale factor dictionary to convert to concentration units -
SF = 1e9
scale_factor_dict = Dict()
scale_factor_dict[:II] = (1.4e-6)*SF*(1/100)
scale_factor_dict[:V] = (2e-8)*SF*(1/100)
scale_factor_dict[:VII] = (1e-8)*SF*(1/100)
scale_factor_dict[:VIII] = (7e-10)*SF*(1/100)
scale_factor_dict[:IX] = (9e-8)*SF*(1/100)
scale_factor_dict[:X] = (1.6e-7)*SF*(1/100)
scale_factor_dict[:XI] = (1e-8)*SF*(1/100)
scale_factor_dict[:XII] = (1e-8)*SF*(1/100)
scale_factor_dict[:AT] = (3.4e-6)*SF*(1/100)
scale_factor_dict[:PC] = (6.3e-8)*SF*(1/100)
scale_factor_dict[:TFPI] = 1.0
scale_factor_dict[:Lagtime] = 1.0
scale_factor_dict[:Peak] = 1.0
scale_factor_dict[Symbol("T.Peak")] = 1.0
scale_factor_dict[:Max] = 1.0
scale_factor_dict[:AUC] = 1.0
scale_factor_dict[:subjid] = 1
scale_factor_dict[:visitid] = 1


# initialize -
transformed_df = DataFrame()
for rᵢ ∈ 1:number_of_records


    tmp = Vector()
    for fᵢ ∈ 1:number_of_fields
        field_symbol = data_col_symbols[fᵢ]
        value = scale_factor_dict[field_symbol]*full_df[rᵢ,field_symbol]
        push!(tmp, value)
    end

    # create new record -
    transformed_tuple = (
        subjid = tmp[1], 
        visitid = tmp[2],
        II = tmp[3],
        V = tmp[4],
        VII = tmp[5],
        VIII = tmp[6],
        IX = tmp[7],
        X = tmp[8],
        XI = tmp[9],
        XII = tmp[10],
        AT = tmp[11],
        PC = tmp[12],
        TFPI = tmp[13],
        Lagtime = tmp[14],
        Peak = tmp[15],
        TPeak = tmp[16],
        Max = tmp[17],
        AUC = tmp[18]
    )
    push!(transformed_df, transformed_tuple)
end

# dump sample data to disk -
path_to_transformed_data = joinpath(_PATH_TO_DATA, "Training-Thrombin-TF-Transformed-w-Labels.csv")
CSV.write(path_to_transformed_data, transformed_df)