# include -
include("Include.jl")

# distance -
function d(synthetic::Array{Float64,1}, actual::Array{Float64,1}; bandwidth::Float64 = 0.01)::Float64

	# setup -
	l = bandwidth;
	xᵢ = synthetic
	x̄ⱼ = actual

	# compute the delta -
	δ = xᵢ .- x̄ⱼ

	# return -
	return exp(-(1/2*l^2)*(norm(δ)^2))
end;

# load the data -
df_actual_data_full = CSV.read(joinpath(_PATH_TO_DATA, "Training-Thrombin-TF-Transformed-w-Labels.csv"), DataFrame)
df_actual_data = CSV.read(joinpath(_PATH_TO_DATA, "Training-Thrombin-TF-Transformed-w-Labels.csv"), DataFrame)[!,3:end]
P = length(df_actual_data[:,:II]);
DA = zeros(P,P)
rng_compare = 1:16
for i ∈ 1:P

    # get the ith record -
    xᵢ = Vector(df_actual_data[i, rng_compare])
    
    for j ∈ 1:P

        # get the ith record -
        xⱼ = Vector(df_actual_data[j, rng_compare])
        DA[i,j] = d(xᵢ,xⱼ; bandwidth=0.001)
    end
end
