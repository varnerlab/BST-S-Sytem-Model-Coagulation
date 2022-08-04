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

# Load the synthetic data 10K as a DataFrame, along with the actual data -
df_synthetic_data = CSV.read(joinpath(_PATH_TO_DATA, "Training-Synthetic-Thrombin-TF-10K.csv"), DataFrame)
df_actual_data = CSV.read(joinpath(_PATH_TO_DATA, "Training-Thrombin-TF-Transformed-w-Labels.csv"), DataFrame)[!,3:end]

# which patient records are similar to the synthetic data?
P = length(df_actual_data[:,:II]);
𝒫 = length(df_synthetic_data[:,:II]);
DS = zeros(𝒫,P)
rng_compare = 1:16
for i ∈ 1:𝒫

    # get the ith record -
    xᵢ = Vector(df_synthetic_data[i, rng_compare])
    
    for j ∈ 1:P

        # get the ith record -
        xⱼ = Vector(df_actual_data[j, rng_compare])
        DS[i,j] = d(xᵢ,xⱼ; bandwidth=0.001)
    end
end

# setup plot for synthetic data -
idx_synth = sortperm(DS[:,1]; rev=true)
N = 30
idx_set_shown = Set{Int64}()
for p ∈ 1:N
	
	sim_filename = "SIM-Synthetic-1k-P$(idx_synth[p]).csv"
	pset_filename = "PSET-Synthetic-1k-P$(idx_synth[p]).csv"

	# load the pset to check the error -
	pset_data = CSV.read(joinpath(_PATH_TO_SYNTHETIC_ENSEMBLE, pset_filename), DataFrame)
	synth_data = CSV.read(joinpath(_PATH_TO_SYNTHETIC_ENSEMBLE, sim_filename), DataFrame)
	
	if (pset_data[1,:parameters]<1.5)
		plot!(synth_data[:,:T], synth_data[:,:FIIa], label="",
			c=colorant"#89CCE2", bg="aliceblue", background_color_outside="white", framestyle = :box)
		push!(idx_set_shown, idx_synth[p])
	end
end

# setup plot for actual data -
actual_data = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "SIM-Actual-P1.csv"), DataFrame)
plot!(actual_data[:,:T], actual_data[:,:FIIa], c=:red, lw=2, label="FIIa actual", fg_legend = :transparent, ylim=(0.0, 320.0))
xlabel!("Time (min)", fontsize=18)
ylabel!("FIIa concentration (nM)", fontsize=18)
gui()

# dump fig file to disk -
# filename = "FIIa-Sort-All-N30.pdf"
# savefig(joinpath(_PATH_TO_BASE_FIGS, filename))