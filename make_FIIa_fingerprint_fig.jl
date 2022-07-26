# include -
include("Include.jl")

# load *actual* simulation data -
actual_data = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "SIM-Actual-P1.csv"), DataFrame)

U = actual_data[:,:FIIa]
Z = zeros(length(actual_data[:,:T]))

plot(actual_data[:,:T], actual_data[:,:FIIa], c=colorant"#EF4035", lw=2, label="", ylim=(0.0,260.0))
plot!(actual_data[:,:T], Z, fillrange=U, c=colorant"#89CCE2", fillalpha = 0.4, lw=2, label="")

xlabel!("Time (min)", fontsize=18)
ylabel!("FIIa concentration (nM)", fontsize=18)

# dump fig file to disk -
filename = "FIIa-Fingerprint-Fig.pdf"
savefig(joinpath(_PATH_TO_BASE_FIGS, filename))