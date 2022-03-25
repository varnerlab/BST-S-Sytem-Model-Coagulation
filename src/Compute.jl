function τpeak(T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # find the index of the max value -
    idx = argmax(X)
    return T[idx]
end

function τlag(rule::Function, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(rule, X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end

function auc(T::Array{Float64,1}, X::Array{Float64,1})::Float64    
    return integrate(T, X)
end

function max_FIIa_rate(T::Array{Float64,1}, X::Array{Float64,1})::Float64
    
    # intialize -
    N = 600 # sim to 6 min
    diff_array = Array{Float64,1}(undef, N - 1)
    ΔT = T[2] - T[1]

    # compute the ΔX/ΔT -
    for i ∈ 1:(N-1)
        
        # compute diff -
        ΔX = (X[i+1] - X[i])*(1/ΔT)
        diff_array[i] = ΔX
    end
    
    # find the max -
    return maximum(diff_array)
end

function model_output_vector(T::Array{Float64,1}, X::Array{Float64,1})::Array{Float64,1}

    # initialize -
    output_vector = Array{Float64,1}()

    # Lagtime,Peak,T.Peak,Max,AUC
    lagtime = τlag(x->x>10.0, T, X)
    FIIa_peak = maximum(X)
    time_to_peak = τpeak(T, X)
    FIIa_rate = max_FIIa_rate(T, X)
    area_under_curve = auc(T, X)

    # package -
    push!(output_vector, lagtime)
    push!(output_vector, FIIa_peak)
    push!(output_vector, time_to_peak)
    push!(output_vector, FIIa_rate)
    push!(output_vector, area_under_curve)

    # return -
    return output_vector
end