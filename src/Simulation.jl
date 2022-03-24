function evaluate(model::Dict{String,Any}; 
    tspan::Tuple{Float64,Float64} = (0.0,20.0), Δt::Float64 = 0.01)

    # get stuff from model -
    xₒ = model["initial_condition_vector"]

    # setup the solver -
    prob = ODEProblem(balances, xₒ, tspan, model; saveat = Δt)
    soln = solve(prob)

    # get the results from the solver -
    T = soln.t
    X = soln.u

    # return -
    return (T,X)
end

function loss(κ::Array{Float64,1}, Y::Array{Float64,1},  model::Dict{String,Any})

    # get the G matrix -
    G = model["G"]

    # κ map
    # 1 - 9 : α vector
    model["α"] = κ[1:9]
    idx = indexin(model, "FVIIa")   # 10
    G[idx, 4] = κ[10]
    idx = indexin(model, "AT")      # 11
    G[idx, 9] = κ[11]
    idx = indexin(model, "TFPI")    # 12
    G[idx, 1] = κ[12]

    # run the model -
    (T,U) = evaluate(model)
    Xₘ = hcat(U...)

    # compute the model ouput -
    Yₘ = model_output_vector(T,Xₘ[9,:]) # properties of the Thrombin curve 
    ϵ = ((Y .- Yₘ)./Y).^2

    # return -
    return ϵ
end