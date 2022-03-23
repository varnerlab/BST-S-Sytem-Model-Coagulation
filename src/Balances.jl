function balances(dx, x, p, t)

    # grab data from the p dictionary -
    α = p["α"]
    G = p["G"]
    S = p["S"]
    number_of_dynamic_states = p["number_of_dynamic_states"]
    
    # get the static factors (and thier values)
    static_factors_array = p["static_factors_array"]
    
    # build the "state" array (dynamic | static)
    state_array = vcat(x,static_factors_array)

    # compute the kinetics - powerlaw
    rV = powerlaw(state_array,α,G)

    # compute the rhs -> store in a temp vector
    tmp = S*rV

    # populate the dx vector with the tmp vector -
    for i ∈ 1:number_of_dynamic_states
        dx[i] = tmp[i]
    end
end