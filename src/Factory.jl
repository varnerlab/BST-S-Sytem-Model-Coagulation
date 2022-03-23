import Base.+

function +(buffer::Array{String,1}, line::String)
    push!(buffer, line)
end

function parse_species_record(buffer::Array{String,1})::Array{String,1}

    # initialize -
    species_symbol_array = Array{String,1}()

    # main -
    for line ∈ buffer

        # split around the ,
        tmp_array = String.(split(line,","))

        for species ∈ tmp_array
            
            species_symbol_string =  species |> lstrip |> rstrip
            if (in(species_symbol_string, species_symbol_array) == false)
               push!(species_symbol_array,species_symbol_string)
            end
        end
    end

    return species_symbol_array
end

function parse_structure_section(buffer::Array{String,1})::Array{Dict{String,Any}}

    # initialize -
    record_array = Array{Dict{String,Any},1}()

    for line ∈ buffer
        
        # split around the ,
        record_components = String.(split(line,","))
        name = record_components[1]
        left_phrase = record_components[2]
        right_phrase = record_components[3]

        # add -
        tmp_dict = Dict{String,Any}()
        tmp_dict["id"] = name
        tmp_dict[left_phrase] = -1.0
        tmp_dict[right_phrase] = 1.0

        # store -
        push!(record_array,tmp_dict)
    end

    # return -
    return record_array
end

function parse_rate_section(buffer::Array{String,1})::Array{Dict{String,Any}}

    # initialize -
    record_array = Array{Dict{String,Any},1}()

    for record ∈ buffer
        
        # slpit around ::
        record_components = String.(split(record,"::"))

        # process each component -
        name = record_components[1]
        factors = record_components[2]

        # factors are a set of stuff, so we need to split around the ,
        factor_list = factors[2:end-1] # this get's rid of the {}
        factor_list_components = String.(split(factor_list,","))

        # ok, let's rock ...
        tmp_dict = Dict{String,Any}()
        tmp_dict["name"] = name
        for factor ∈ factor_list_components
            tmp_dict[factor] = 1.0
        end

        # grab -
        push!(record_array, tmp_dict)
    end

    # return -
    return record_array
end

function build_stoichiometric_matrix(list_of_dynamic_species::Array{String,1},
    reactions::Array{Dict{String,Any}})::Array{Float64,2}

    # initialize -
    ℳ = length(list_of_dynamic_species)
    ℛ = length(reactions)
    S = zeros(ℳ,ℛ)

    for species_index ∈ 1:ℳ
        
        # get the species -
        species_symbol = list_of_dynamic_species[species_index]

        for reaction_index = 1:ℛ
            
            reaction_dictionary = reactions[reaction_index]
            if (haskey(reaction_dictionary, species_symbol) == true)
                S[species_index,reaction_index] = reaction_dictionary[species_symbol]
            end
        end
    end

    # return -
    return S
end

function build_exponent_matrix(list_of_factors::Array{String,1}, 
    reaction_factors::Array{Dict{String,Any}})

    # initialize -
    ℳ = length(list_of_factors)
    ℛ = length(reaction_factors)
    factor_matrix = zeros(ℳ,ℛ)

    # main -
    for (i,f) ∈ enumerate(list_of_factors)        
        for (j,d) ∈ enumerate(reaction_factors)

            # ok: so in this d, to we have f?
            if (haskey(d,f) == true)
                factor_matrix[i,j] = d[f]
            end
        end
    end

    # return -
    return factor_matrix
end

function extract_model_section(file_buffer_array::Array{String,1},
    start_section_marker::String,end_section_marker::String)::Array{String,1}

    # initialize -
    section_buffer = String[]

    # find the SECTION START AND END -
    section_line_start = 1
    section_line_end = 1
    for (index, line) in enumerate(file_buffer_array)

        if (occursin(start_section_marker, line) == true)
            section_line_start = index
        elseif (occursin(end_section_marker, line) == true || length(line) == index)
            section_line_end = index
        end
    end

    for line_index = (section_line_start+1):(section_line_end-1)
        line_item = file_buffer_array[line_index]
        push!(section_buffer, line_item)
    end

    # return -
    return section_buffer
end

function build_default_model_dictionary(model_buffer::Array{String,1})::Dict{String,Any}

    # initialize -
    model_dict = Dict{String,Any}()

    # get the sections of the model file -
    dynamic_section = extract_model_section(model_buffer,"#pragma::dynamic","#dynamic::end")
    static_section = extract_model_section(model_buffer,"#pragma::static","#static::end")
    structure_section = extract_model_section(model_buffer,"#pragma::structure","#structure::end")
    rate_section = extract_model_section(model_buffer,"#pragma::rate","#rate::end")

    # get list of dynamic species -
    list_of_dynamic_species = parse_species_record(dynamic_section)
    number_of_dynamic_states = length(list_of_dynamic_species)

    # get list of static species -
    list_of_static_species = parse_species_record(static_section)
    number_of_static_states = length(list_of_static_species)
    static_factors_array = zeros(number_of_static_states)

    # total species list -
    total_species_list = vcat(list_of_dynamic_species,list_of_static_species)

    # build the stoichiometric (connectivity) array -
    structure_dict_array = parse_structure_section(structure_section)
    S = build_stoichiometric_matrix(list_of_dynamic_species, structure_dict_array)

    # build rate exponent array -
    rate_dict_array = parse_rate_section(rate_section)
    G = build_exponent_matrix(total_species_list,rate_dict_array)
    
    # build the rate constant array -
    α = ones(length(rate_dict_array))

    # populate -
    model_dict["number_of_dynamic_states"] = number_of_dynamic_states
    model_dict["number_of_static_states"] = number_of_static_states
    model_dict["list_of_dynamic_species"] = list_of_dynamic_species
    model_dict["list_of_static_fators"] = list_of_static_species
    model_dict["total_species_list"] = total_species_list
    model_dict["static_factors_array"] = static_factors_array
    model_dict["S"] = S
    model_dict["G"] = G
    model_dict["α"] = α

    # return -
    return model_dict
end

function read_model_file(path_to_file::String)::Array{String,1}

    # initialize -
    model_file_buffer = String[]
    model_buffer = Array{String,1}()

    # Read in the file -
    open("$(path_to_file)", "r") do file
        for line in eachline(file)
            +(model_file_buffer,line)
        end
    end

    # process -
    for line ∈ model_file_buffer
        
        # skip comments and empty lines -
        if (occursin("//", line) == false && 
            isempty(line) == false)
        
            # grab -
            push!(model_buffer,line)
        end
    end

    # return -
    return model_buffer
end