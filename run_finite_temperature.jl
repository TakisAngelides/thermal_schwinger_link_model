using ITensors
using LinearAlgebra
using HDF5
using Statistics
using SparseArrays
using Dates
include("utilities.jl")
ITensors.disable_warn_order()

# Input arguments for file and opening the results h5 file
println("Now getting the input arguments for the file and opening the results h5 file ", now())
flush(stdout)
path_to_project_number = ARGS[1] 
path_to_inputs_h5 = "$(path_to_project_number)/inputs.h5" # Where the h5 file containing the inputs is
job_id = ARGS[2] # The job id to get the inputs within this h5 file for this specific job
path_to_results = "$(path_to_project_number)/HDF5/$(job_id).h5" # The path to the h5 file to save all results of states
results_file = h5open(path_to_results, "w")
println("Finished getting the input arguments for the file and opening the results h5 file ", now())
flush(stdout)

# Build a dictionary of the inputs
println("Now building a dictionary of the inputs ", now())
flush(stdout)
inputs_file = h5open(path_to_inputs_h5, "r")
group = inputs_file["$(job_id)"]
group_attributes = attributes(group)
inputs = Dict()
for key in keys(group_attributes)
    inputs[key] = read(group_attributes[key])
end
close(inputs_file)
println("The inputs are ", inputs)
flush(stdout)
println("Finished building a dictionary of the inputs ", now())
flush(stdout)

# Get spin S number and create operators
println("Now getting the spin S dependent functions ", now())
flush(stdout)
S = inputs["S"]
ITensors.space(::SiteType"Spin") = Int(2*S+1)
function ITensors.op(::OpName"Sz_sp", ::SiteType"Spin", s::Index)

    d = dim(s)
    S = div(d-1, 2) 

    return diagm(S:-1:-S)

end

function ITensors.op(::OpName"Sz2_sp", ::SiteType"Spin", s::Index)

    d = dim(s) 
    S = div(d-1, 2)
    
    return diagm((S:-1:-S).^2)

end

function ITensors.op(::OpName"Sz", ::SiteType"Spin", s::Index)

    d = dim(s)
    S = div(d-1, 2)
    tmp = diagm(S:-1:-S)

    return ITensor(tmp, s', dag(s))

end

function ITensors.op(::OpName"Sz2", ::SiteType"Spin", s::Index)

    d = dim(s) 
    S = div(d-1, 2)
    tmp = diagm((S:-1:-S).^2)

    return ITensor(tmp, s', dag(s))

end

function ITensors.op(::OpName"S+", ::SiteType"Spin", s::Index)

    """
    See definition: https://easyspin.org/easyspin/documentation/spinoperators.html
    """

    d = dim(s) 
    S = div(d-1, 2)
    m = zeros(Float64, d, d)
    Sz_vec = S:-1:-S
    for row in 1:d
        for col in 1:d

            if col == row + 1

                m_row, m_col = Sz_vec[row], Sz_vec[col]

                m[row, col] = sqrt(S*(S+1) - m_row*m_col)

            end

        end
    end

    return ITensor(m, s', dag(s))

end

function ITensors.op(::OpName"S-", ::SiteType"Spin", s::Index)

    d = dim(s)
    S = div(d-1, 2)
    m = zeros(Float64, d, d)
    Sz_vec = S:-1:-S
    for row in 1:d
        for col in 1:d

            if row == col + 1

                m_row, m_col = Sz_vec[row], Sz_vec[col]

                m[row, col] = sqrt(S*(S+1) - m_row*m_col)

            end

        end
    end

    return ITensor(m, s', dag(s))

end
println("Finished getting the spin S dependent functions ", now())
flush(stdout)

# Get the initial state
println("Now getting the initial state ", now())
flush(stdout)
N = inputs["N"]
sites = get_sites(N)
rho = MPO(sites, "Id")
cutoff = inputs["cutoff"]
gauss_law_projector_mpo_no_charges = get_gauge_invariant_subspace_projector_MPO(S, sites, [0 for _ in 1:N])
rho = apply(gauss_law_projector_mpo_no_charges, apply(rho, hermitian_conjugate_mpo(gauss_law_projector_mpo_no_charges); cutoff = cutoff); cutoff = cutoff)
rho /= tr(rho)
orthogonalize!(rho, 1)
println("Finished getting the initial state ", now())
flush(stdout)

# Evolution to get the thermal state
println("Now starting the evolution ", now())
flush(stdout)
function evolve(rho, inputs, S, results_file)
    
    # Define the Hamiltonian to track the energy
    x = inputs["x"]
    mg = inputs["mg"]
    H = get_Hamiltonian_MPO(x, mg, sites)

    # List to write to file after including the states rho and corresponding beta
    current_beta = 0.0
    beta_list = Float64[current_beta]
    energy_list = ComplexF64[inner(rho, H)]
    save_steps_list = inputs["ssl"]
    counter = 1
    link_dims_list = zeros(Int64, length(save_steps_list) + 1, length(rho)-1)
    link_dims_list[1, :] = linkdims(rho)
    greens_fn_flag = parse(Bool, inputs["gff"])
    if greens_fn_flag
        greens_functions_list = zeros(ComplexF64, length(save_steps_list) + 1, div(N, 2))
        greens_functions_distances_list = zeros(Int64, length(save_steps_list) + 1, div(N, 2))
        greens_functions, idxs = get_greens_function_expectation_values(rho, sites)
        greens_functions_list[1, :] = greens_functions
        greens_functions_distances_list[1, :] = idxs
    end
    counter += 1

    # Prepare the gates for the integration scheme: e^(-db*H_o/4)e^(-db*H_e/2)e^(-db*H_o/4) I e^(-db*H_o/4)e^(-db*H_e/2)e^(-db*H_o/4)
    delta_beta = inputs["db"]
    a = -delta_beta/8
    odd_gates = get_odd_gates(sites, a, x, mg, S) # odd/8
    even_gates = get_even_gates(sites, 2*a, x, mg, S) # even/4
    double_step_odd_gates = get_odd_gates(sites, 2*a, x, mg, S) # odd/4 (multiply last odd into next time step)

    # Start the evolution
    t = time()
    max_steps = inputs["ms"]
    for step in 1:max_steps
    
        if step == 1
            apply_odd_gates!(odd_gates, rho; cutoff = cutoff)
            apply_even_gates!(even_gates, rho; cutoff = cutoff)
        elseif step == max_steps
            apply_odd_gates!(double_step_odd_gates, rho; cutoff = cutoff)
            apply_even_gates!(even_gates, rho; cutoff = cutoff)
            apply_odd_gates!(odd_gates, rho; cutoff = cutoff)
        else
            apply_odd_gates!(double_step_odd_gates, rho; cutoff = cutoff)
            apply_even_gates!(even_gates, rho; cutoff = cutoff)    
        end

        # Fix positivity and hermiticity with rho -> rho dagger * rho
        rho /= tr(rho) # this is also required for stability after applying the above fix of positivity and hermiticity

        # Increment beta
        current_beta += delta_beta

        # Write data to lists
        if step in save_steps_list
            rho_to_save = apply(hermitian_conjugate_mpo(rho), rho; cutoff = cutoff) # fixed positivity and doubles the beta reached
            rho_to_save /= tr(rho_to_save) # fixes the trace to 1
            link_dims_list[counter, :] = linkdims(rho_to_save)
            write(results_file, "rho_$(step)", rho_to_save)
            push!(beta_list, 2*current_beta) # we save twice the beta because of the first line in this if statement rho(beta) -> rho^dagger(beta/2)*rho(beta/2)
            push!(energy_list, inner(rho_to_save, H))
            if greens_fn_flag
                greens_functions, greens_functions_distances = get_greens_function_expectation_values(rho_to_save, sites)
                greens_functions_list[counter, :] = greens_functions
                greens_functions_distances_list[counter, :] = greens_functions_distances
            end
            counter += 1
        end
        println("Step $(step), t = $(time() - t), linkdims = $(linkdims(rho))")
        flush(stdout)
        t = time()
        
    end

    # Write tracked observables to file
    write(results_file, "beta_list", beta_list)
    write(results_file, "link_dims_list", link_dims_list)
    if greens_fn_flag
        write(results_file, "greens_functions_list", greens_functions_list)
        write(results_file, "greens_functions_distances_list", greens_functions_distances_list)
    end
    
end
evolve(rho, inputs, S, results_file)
println("Finished the evolution ", now())
flush(stdout)
