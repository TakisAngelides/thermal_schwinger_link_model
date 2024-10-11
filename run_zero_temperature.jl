using ITensors
using LinearAlgebra
using HDF5
using SparseArrays
using Plots
using LaTeXStrings
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

# Build the operators based on spin S
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

# Function to perform constrained DMRG with the constrain being Gauss law in zero external charges sector
println("Now getting ground state with constrained DMRG ", now())
flush(stdout)
function get_DMRG_results(inputs, results_file)

    N = inputs["N"]
    x = inputs["x"]
    mg = inputs["mg"]
    cutoff = inputs["cutoff"]
    nsweeps = inputs["ns"]
    energy_tol = inputs["et"]
    lambd = inputs["lambd"]
    greens_fn_flag = parse(Bool, inputs["gff"])

    # Get the Hamiltonian
    sites = get_sites(N)
    gauss_law_squared_sum_opsum = get_gauss_law_squared_sum_opsum(N, [0 for _ in 1:N])
    println("Now preparing the Hamiltonian MPO for the constrained DMRG ", now())
    H = MPO(get_Hamiltonian_opsum(x, mg, sites) + lambd*gauss_law_squared_sum_opsum, sites)
    println("Finished preparing the Hamiltonian MPO for the constrained DMRG and it has linkdims = ", linkdims(H), " ", now())

    # Get the initial state for the DMRG and define the observer 
    psi0 = randomMPS(sites, linkdims = 2)
    obs = DMRGObserver(;energy_tol = energy_tol)

    # Get the ground state and energy with DMRG
    gs_energy, gs = dmrg(H, psi0; nsweeps, cutoff, observer=obs, outputlevel = 1)

    # Measure penalty term and total Z and greens functions
    penalty_term_mpo = MPO(gauss_law_squared_sum_opsum, sites)
    penalty_term = inner(gs', penalty_term_mpo, gs)
    total_z = sum(expect(gs, "Z"; sites = 1:2:2*N-1))
    if greens_fn_flag
        greens_functions, greens_functions_distances = get_greens_function_expectation_values(gs, sites)
    end

    # Write results to file
    write(results_file, "energy", gs_energy)
    write(results_file, "total_z", total_z)
    write(results_file, "penalty_term", penalty_term)
    write(results_file, "gs", gs)
    if greens_fn_flag
        write(results_file, "greens_functions", greens_functions)
        write(results_file, "greens_functions_distances", greens_functions_distances)
    end
    close(results_file)

end

get_DMRG_results(inputs, results_file)
println("Finished getting ground state with constrained DMRG ", now())
flush(stdout)
