using ITensors
using LinearAlgebra
using HDF5
using Statistics
using SparseArrays
include("utilities.jl")
ITensors.disable_warn_order()

inputs = Dict()
inputs["S"] = 1/2
inputs["x"] = 0.7
inputs["mg"] = 0.3
inputs["N"] = 4
inputs["ms"] = 10
inputs["ssl"] = [10]
inputs["cutoff"] = 1e-11
inputs["db"] = 0.01

# Get spin S number and create operators
S = inputs["S"]
ITensors.space(::SiteType"Spin") = Int(2*S+1)
function ITensors.op(::OpName"Sz_sp", ::SiteType"Spin", s::Index)

    d = dim(s)
    # S = div(d-1, 2) 

    return diagm(S:-1:-S)

end

function ITensors.op(::OpName"Sz2_sp", ::SiteType"Spin", s::Index)

    d = dim(s) 
    # S = div(d-1, 2)
    
    return diagm((S:-1:-S).^2)

end

function ITensors.op(::OpName"Sz", ::SiteType"Spin", s::Index)

    d = dim(s)
    # S = div(d-1, 2)
    tmp = diagm(S:-1:-S)

    return ITensor(tmp, s', dag(s))

end

function ITensors.op(::OpName"Sz2", ::SiteType"Spin", s::Index)

    d = dim(s) 
    # S = div(d-1, 2)
    tmp = diagm((S:-1:-S).^2)

    return ITensor(tmp, s', dag(s))

end

function ITensors.op(::OpName"S+", ::SiteType"Spin", s::Index)

    """
    See definition: https://easyspin.org/easyspin/documentation/spinoperators.html
    """

    d = dim(s) 
    # S = div(d-1, 2)
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
    # S = div(d-1, 2)
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

# Get the sites 
N = inputs["N"]
sites = get_sites(N)

# Get the initial state
rho = MPO(sites, "Id")
cutoff = inputs["cutoff"]
gauss_law_projector_mpo_no_charges = get_gauge_invariant_subspace_projector_MPO(S, sites, [0 for _ in 1:N])
rho = apply(gauss_law_projector_mpo_no_charges, apply(rho, hermitian_conjugate_mpo(gauss_law_projector_mpo_no_charges); cutoff = cutoff); cutoff = cutoff)
rho /= tr(rho)
orthogonalize!(rho, 1)

function evolve(rho, inputs, S)
    
    # List to write to file after including the states rho and corresponding beta
    current_beta = 0.0
    beta_list = Float64[current_beta]
    rho_list = MPO[rho]
    save_steps_list = inputs["ssl"]
    link_dims = zeros(Int64, length(save_steps_list) + 1, length(rho)-1)
    link_dims[1, :] = linkdims(rho)
    counter = 2 # where to put the next element in link_dims, incremented with every addition

    # Prepare the gates for the integration scheme: e^(-db*H_o/4)e^(-db*H_e/2)e^(-db*H_o/4) I e^(-db*H_o/4)e^(-db*H_e/2)e^(-db*H_o/4)
    delta_beta = inputs["db"]
    x = inputs["x"]
    mg = inputs["mg"]
    a = -delta_beta/4
    odd_gates = get_odd_gates(sites, a, x, mg, S) # odd/8
    even_gates = get_even_gates(sites, 2*a, x, mg, S) # even/4
    double_step_odd_gates = get_odd_gates(sites, 2*a, x, mg, S) # odd/4 (multiply last odd into next time step)

    # Define the Hamiltonian to track the energy
    H = get_Hamiltonian_MPO(x, mg, sites, S)

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
            push!(rho_list, apply(hermitian_conjugate_mpo(rho), rho; cutoff = cutoff))
            push!(beta_list, 2*current_beta)
            link_dims[counter, :] = linkdims(rho)
            counter += 1
        end
        println("Step $(step), t = $(time() - t), linkdims = $(linkdims(rho)), energy = $(inner(rho, H))")
        flush(stdout)
        t = time()
        
    end

    return rho_list, beta_list, link_dims

end
rho_list, beta_list, link_dims = evolve(rho, inputs, S)

# ------------------------------------------------------------------------------------------------------------------------

# Check that the evolved state corresponds to exp(-beta * H) normalized

# rho = rho_list[end]
# rho_m = mixed_mpo_to_matrix(rho)
# rho_m /= tr(rho_m)
# rho_m = sparse(rho_m)

# x, mg = inputs["x"], inputs["mg"]
# H = get_Hamiltonian_MPO(x, mg, sites)
# H_m = mixed_mpo_to_matrix(H)
# beta = beta_list[end]
# rho_m_sp = exp(-beta*H_m)
# gauss_law_projector_mpo_no_charges_m = mixed_mpo_to_matrix(gauss_law_projector_mpo_no_charges)
# rho_m_sp = gauss_law_projector_mpo_no_charges_m * rho_m_sp * gauss_law_projector_mpo_no_charges_m'
# rho_m_sp /= tr(rho_m_sp)
# rho_m_sp = sparse(rho_m_sp)

# # e1 = real(eigen(rho_m).values)
# # e2 = real(eigen(rho_m_sp).values)

# # display(eigen(rho_m).values)
# # display(eigen(rho_m_sp).values)

# # for i in 1:length(e1)
# #     if abs(e1[i] - e2[i]) > 1e-10
# #         println(i, ": ", abs(e1[i] - e2[i]))
# #     end
# # end

# println(norm(rho_m_sp - rho_m))

# ------------------------------------------------------------------------------------------------------------------------

# Check that the odd and even gates together sum up to the Hamiltonian

# x, mg = inputs["x"], inputs["mg"]
# H = get_Hamiltonian_MPO(x, mg, sites)
# H_m = mixed_mpo_to_matrix(H)

# odd, even = get_odd_gates_opsum(sites, x, mg, S), get_even_gates_opsum(sites, x, mg, S)
# total = odd + even
# H_m_1 = MPO(total, sites)
# H_m_1 = mixed_mpo_to_matrix(H_m_1)

# println(norm(H_m_1 - H_m))

# ------------------------------------------------------------------------------------------------------------------------

# Check that the Hamiltonian is Hermitian

# x, mg = inputs["x"], inputs["mg"]
# H = get_Hamiltonian_MPO(x, mg, sites)
# H_m = mixed_mpo_to_matrix(H)
# println(isapprox(H_m', H_m))

# ------------------------------------------------------------------------------------------------------------------------

# Check that the ground state of the Hamiltonian has an energy equivalent to the energy we get from getting a rho at
# very larger beta

# function get_DMRG_results(cutoff, N, x, mg, nsweeps, energy_tol, dist, lambd)

#     external_charges = get_external_charges_list(dist)
#     H = MPO(get_Hamiltonian_opsum(x, mg, sites) + lambd*get_gauss_law_squared_sum_opsum(N, external_charges), sites)

#     psi0 = randomMPS(sites, linkdims = 2)
#     obs = DMRGObserver(;energy_tol = energy_tol)
    
#     gs_energy, gs = dmrg(H, psi0; nsweeps, cutoff, observer=obs, outputlevel = 1)

#     penalty_term_mpo = MPO(get_gauss_law_squared_sum_opsum(N, external_charges), sites)
#     penalty_term = inner(gs', penalty_term_mpo, gs)
#     total_z = sum(expect(gs, "Z"; sites = 1:2:2*N-1))

#     return gs_energy, gs, penalty_term, total_z

# end

# x, mg = inputs["x"], inputs["mg"]
# nsweeps, energy_tol, dist, lambd = 1000, 1e-16, 0, 100 # make sure the nsweeps are large enough
# gs_energy, gs, penalty_term, total_z = get_DMRG_results(cutoff, N, x, mg, nsweeps, energy_tol, dist, lambd)

# H = get_Hamiltonian_MPO(x, mg, sites)
# rho = rho_list[end]
# rho /= tr(rho)
# println(inner(rho, H) - gs_energy)

# ------------------------------------------------------------------------------------------------------------------------

# Check the Green's function 

# function get_DMRG_results(cutoff, N, x, mg, nsweeps, energy_tol, dist, lambd)

#     external_charges = get_external_charges_list(dist)
#     H = MPO(get_Hamiltonian_opsum(x, mg, sites) + lambd*get_gauss_law_squared_sum_opsum(N, external_charges), sites)

#     psi0 = randomMPS(sites, linkdims = 2)
#     obs = DMRGObserver(;energy_tol = energy_tol)
    
#     gs_energy, gs = dmrg(H, psi0; nsweeps, cutoff, observer=obs, outputlevel = 1)

#     penalty_term_mpo = MPO(get_gauss_law_squared_sum_opsum(N, external_charges), sites)
#     penalty_term = inner(gs', penalty_term_mpo, gs)
#     total_z = sum(expect(gs, "Z"; sites = 1:2:2*N-1))

#     return gs_energy, gs, penalty_term, total_z

# end

# x, mg = inputs["x"], inputs["mg"]
# nsweeps, energy_tol, dist, lambd = 1000, 1e-18, 0, 100 # make sure the nsweeps are large enough
# gs_energy, gs, penalty_term, total_z = get_DMRG_results(cutoff, N, x, mg, nsweeps, energy_tol, dist, lambd)

# rho = rho_list[end]

# # g = get_greens_function_opsum(3, 5)
# # g = MPO(g, sites)
# # println(inner(rho, g))
# # println(inner(gs', g, gs))

# res1, idxs1 = get_greens_function_expectation_values(gs, sites)
# res2, idxs2 = get_greens_function_expectation_values(rho, sites)

# # for i in 1:length(res1)
# #     println(res1[i], " ", res2[i])
# # end

# println(idxs1)

# ------------------------------------------------------------------------------------------------------------------------

# N = 8
# sites = get_sites(N)
# opsum = get_greens_function_opsum(5, 11)
# println(opsum)
# g = MPO(opsum, sites)

# ------------------------------------------------------------------------------------------------------------------------

println("Finished.")
