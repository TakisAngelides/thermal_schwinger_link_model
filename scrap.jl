using ITensors
using LinearAlgebra
using HDF5
using Statistics
using SparseArrays
using Plots
using StatsPlots
include("utilities.jl")
ITensors.disable_warn_order()

# --------------------------------------------------------------------------------------------------------------------------------

# f = h5open("DAGS/2/HDF5/31.h5", "r")

# fi = h5open("DAGS/2/inputs.h5", "r")
# group = fi["31"]
# group_attributes = attributes(group)
# inputs = Dict()
# for key in keys(group_attributes)
#     inputs[key] = read(group_attributes[key])
# end
# println(inputs)

# S = inputs["S"]
# ITensors.space(::SiteType"Spin") = 2*S+1
# function ITensors.op(::OpName"Sz_sp", ::SiteType"Spin", s::Index)

#     d = dim(s)
#     S = div(d-1, 2) 

#     return diagm(S:-1:-S)

# end

# function ITensors.op(::OpName"Sz2_sp", ::SiteType"Spin", s::Index)

#     d = dim(s) 
#     S = div(d-1, 2)
    
#     return diagm((S:-1:-S).^2)

# end

# function ITensors.op(::OpName"Sz", ::SiteType"Spin", s::Index)

#     d = dim(s)
#     S = div(d-1, 2)
#     tmp = diagm(S:-1:-S)

#     return ITensor(tmp, s', dag(s))

# end

# function ITensors.op(::OpName"Sz2", ::SiteType"Spin", s::Index)

#     d = dim(s) 
#     S = div(d-1, 2)
#     tmp = diagm((S:-1:-S).^2)

#     return ITensor(tmp, s', dag(s))

# end

# function ITensors.op(::OpName"S+", ::SiteType"Spin", s::Index)

#     """
#     See definition: https://easyspin.org/easyspin/documentation/spinoperators.html
#     """

#     d = dim(s) 
#     S = div(d-1, 2)
#     m = zeros(Float64, d, d)
#     Sz_vec = S:-1:-S
#     for row in 1:d
#         for col in 1:d

#             if col == row + 1

#                 m_row, m_col = Sz_vec[row], Sz_vec[col]

#                 m[row, col] = sqrt(S*(S+1) - m_row*m_col)

#             end

#         end
#     end

#     return ITensor(m, s', dag(s))

# end

# function ITensors.op(::OpName"S-", ::SiteType"Spin", s::Index)

#     d = dim(s)
#     S = div(d-1, 2)
#     m = zeros(Float64, d, d)
#     Sz_vec = S:-1:-S
#     for row in 1:d
#         for col in 1:d

#             if row == col + 1

#                 m_row, m_col = Sz_vec[row], Sz_vec[col]

#                 m[row, col] = sqrt(S*(S+1) - m_row*m_col)

#             end

#         end
#     end

#     return ITensor(m, s', dag(s))

# end

# gs = read(f, "gs", MPS)
# sites = siteinds(gs)

# greens_fns, d = get_greens_function_expectation_values(gs, sites)

# scatter(d, real(greens_fns))
# savefig("test.png")

# --------------------------------------------------------------------------------------------------------------------------------

# Opening up a beta = inf, T = 0 state from results

project = 6
job_id = 7

f = h5open("DAGS/$(project)/HDF5/$(job_id).h5", "r")
fi = h5open("DAGS/$(project)/inputs.h5", "r")
group = fi["$(job_id)"]
group_attributes = attributes(group)
inputs = Dict()
for key in keys(group_attributes)
    inputs[key] = read(group_attributes[key])
end
N, x, mg, S, lambda, energy_tol, cutoff, iters = inputs["N"], inputs["x"], inputs["mg"], inputs["S"], inputs["lambd"], inputs["et"], inputs["cutoff"], inputs["ns"]
title = "N = $(N), x = $(x), mg = $(mg), S = $(S), lambda = $(lambda)\n energy_tol = $(energy_tol), cutoff = $(cutoff), iters = $(iters), "
# println(inputs)

ITensors.space(::SiteType"Spin") = 2*S+1
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

energy = read(f, "energy")
gs = read(f, "gs", MPS)
penalty_term = read(f, "penalty_term")
total_z = read(f, "total_z")

function plot_string_antistring_statistics()

    function count_string_antistring(s)

        string_dict = Dict()
        antistring_dict = Dict()

        t = s[1]
        len = 1

        for i in 2:length(s)
            if s[i] == t
                len += 1
            else
                if t == 1
                    if haskey(string_dict, len)
                        string_dict[len] += 1
                    else
                        string_dict[len] = 1
                    end
                else 
                    if haskey(antistring_dict, len)
                        antistring_dict[len] += 1
                    else
                        antistring_dict[len] = 1
                    end
                end
                len = 1
                t = s[i]
            end
        end
        return string_dict, antistring_dict
    end

    function get_string_antistring_count_statistics(state, samples)

        string_dict_total = Dict()
        antistring_dict_total = Dict()

        l_state = length(state)

        for sample_number in 1:samples

            s = (collect(sample(state)[1:2:l_state]).-1.5).*2
            string_dict, antistring_dict = count_string_antistring(s)
            string_dict_total = mergewith(+, string_dict_total, string_dict)
            antistring_dict_total = mergewith(+, antistring_dict_total, antistring_dict)

        end

        return string_dict_total, antistring_dict_total

    end

    string_dict_total, antistring_dict_total = get_string_antistring_count_statistics(gs, 1000)

    keys_to_plot = sort(collect(union(keys(string_dict_total), keys(antistring_dict_total))))

    values1 = [get(string_dict_total, k, 0) for k in keys_to_plot]
    values2 = [get(antistring_dict_total, k, 0) for k in keys_to_plot]
    values = zeros(length(values1), 2)
    values[:, 1] = values1
    values[:, 2] = values2
    values ./= maximum(values) 

    # bar(keys_to_plot, values1, color = :black)
    # bar!(keys_to_plot, values2, alpha = 0.3, color = :yellow)
    groupedbar(keys_to_plot, values)
    savefig("Scrap_Plots/p_$(project)_id_$(job_id)_string_length_statistics.png")

end

plot_string_antistring_statistics()

# z_config = expect(gs, "Z"; sites = 1:2:length(gs))
# N = inputs["N"]
# staggering = [(-1)^(n-1) for n in 1:N]
# q_config = (z_config .+ staggering) ./2
# scatter(z_config, titlefont = 10)
# title!(title)
# ylabel!("Z expectation value")
# xlabel!("Site number")
# savefig("Scrap_Plots/p_$(project)_id_$(job_id)_zconfig.png")
# scatter(q_config, titlefont = 10)
# title!(title)
# ylabel!("Charge expectation value")
# xlabel!("Site number")
# savefig("Scrap_Plots/p_$(project)_id_$(job_id)_qconfig.png")

# sites = siteinds(gs)
# mpo = MPO(sites, "Id")
# z_config = []
# for i in 1:2:length(mpo)
#     opsum = OpSum()
#     opsum += "Z",i
#     mpo_z = MPO(opsum, sites)
#     push!(z_config, inner(mpo, mpo_z)/tr(mpo))
# end
# println(sum(z_config))

# gs = read(f, "gs", MPS)
# sites = siteinds(gs)

# greens_fns, d = get_greens_function_expectation_values(gs, sites)

# scatter(d, real(greens_fns))
# savefig("test.png")

# --------------------------------------------------------------------------------------------------------------------------------
