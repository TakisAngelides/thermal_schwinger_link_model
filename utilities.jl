function get_Hamiltonian_MPO(x::Real, mg::Real, sites::Vector{Index{Int64}})

    """
    Given a lattice, get the spin S QLM Hamiltonian in MPO representation on it (not including the penalty term for the Gauss law)
    """

    # Compute the length
    Nsites = length(sites)
    N = div(Nsites + 1, 2)
    link_dof = dim(sites[2])
    S = div(link_dof - 1, 2)

    # Prepare the Hamiltonian
    terms = OpSum()

    for i = 1:2:Nsites

        # Compute the site index
        n = Int(round((i + 1) / 2))

        # The staggered mass term
        terms += mg * sqrt(x) * (-1)^n, "Z", i

        if n < N
            # The electric energy operator
            terms += "Sz2", i + 1
            
            # The kinetic energy part
            terms += -x / (sqrt(S * (S + 1))), "S+", i, "S+", i + 1, "S-", i + 2
            terms += -x / (sqrt(S * (S + 1))), "S-", i, "S-", i + 1, "S+", i + 2
        end
    end

    return MPO(terms, sites)
end

function get_external_charges_list(dist)

    external_charges = [0 for _ in 1:N] # negative charges on odd sites using the convention Q_n = (Z_n + (-1)^n) / 2
    
    if dist == 0
        
        return external_charges

    else
    
        q1 = Int(ceil(N / 2 - dist / 2))
        q2 = q1 + dist
        if iseven(q1)
            external_charges[q1] = 1
            external_charges[q2] = -1
        else
            external_charges[q1] = -1
            external_charges[q2] = 1
        end

        return external_charges

    end

end

function get_gauge_invariant_subspace_projector_MPO(electric_field_cutoff, sites, external_charges)

    n = length(sites) # this is 2N-1 where N is the number of matter sites, i.e. we dont have a link to the right of the last matter site

    max_link_site_dim = 2 * electric_field_cutoff + 1 # this will be the bond dimension of the MPO
    links = [Index(max_link_site_dim, "Link,l=$j") for j in 1:n-1] # create the links for the MPO
    mpo = MPO(sites)
    Sz_vec = diag(op("Sz", siteinds("Spin", 1)))

    for site_idx in 1:n

        # case of matter sites
        if isodd(site_idx)

            matter_idx = div(site_idx + 1, 2) # this ranges from 1 to N and we need this to evaluate properly Q_n = ( Z_n + (-1)^n ) / 2

            if site_idx == 1

                u, d, r = prime(sites[site_idx]), dag(sites[site_idx]), links[site_idx] # indices up, down and right
                mpo[site_idx] = ITensor(Int64, u, d, r)

                for uval in 1:2
                    for dval in 1:2
                        for rval in 1:max_link_site_dim

                            # Convert index values to the physical quantities based on q_n = (Z_n + (-1)^(matter_idx))/2
                            if isodd(matter_idx)
                                if uval == 1
                                    q = 0
                                else
                                    q = -1
                                end
                            else
                                if uval == 1
                                    q = 1
                                else
                                    q = 0
                                end
                            end
                            L_n = Sz_vec[rval] # this will now range from -electric_field_cutoff to +electric_field_cutoff

                            # Check for Gauss law
                            if (L_n == q + external_charges[matter_idx]) && (uval == dval)
                                mpo[site_idx][u=>uval, d=>dval, r=>rval] = 1 # d => uval makes this a delta function on the physical legs 
                            else
                                mpo[site_idx][u=>uval, d=>dval, r=>rval] = 0
                            end

                        end
                    end
                end

            elseif site_idx == n

                u, d, l = prime(sites[site_idx]), dag(sites[site_idx]), links[site_idx-1] # indices up, down and left
                mpo[site_idx] = ITensor(Int64, u, d, l)

                for uval in 1:2
                    for dval in 1:2
                        for lval in 1:max_link_site_dim

                            # Convert index values to the physical quantities based on q_n = (Z_n + (-1)^(matter_idx))/2
                            if isodd(matter_idx)
                                if uval == 1
                                    q = 0
                                else
                                    q = -1
                                end
                            else
                                if uval == 1
                                    q = 1
                                else
                                    q = 0
                                end
                            end
                            L_nminus1 = Sz_vec[lval] # this will now range from -electric_field_cutoff to +electric_field_cutoff

                            # Check for Gauss law
                            if (-L_nminus1 == q + external_charges[matter_idx]) && (uval == dval)
                                mpo[site_idx][u=>uval, d=>dval, l=>lval] = 1
                            else
                                mpo[site_idx][u=>uval, d=>dval, l=>lval] = 0
                            end

                        end
                    end
                end

            else

                u, d, l, r = prime(sites[site_idx]), dag(sites[site_idx]), links[site_idx-1], links[site_idx] # indices up, down and left
                mpo[site_idx] = ITensor(Int64, u, d, l, r)

                for uval in 1:2
                    for dval in 1:2
                        for lval in 1:max_link_site_dim
                            for rval in 1:max_link_site_dim

                                # Convert index values to the physical quantities based on q_n = (Z_n + (-1)^(matter_idx))/2
                                if isodd(matter_idx)
                                    if uval == 1
                                        q = 0
                                    else
                                        q = -1
                                    end
                                else
                                    if uval == 1
                                        q = 1
                                    else
                                        q = 0
                                    end
                                end
                                L_nminus1 = Sz_vec[lval] # this will now range from -electric_field_cutoff to +electric_field_cutoff
                                L_n = Sz_vec[rval]

                                # Check for Gauss law L_n - L_n-1 = q_n + Q_n
                                if (L_n == L_nminus1 + q + external_charges[matter_idx]) && (uval == dval)
                                    mpo[site_idx][u=>uval, d=>dval, l=>lval, r=>rval] = 1
                                else
                                    mpo[site_idx][u=>uval, d=>dval, l=>lval, r=>rval] = 0
                                end

                            end
                        end
                    end
                end

            end

            # case of gauge fields
        else

            u, d, l, r = prime(sites[site_idx]), dag(sites[site_idx]), links[site_idx-1], links[site_idx] # indices up, down and left
            mpo[site_idx] = ITensor(Int64, u, d, l, r)

            for uval in 1:max_link_site_dim
                for dval in 1:max_link_site_dim
                    for lval in 1:max_link_site_dim
                        for rval in 1:max_link_site_dim
                            if (uval == dval) && (lval == rval) && (uval == lval)
                                mpo[site_idx][u=>uval, d=>dval, l=>lval, r=>rval] = 1 # this is just delta on the physical legs and a delta on the virtual legs and a delta between one virtual and one physical leg
                            else
                                mpo[site_idx][u=>uval, d=>dval, l=>lval, r=>rval] = 0
                            end
                        end
                    end
                end
            end

        end

    end

    return mpo

end

function mixed_mpo_to_matrix(mpo)

    mixed_mpo_dims = []
    for element in siteinds(mpo)
        push!(mixed_mpo_dims, Int(sqrt(dim(element))))
    end
    final_dim = prod(mixed_mpo_dims)
    a = contract(mpo)
    a = Array(a, inds(a; :plev => 1)..., inds(a; :plev => 0)...)
    a = reshape(a, final_dim, final_dim)

    return a

end

function get_sites(N::Int)
    
    Ntotal = 2 * N - 1
    s = Vector{Index{Int64}}(undef, Ntotal)
    # Make an array of 'site' indices alternating between spin 1/2 and spin 1
    for i = 1:Ntotal
        if isodd(i)
            s[i] = siteind("S=1/2", i)
        else
            s[i] = siteind("Spin", i)
        end
    end
    return s
end

function get_charge_config_rho(rho, sites)

    """
    This is always assuming the rho is an MPO with a matter site, then a link site and so on ending on a matter site.
    """

    len = length(sites)
    N = div(len + 1, 2)
    charge_config = []
    for idx in 1:N
        opsum = OpSum()
        opsum += 0.5*(-1)^idx, "Id", 1
        opsum += 0.5, "Z", 2 * idx - 1
        mpo = MPO(opsum, sites)
        push!(charge_config, real(tr(apply(mpo, rho))))
    end

    return charge_config

end

function get_electric_field_config(rho, sites)

    N = div(length(sites) + 1, 2)
    electric_field_config = []
    for idx in 1:N-1
    
        opsum = OpSum()
        opsum += "Sz",2 * idx - 1
        mpo = MPO(opsum, sites)
        push!(electric_field_config, real(tr(apply(rho, mpo))))
    
    end

    return electric_field_config

end

function get_which_canonical_form(mps)

    N = length(mps)
    canonical_form::Array{String} = []

    for site in 1:N

        mps_site = mps[site]
        
        a = mps_site
        adag = dag(mps_site)
        if site != 1
            adag_idx = commonind(a, mps[site-1])
            replaceind!(adag, adag_idx, prime(adag_idx))
        end
        res = a*adag
        inds_res = inds(res)
        res = ITensors.Array(res, inds_res...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_right = isapprox(res, I(l))

        mps_site = mps[site]
        a = mps_site
        adag = dag(mps_site)
        if site != N
            adag_idx = commonind(a, mps[site+1])
            replaceind!(adag, adag_idx, prime(adag_idx))
        end
        res = a*adag
        inds_res = inds(res)
        res = ITensors.Array(res, inds_res...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_left = isapprox(res, I(l))

        if is_left
            if is_right 
                push!(canonical_form, "L/R")
            else
                push!(canonical_form, "L")
            end
        elseif is_right
            push!(canonical_form, "R")
        else
            push!(canonical_form, "N")
        end

    end

    return canonical_form

end

function get_odd_gates(sites, a, x, mg, S)

    l = length(sites)

    gates = []

    for (i_idx, i) in enumerate(1:4:l-2)

        physical_index = div(i + 1, 2)

        if i_idx == length(1:4:l-2)

            h = (-x/sqrt(S*(S+1))) * op("S+", sites[i]) * op("S+", sites[i+1]) * op("S-", sites[i+2])
            h += (-x/sqrt(S*(S+1))) * op("S-", sites[i]) * op("S-", sites[i+1]) * op("S+", sites[i+2]) 
            h += (mg*sqrt(x)*(-1)^(physical_index)) * op("Z", sites[i]) * op("I", sites[i+1]) * op("I", sites[i+2]) 
            h += op("I", sites[i]) * op("Sz2", sites[i+1]) * op("I", sites[i+2])
            h += (mg*sqrt(x)*(-1)^(physical_index + 1)) * op("I", sites[i]) * op("I", sites[i+1]) * op("Z", sites[i+2])  
            G = exp(a*h)    

        else

            h = (-x/sqrt(S*(S+1))) * op("S+", sites[i]) * op("S+", sites[i+1]) * op("S-", sites[i+2])
            h += (-x/sqrt(S*(S+1))) * op("S-", sites[i]) * op("S-", sites[i+1]) * op("S+", sites[i+2]) 
            h += (mg*sqrt(x)*(-1)^(physical_index)) * op("Z", sites[i]) * op("I", sites[i+1]) * op("I", sites[i+2]) 
            h += op("I", sites[i]) * op("Sz2", sites[i+1]) * op("I", sites[i+2]) 
            G = exp(a*h)     

        end

        push!(gates, G)

    end

    return gates

end

function get_even_gates(sites, a, x, mg, S)

    l = length(sites)

    gates = []

    for i in 3:4:l-2

        physical_index = div(i + 1, 2)

        h = (-x/sqrt(S*(S+1))) * op("S+", sites[i]) * op("S+", sites[i+1]) * op("S-", sites[i+2])
        h += (-x/sqrt(S*(S+1))) * op("S-", sites[i]) * op("S-", sites[i+1]) * op("S+", sites[i+2]) 
        h += (mg*sqrt(x)*(-1)^(physical_index)) * op("Z", sites[i]) * op("I", sites[i+1]) * op("I", sites[i+2]) 
        h += op("I", sites[i]) * op("Sz2", sites[i+1]) * op("I", sites[i+2]) 
        G = exp(a*h)     

        push!(gates, G)

    end

    return gates

end

function get_odd_gates_opsum(sites, x, mg, S)

    l = length(sites)

    h = OpSum()

    for (i_idx, i) in enumerate(1:4:l-2)

        physical_index = div(i + 1, 2)

        if i_idx == length(1:4:l-2)

            h += (-x/sqrt(S*(S+1))),"S+",i,"S+",i+1,"S-",i+2
            h += (-x/sqrt(S*(S+1))),"S-",i,"S-",i+1,"S+",i+2 
            h += (mg*sqrt(x)*(-1)^(physical_index)),"Z",i 
            h += "Sz2",i+1 
            h += (mg*sqrt(x)*(-1)^(physical_index + 1)),"Z",i+2   

        else

            h += (-x/sqrt(S*(S+1))),"S+",i,"S+",i+1,"S-",i+2
            h += (-x/sqrt(S*(S+1))),"S-",i,"S-",i+1,"S+",i+2 
            h += (mg*sqrt(x)*(-1)^(physical_index)),"Z",i 
            h += "Sz2",i+1 

        end

    end

    return h

end

function get_even_gates_opsum(sites, x, mg, S)

    l = length(sites)

    h = OpSum()

    for i in 3:4:l-2

        physical_index = div(i + 1, 2)

        h += (-x/sqrt(S*(S+1))),"S+",i,"S+",i+1,"S-",i+2
        h += (-x/sqrt(S*(S+1))),"S-",i,"S-",i+1,"S+",i+2 
        h += (mg*sqrt(x)*(-1)^(physical_index)),"Z",i 
        h += "Sz2",i+1 

    end

    return h

end

function get_exp_Hz_MPO(x, mg, sites, b)

    """
    The 'a' input here will take care of the beta prefactors
    """

    n = length(sites)
    links = [Index(1, "Link,l=$i") for i in 1:n-1]
    mpo = MPO(sites)
    link_dof = dim(sites[2])
    Sz2_sp = op("Sz2_sp", sites[2])

    for site_idx in 1:n

        # case of matter sites
        if isodd(site_idx)

            matter_idx = div(site_idx + 1, 2) # this ranges from 1 to N

            if site_idx == 1

                u, d, r = prime(sites[site_idx]), dag(sites[site_idx]), links[site_idx] # indices up, down and right
                mpo[site_idx] = ITensor(Float64, u, d, r)

                mpo[site_idx][u=>1, d=>1, r=>1] = exp(b*sqrt(x)*mg*(-1)^(matter_idx))
                mpo[site_idx][u=>1, d=>2, r=>1] = 0.0 
                mpo[site_idx][u=>2, d=>1, r=>1] = 0.0
                mpo[site_idx][u=>2, d=>2, r=>1] = exp(-b*sqrt(x)*mg*(-1)^(matter_idx))

            elseif site_idx == n

                u, d, l = prime(sites[site_idx]), dag(sites[site_idx]), links[site_idx-1] # indices up, down and left
                mpo[site_idx] = ITensor(Float64, u, d, l)

                mpo[site_idx][u=>1, d=>1, l=>1] = exp(b*sqrt(x)*mg*(-1)^(matter_idx))
                mpo[site_idx][u=>1, d=>2, l=>1] = 0.0
                mpo[site_idx][u=>2, d=>1, l=>1] = 0.0
                mpo[site_idx][u=>2, d=>2, l=>1] = exp(-b*sqrt(x)*mg*(-1)^(matter_idx))

            else

                u, d, l, r = prime(sites[site_idx]), dag(sites[site_idx]), links[site_idx-1], links[site_idx] # indices up, down and left
                mpo[site_idx] = ITensor(Float64, u, d, l, r)

                mpo[site_idx][u=>1, d=>1, l=>1, r=>1] = exp(b*sqrt(x)*mg*(-1)^(matter_idx))
                mpo[site_idx][u=>1, d=>2, l=>1, r=>1] = 0.0
                mpo[site_idx][u=>2, d=>1, l=>1, r=>1] = 0.0
                mpo[site_idx][u=>2, d=>2, l=>1, r=>1] = exp(-b*sqrt(x)*mg*(-1)^(matter_idx))

            end

        # case of gauge fields
        else

            u, d, l, r = prime(sites[site_idx]), dag(sites[site_idx]), links[site_idx-1], links[site_idx] # indices up, down and left
            mpo[site_idx] = ITensor(Float64, u, d, l, r)

            for uval in 1:link_dof
                for dval in 1:link_dof
                    
                    if uval == dval
                        mpo[site_idx][u=>uval, d=>dval, l=>1, r=>1] = exp(b*Sz2_sp[uval, dval])
                    else
                        mpo[site_idx][u=>uval, d=>dval, l=>1, r=>1] = 0.0
                    end
                           
                end
            end

        end

    end

    return mpo

end

function apply_odd_gates!(odd, mpo; cutoff = 0, maxdim = 1000)

    l = length(mpo)

    for (n_idx, n) in enumerate(1:4:l-2) # n is the left most site the gate acts on, the gates are always 3 site in span

        gate = odd[n_idx]

        t = replaceprime(gate''*prime(prod(mpo[n:n+2]); :tags => "Site")*gate, 3, 1)

        # Separate the t tensor into individual MPO site tensors
        for idx in n:n+1
            t_indices = inds(t)
            U_indices = filter(i -> hastags(i, "n=$(idx)") || hastags(i, "Link,l=$(idx-1)"), t_indices)
            mpo[idx], S, V = ITensors.svd(t, U_indices; cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)")
            t = S*V
        end
        mpo[n+2] = t
        
        if n+4 <= l
            
            # Extra SVD as required by ATD DMRG
            t = mpo[n+2]*mpo[n+3]
            mpo[n+2], S, V = ITensors.svd(t, uniqueinds(t, mpo[n+3]); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(n+2)", righttags = "Link,l=$(n+2)")
            mpo[n+3] = S*V

            t = mpo[n+3]*mpo[n+4]
            mpo[n+3], S, V = ITensors.svd(t, uniqueinds(t, mpo[n+4]); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(n+3)", righttags = "Link,l=$(n+3)")
            mpo[n+4] = S*V

        end

    end

end

function apply_even_gates!(even, mpo; cutoff = 0, maxdim = 1000)

    l = length(mpo)

    t = mpo[l-1]*mpo[l]
    U, S, mpo[l] = ITensors.svd(t, uniqueinds(t, mpo[l]); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(l-1)", righttags = "Link,l=$(l-1)")
    mpo[l-1] = U*S

    t = mpo[l-2]*mpo[l-1]
    U, S, mpo[l-1] = ITensors.svd(t, uniqueinds(t, mpo[l-1]); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(l-2)", righttags = "Link,l=$(l-2)")
    mpo[l-2] = U*S

    for (n_idx, n) in enumerate(l-4:-4:3) # n is the left most site the gate acts on, the gates are always 3 site in span

        gate = even[end-n_idx+1]

        t = replaceprime(gate''*prime(prod(mpo[n:n+2]); :tags => "Site")*gate, 3, 1)

        # Separate the t tensor into individual MPO site tensors
        for idx in n+2:-1:n+1
            t_indices = inds(t)
            V_indices = filter(i -> hastags(i, "n=$(idx)") || hastags(i, "Link,l=$(idx)"), t_indices)
            U, S, mpo[idx] = ITensors.svd(t, uniqueinds(t_indices, V_indices); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(idx-1)", righttags = "Link,l=$(idx-1)")
            t = U*S
        end
        mpo[n] = t
            
        # Extra SVD as required by ATD DMRG
        t = mpo[n-1]*mpo[n]
        U, S, mpo[n] = ITensors.svd(t, uniqueinds(t, mpo[n]); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(n-1)", righttags = "Link,l=$(n-1)")
        mpo[n-1] = U*S

        t = mpo[n-2]*mpo[n-1]
        U, S, mpo[n-1] = ITensors.svd(t, uniqueinds(t, mpo[n-1]); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(n-2)", righttags = "Link,l=$(n-2)")
        mpo[n-2] = U*S

    end

end

function get_product_state_MPO(sites, state)

    """
    Prepare the mpo = |state><state| where |state> is a basis state given as a list of integers 1 and 2 e.g. state = [1,1,2,1] which would be the state |0010>
    """

    N = length(sites)
    mpo = MPO(sites)

    for i=1:N

        if i == 1 || i == N

            s, sp = inds(mpo[i]; :tags => "Site")
            l = inds(mpo[i]; :tags => "Link")[1]
            mpo[i][s => state[i], sp => state[i], l => 1] = 1.0

        else

            s, sp = inds(mpo[i]; :tags => "Site")
            l1, l2 = inds(mpo[i]; :tags => "Link")
            mpo[i][s => state[i], sp => state[i], l1 => 1, l2 => 1] = 1.0

        end

    end

    return mpo

end

function get_MPO_site_canonical_form(mpo_site, which_site, mpo_site_index)

    mpo_site_dag = mpo_site'
    noprime!(mpo_site_dag; :plev => 2)
    mpo_site_dag = dag(mpo_site_dag'; :tags => "Site")
    tmp = mpo_site_dag * mpo_site

    # Checking for LCF
    if which_site == "first"

        tmp_l = tmp * dag(delta(inds(tmp; :tags => "Site")))
        res = ITensors.Array(tmp_l, inds(tmp_l)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_left = isapprox(res, I(l))

    elseif which_site == "last"

        tmp_l = tmp * dag(delta(inds(tmp; :tags => "Site")))
        tmp_l = tmp_l * dag(delta(inds(tmp_l; :tags => "Link")))
        res = ITensors.Array(tmp_l, inds(tmp_l)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_left = isapprox(res, I(l))
        
    else

        tmp_l = tmp * dag(delta(inds(tmp; :tags => "Site")))
        tmp_l = tmp_l * dag(delta(inds(tmp_l; :tags => "Link,l=$(mpo_site_index-1)")))
        res = ITensors.Array(tmp_l, inds(tmp_l)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_left = isapprox(res, I(l))

    end

    # Checking for RCF
    if which_site == "first"

        tmp = tmp * dag(delta(inds(tmp; :tags => "Site")))
        tmp = tmp * dag(delta(inds(tmp; :tags => "Link")))
        res = ITensors.Array(tmp, inds(tmp)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_right = isapprox(res, I(l))

    elseif which_site == "last"

        tmp = tmp * dag(delta(inds(tmp; :tags => "Site")))
        res = ITensors.Array(tmp, inds(tmp)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_right = isapprox(res, I(l))
        
    else

        tmp = tmp * dag(delta(inds(tmp; :tags => "Site")))
        tmp = tmp * dag(delta(inds(tmp; :tags => "Link,l=$(mpo_site_index)")))
        res = ITensors.Array(tmp, inds(tmp)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_right = isapprox(res, I(l))

    end

    if is_left
        if is_right 
            return "L/R"
        else
            return "L"
        end
    elseif is_right
        return "R"
    else
        return "N"
    end

end

function get_MPO_canonical_form(mpo)
    
    res = String[]
    n = length(mpo)
    for (mpo_site_index, mpo_site) in enumerate(mpo)
        if mpo_site_index == 1
            push!(res, get_MPO_site_canonical_form(mpo_site, "first", mpo_site_index))
        elseif mpo_site_index == n
            push!(res, get_MPO_site_canonical_form(mpo_site, "last", mpo_site_index))
        else
            push!(res, get_MPO_site_canonical_form(mpo_site, "none", mpo_site_index))
        end
    end
    return res
end

function is_MPO_hermitian(mpo; tol = 1e-14)

    return norm(mpo - dag(swapprime(mpo, 0, 1))) < tol

end

function is_MPO_positive(mpo, sites; tol = 1e-14)

    """
    We check this by using DMRG to find the ground state energy of the "density matrix" MPO
    """

    dmrg_tol = tol
    cutoff = tol
    psi0 = randomMPS(sites)
    obs = DMRGObserver(;energy_tol = dmrg_tol)
    nsweeps = 100
    dmrg_energy, _ = dmrg(mpo, psi0; nsweeps, cutoff, observer=obs, outputlevel = 0)

    return abs(dmrg_energy) < tol

end

function mpo_to_matrix(mpo)

    n = length(mpo)
    a = contract(mpo)
    a = Array(a, inds(a; :plev => 1)..., inds(a; :plev => 0)...)
    a = reshape(a, 2^n, 2^n)

    return a

end

function my_kron(A, B)
    
    m, n = size(A)
    p, q = size(B)

    C = zeros(ComplexF64, m * p, n * q)

    for i in 1:p
        for j in 1:q
            C[(i-1)*m+1 : i*m, (j-1)*n+1 : j*n] = A * B[i, j]
        end
    end

    return C
end

function get_op(ops, positions, N; reverse_flag = true)

    op_dict = Dict("X" => sparse([0 1; 1 0]), "Y" => sparse([0 -1im; 1im 0]), "Z" => sparse([1 0; 0 -1]))
    zipped = TupleTools.sort(Tuple(zip(1:length(ops), positions, ops)); by = x -> x[2])
    old_positions = [element[2] for element in zipped] 
    old_ops = [element[3] for element in zipped]

    positions = []
    ops = []

    if length(Set(old_positions)) != length(old_positions) # case where we have duplicate positions
        
        flag = false

        for (idx, pos) in enumerate(old_positions)

            if flag

                flag = false
                continue

            end

            if idx != length(old_positions)

                if pos != old_positions[idx+1]

                    push!(positions, pos)
                    push!(ops, op_dict[old_ops[idx]])

                else

                    push!(positions, pos)
                    push!(ops, op_dict[old_ops[idx]]*op_dict[old_ops[idx+1]])
                    flag = true

                end

            else

                push!(positions, pos)
                push!(ops, op_dict[old_ops[idx]])

            end

        end

    else

        for (idx, pos) in enumerate(old_positions)

            push!(positions, pos)
            push!(ops, op_dict[old_ops[idx]])
        
        end

    end

    eye(n) = sparse(I, n, n)

    res = eye(1)

    for (i, pos) in enumerate(positions)

        if i == 1
            how_many_I_before = pos-1
        else
            how_many_I_before = pos - positions[i-1] - 1
        end

        pos = positions[i]
        op = ops[i]
    
        if reverse_flag
            res = my_kron(res, eye(2^how_many_I_before))
            res = my_kron(res, op)
        else
            res = kron(res, eye(2^how_many_I_before))
            res = kron(res, op)
        end

    end

    if reverse_flag
        res = my_kron(res, eye(2^(N - positions[end])))
    else
        res = kron(res, eye(2^(N - positions[end])))
    end

    return res

end

function hermitian_conjugate_mpo(mpo)

    return dag(swapprime(mpo, 0, 1))

end

function transpose_mpo(mpo)

    return swapprime(dag(conj(mpo)), 0 => 1)

end

function get_entanglement_entropy_mpo(rho, trace_indices, sites; tol = 1e-12)

    N = length(sites) - length(trace_indices)

    tmp = []
    for (idx, element) in enumerate(rho)
        if idx in trace_indices
            push!(tmp, element * delta(dag(sites[idx]'), sites[idx]))
        else
            push!(tmp, element)
        end
    end

    a = contract(tmp)
    a = Array(a, inds(a; :plev => 1)..., inds(a; :plev => 0)...)
    a = reshape(a, 2^N, 2^N)

    evals, _ = eigen(a)

    ee = sum(-real(eval)*log(real(eval)) for eval in evals if real(eval) >= tol)

    return ee

end

function get_entanglement_entropy_matrix(N, rho_m, keep_indices; tol = 1e-12)

    a = partial_trace(rho_m, keep_indices)

    a = reshape(a, 2^(div(N, 2)), 2^(div(N, 2)))

    evals, _ = eigen(a)

    ee2 = sum(-real(eval)*log(real(eval)) for eval in evals if real(eval) >= tol)

    return ee2

end

function check_zeroq(n, N)

    """
    Given an integer n for N qubits, it converts the integer n to a bit representation and checks whether
    the bit representation has equal number of 0s and 1s
    """

    return sum((digits(n, base=2, pad = N).*2).-1) == 0
end

function project_zeroq(M)

    """
    Given a matrix M, it projects it to the zero total charge subspace.
    The basis of M is assumed to be going from 1 to 2^N for N qubits, i.e. the first state is |00...00><00..00| and the last |11..11><11..11|
    so if the matrix is already somehow projected, e.g. projected into the Gauge invariant subspace this function probably will not work on that.
    """

    nrow, ncol = size(M)
    n = Int(log2(nrow))
    new_nrow = binomial(n, div(n, 2))
    res = zeros(ComplexF64, new_nrow, new_nrow)

    row_count = 0
    for row in 1:nrow

        if !(check_zeroq(row-1, n))
            continue
        else
            row_count += 1
        end
        # println(row_count, " ", digits(row-1, base = 2, pad = 4))
        col_count = 0
        for col in 1:ncol

            if !(check_zeroq(col-1, n))
                continue
            else
                col_count += 1
                # println("row_count, col_count: ", (row_count, col_count), ", row, col: ", (digits(row-1, base = 2, pad = n), digits(col-1, base = 2, pad = n)))
                res[row_count, col_count] = M[row, col]
            end

        end

    end

    return res

end

function get_entanglement_entropy_reduced_matrix(N, rho_m; tol = 1e-12)

    # If you have 4 qubits 1234 it computes the von Neumann entropy by partitioning 
    # the system in half 12 34 so it always assumes you have even number of sites

    dimres = 2^(div(N, 2))
    
    res = zeros(ComplexF64, dimres, dimres)

    zero_q_list = [join(digits(i-1, base = 2, pad = N)) for i in 1:2^N if check_zeroq(i-1, N)] 

    for row in 1:dimres
        for col in 1:dimres

            for trace_idx in 1:dimres 

                bigrow = join(vcat(digits(row-1, base=2, pad = div(N,2)), digits(trace_idx-1, base=2, pad = div(N,2)))) # this is bin(row))bin(trace_idx)
                bigcol = join(vcat(digits(col-1, base=2, pad = div(N,2)), digits(trace_idx-1, base=2, pad = div(N,2)))) # this is bin(col))bin(trace_idx)

                if !(bigrow in zero_q_list) || !(bigcol in zero_q_list)
                    continue
                else
                    bigrow_idx = findfirst(x -> x == bigrow, zero_q_list)
                    bigcol_idx = findfirst(x -> x == bigcol, zero_q_list)
                    # println("r_a = ", row, ", c_a = ", col, " ", bigrow, " ", bigcol, " r_a r_b = ", bigrow_idx, ", c_a c_b = ", bigcol_idx, " value = ", real(rho_m[bigrow_idx, bigcol_idx]))
                    # println(bigrow_idx, " ", bigcol_idx)
                    res[row, col] += rho_m[bigrow_idx, bigcol_idx]
                end

            end
                
        end

    end

    evals, _ = eigen(res)

    ee2 = sum(-real(eval)*log(real(eval)) for eval in evals if real(eval) >= tol)

    return ee2

end

function swap(i, j, N)

    res = sparse(I, 2^N, 2^N)

    idx1 = min(i, j)
    idx2 = max(i, j)

    local_swap = sparse([[1 0 0 0]; [0 0 1 0]; [0 1 0 0]; [0 0 0 1]])
    full_swap(idx) = kron(sparse(I, 2^(idx-1), 2^(idx-1)), kron(local_swap, sparse(I, 2^(N-idx-1), 2^(N-idx-1))))

    for k in idx1:idx2-1

        res *= full_swap(k)

    end

    if idx2-idx1 > 1

        for k in reverse(idx1:idx2-2)

            res *= full_swap(k)

        end

    end

    return res

end

function decimal_to_padded_binary_list(decimal, bit_length)
    binary_list = Int[]

    while decimal > 0 || length(binary_list) < bit_length
        pushfirst!(binary_list, (decimal % 2) + 1)
        decimal = div(decimal, 2)
    end

    # Pad with leading zeros if needed
    while length(binary_list) < bit_length
        pushfirst!(binary_list, 0)
    end

    return binary_list 
end

function get_gauss_law_MPO(matter_idx, sites, external_charges)

    """
    n here is the matter site
    """

    len = length(sites)
    N = div(len + 1, 2)
    n = 2 * matter_idx - 1

    opsum = OpSum()
    if matter_idx == 1
        opsum += "Sz", n + 1
        opsum -= 0.5, "Z", n
        opsum -= 0.5 * (-1)^(matter_idx), "Id", n
        opsum -= external_charges[matter_idx], "Id", n
    elseif matter_idx == N
        opsum -= "Sz", n - 1
        opsum -= 0.5, "Z", n
        opsum -= 0.5 * (-1)^(matter_idx), "Id", n
        opsum -= external_charges[matter_idx], "Id", n
    else
        opsum += "Sz", n + 1
        opsum -= "Sz", n - 1
        opsum -= 0.5, "Z", n
        opsum -= 0.5 * (-1)^(matter_idx), "Id", n
        opsum -= external_charges[matter_idx], "Id", n
    end

    return MPO(opsum, sites)

end

function get_gauss_law_squared_sum_MPO(sites, external_charges)

    p = get_gauss_law_MPO(1, sites, external_charges)
    pdag = hermitian_conjugate_mpo(p)
    pdagp = apply(pdag, p)
    sum = pdagp

    for matter_idx in 2:N
    
        p = get_gauss_law_MPO(matter_idx, sites, external_charges)
        pdag = hermitian_conjugate_mpo(p)
        pdagp = apply(pdag, p)
        sum += pdagp
    
    end

    return sum

end

function get_gauss_law_squared_sum_opsum(N, external_charges)

    opsum = OpSum()

    for i in 1:2:2*N-1 # iterate over the links

        matter_idx = div(i+1, 2)

        if matter_idx == 1

            opsum += "Sz2",i+1
            opsum += -1,"Z",i,"Sz",i+1
            opsum += -(-1)^matter_idx,"Sz",i+1
            opsum += -2*external_charges[matter_idx],"Sz",i+1

            opsum += 0.25,"Z",i,"Z",i
            opsum += 0.5*(-1)^matter_idx,"Z",i
            opsum += external_charges[matter_idx],"Z",i

            opsum += 0.25,"Id",1
            opsum += external_charges[matter_idx]*(-1)^matter_idx,"Id",1

            opsum += external_charges[matter_idx]^2,"Id",1

        elseif matter_idx == N

            opsum += "Sz2",i-1
            opsum += "Sz",i-1,"Z",i
            opsum += (-1)^matter_idx,"Sz",i-1
            opsum += 2*external_charges[matter_idx],"Sz",i-1

            opsum += 0.25,"Z",i,"Z",i
            opsum += 0.5*(-1)^matter_idx,"Z",i
            opsum += external_charges[matter_idx],"Z",i

            opsum += 0.25,"Id",1
            opsum += external_charges[matter_idx]*(-1)^matter_idx,"Id",1

            opsum += external_charges[matter_idx]^2,"Id",1
        else

            opsum += "Sz2",i+1
            opsum += -2,"Sz",i-1,"Sz",i+1
            opsum += -1,"Z",i,"Sz",i+1
            opsum += -(-1)^matter_idx,"Sz",i+1
            opsum += -2*external_charges[matter_idx],"Sz",i+1

            opsum += "Sz2",i-1
            opsum += "Sz",i-1,"Z",i
            opsum += (-1)^matter_idx,"Sz",i-1
            opsum += 2*external_charges[matter_idx],"Sz",i-1

            opsum += 0.25,"Z",i,"Z",i
            opsum += 0.5*(-1)^matter_idx,"Z",i
            opsum += external_charges[matter_idx],"Z",i

            opsum += 0.25,"Id",1
            opsum += external_charges[matter_idx]*(-1)^matter_idx,"Id",1

            opsum += external_charges[matter_idx]^2,"Id",1
        
        end

    end

    return opsum

end

function get_Hamiltonian_opsum(x::Real, mg::Real, sites::Vector{Index{Int64}})

    """
    Given a lattice, get the spin S QLM Hamiltonian in MPO representation on it (not including the penalty term for the Gauss law)
    """

    # Compute the length
    Nsites = length(sites)
    N = div(Nsites + 1, 2)
    link_dof = dim(sites[2])
    S = div(link_dof - 1, 2)

    # Prepare the Hamiltonian
    terms = OpSum()

    for i = 1:2:Nsites

        # Compute the site index
        n = Int(round((i + 1) / 2))

        # The staggered mass term
        terms += mg * sqrt(x) * (-1)^n, "Z", i

        if n < N
            # The electric energy operator
            terms += "Sz2", i + 1
            
            # The kinetic energy part
            terms += -x / (sqrt(S * (S + 1))), "S+", i, "S+", i + 1, "S-", i + 2
            terms += -x / (sqrt(S * (S + 1))), "S-", i, "S-", i + 1, "S+", i + 2
        end
    end

    return terms
end

function get_greens_function_opsum(i, j)

    """
    This function returns the Green's function operator (as an OpSum):

    G(i-j) = S+_i (S+_i+1 Z_i+2 ... Z_j-2 S+_j-1) S-_j + H.c. where i is a matter site

    i+1 is a link site, i+2 is a matter site and so on

    Note: i and j need to be matter sites 

    """

    res1 = OpSum()

    res1 += "S+",i

    for k in i+1:2:j-1
        res1 *= "S+",k
        if k != j-1
            res1 *= "Z",k+1
        end
    end

    res1 *= "S-",j

    res2 = OpSum()

    res2 += "S-",i

    for k in i+1:2:j-1
        res2 *= "S-",k
        if k != j-1
            res2 *= "Z",k+1
        end
    end

    res2 *= "S+",j

    return res1 + res2

end

function get_greens_function_expectation_values(gs::MPS, sites)

    n = length(gs)
    left_indices = 1:2:div(n+1, 2)
    right_indices = [n - left_idx + 1 for left_idx in left_indices] 
    res = []
    idxs = []
    for left_idx in left_indices
        for right_idx in right_indices
            push!(idxs, [left_idx, right_idx, right_idx - left_idx])
            push!(res, inner(gs', MPO(get_greens_function_opsum(left_idx, right_idx), sites), gs))
        end
    end

    return res, idxs

end

function get_greens_function_expectation_values(rho::MPO, sites)

    n = length(rho)
    left_indices = 1:2:div(n+1, 2)
    right_indices = [n - left_idx + 1 for left_idx in left_indices] 
    res = []
    idxs = []
    for left_idx in left_indices
        for right_idx in right_indices
            push!(idxs, [left_idx, right_idx, right_idx - left_idx])
            push!(res, inner(rho, MPO(get_greens_function_opsum(left_idx, right_idx), sites)))
        end
    end

    return res, idxs

end
