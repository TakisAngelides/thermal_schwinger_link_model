using ITensors
using HDF5
using Plots
using Statistics

# Include your utilities script
include("/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/utilities.jl")

function save_distribution_to_hdf5(h5file, N, x, mg, S, beta, total_dict, num_samples)
    group_name = "/N_$(N)/x_$(x)/mg_$(mg)/S_$(S)/beta_$(beta)"
    
    # Create the group if it doesn't exist
    g = haskey(h5file, group_name) ? h5file[group_name] : create_group(h5file, group_name)
    
    # Convert the keys and values to specific types
    x_plot = Int64[convert(Int64, k) for k in keys(total_dict)]
    y_plot = Float64[convert(Float64, v) for v in values(total_dict)]
    
    # Save the data
    write(g, "x_plot", x_plot)
    write(g, "y_plot", y_plot)
end

function collect_distributions(dag_project_number, num_samples)
    # Open the inputs HDF5 file
    inputs_f = h5open("/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/DAGS/$(dag_project_number)/inputs.h5", "r")

    # Open the HDF5 file for writing the string length distributions
    h5out = h5open("/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/DAGS/$(dag_project_number)/string_length_distributions.h5", "w")
    max_x = 0

    # Iterate through each file in the HDF5 directory
    for file_name in readdir("/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/DAGS/$(dag_project_number)/HDF5")
        println(file_name)

        f = h5open("/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/DAGS/$(dag_project_number)/HDF5/$(file_name)", "r")
        group = inputs_f["$(file_name[1:end-3])"]
        group_attributes = attributes(group)
        inputs = Dict()
        for key in keys(group_attributes)
            inputs[key] = read(group_attributes[key])
        end

        N, x, mg, S = inputs["N"], inputs["x"], inputs["mg"], inputs["S"]

        if haskey(inputs, "ns")
            state = read(f, "gs", MPS)
            total_dict = get_distribution_of_lengths(state, num_samples)
            total_dict = Dict{Int, Float64}(total_dict)  # Convert to specific types
            max_x = max(max_x, maximum(collect(keys(total_dict))))
            save_distribution_to_hdf5(h5out, N, x, mg, S, "inf", total_dict, num_samples)

        else
            beta_list = read(f, "beta_list")
            for key in keys(f)
                if occursin("rho_", key)
                    beta = round(beta_list[findfirst(x -> x == parse(Int, (split(key, "_")[end])), inputs["ssl"])], digits = 5)
                    state = read(f, key, MPO)
                    total_dict = get_distribution_of_lengths(state, num_samples)
                    total_dict = Dict{Int, Float64}(total_dict)  # Convert to specific types
                    max_x = max(max_x, maximum(collect(keys(total_dict))))
                    save_distribution_to_hdf5(h5out, N, x, mg, S, beta, total_dict, num_samples)
                end
            end
        end
        close(f)  # Close the HDF5 file after processing
    end

    close(inputs_f)  # Close the inputs HDF5 file
    close(h5out)  # Close the output HDF5 file

    return max_x
end

function recalculate_max_x(hdf5_file_path::String)
    max_x = 0

    # Open the HDF5 file for reading
    h5file = h5open(hdf5_file_path, "r")
    
    # Recursive function to find all x_plot datasets
    function find_max_x(group, max_x)
        for key in keys(group)
            item = group[key]
            if isa(item, HDF5.Group)
                max_x = find_max_x(item, max_x)  # Recursively check nested groups
            elseif isa(item, HDF5.Dataset) && occursin("x_plot", key)
                x_values = read(item)
                max_x = max(max_x, maximum(x_values))
            end
        end
        return max_x
    end

    # Start the recursive search from the root group
    max_x = find_max_x(h5file, max_x)

    close(h5file)  # Close the HDF5 file

    return max_x
end

function plot_distributions_with_fixed_x(max_x, dag_project_number, num_samples)
    
    # Open the HDF5 file for reading
    h5file_path = "/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/DAGS/$(dag_project_number)/string_length_distributions.h5"
    h5file = h5open(h5file_path, "r")

    # Create a directory for the plots if it doesn't exist
    output_dir = "/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/DAGS/$(dag_project_number)/Plots/String_length_distributions"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Iterate through all groups and datasets in the HDF5 file
    for group_name in keys(h5file)
        group = h5file[group_name]
        
        # Extract and process the data from the "x_1" dataset
        for dataset_name in keys(group)
            data = read(group[dataset_name])
            
            # Navigate through the nested dictionaries
            for (mg_key, mg_value) in data
                for (S_key, S_value) in mg_value
                    for (beta_key, beta_value) in S_value
                        # Extract x_plot and y_plot
                        x_plot = beta_value["x_plot"]
                        y_plot = beta_value["y_plot"]

                        # Ensure x_plot covers the full range from 1 to max_x
                        # full_x = 1:max_x
                        # full_y = zeros(Float64, max_x)
                        # for (i, x_val) in enumerate(x_plot)
                        #     full_y[x_val] = y_plot[i]
                        # end

                        println("Plotting for N: $group_name, mg: $mg_key, S: $S_key, beta: $beta_key")  # Debugging
                        
                        # Plotting
                        # x_plot = full_x[1:2:end]
                        y_plot = (y_plot ./ num_samples)
                        # y_plot = (y_plot ./ sum(y_plot))
                        p = plot(x_plot, y_plot; st=:bar, label="", xlabel="String Length", ylabel="Normalized Frequency", texts = round.(y_plot, digits = 6))

                        title!(p, "Distribution for $(mg_key), S=$(S_key), β=$(beta_key)")

                        # Save the plot
                        savefig("$output_dir/$(group_name)_$(mg_key)_$(S_key)_$(beta_key).png")
                    end
                end
            end
        end
    end

    close(h5file)  # Close the HDF5 file
end

# Main Execution
dag_project_number = 8
num_samples = 2000

# Collect distributions and calculate max_x
# max_x = collect_distributions(dag_project_number, num_samples)

# Recalculate max_x based on the actual distributions if needed
max_x = recalculate_max_x("/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/DAGS/$(dag_project_number)/string_length_distributions.h5")

# Plot the distributions with the fixed x-axis range
plot_distributions_with_fixed_x(max_x, dag_project_number, num_samples)
