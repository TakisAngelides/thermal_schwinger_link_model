import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler

# Inputs needed for the finite temperature

# S = inputs["S"]
# x = inputs["x"]
# mg = inputs["mg"]
# N = inputs["N"]
# max_steps = inputs["ms"]
# save_steps_list = inputs["ssl"]
# cutoff = inputs["cutoff"]
# delta_beta = inputs["db"]

# Inputs needed for the zero temperature

# S = inputs["S"]
# x = inputs["x"]
# mg = inputs["mg"]
# N = inputs["N"]
# nsweeps = inputs["nsweeps"]
# energy_tol = inputs["energy_tol"]
# cutoff = inputs["cutoff"]
# lambd = inputs["lambd"]

# Project number of dag
project_number = os.getcwd().strip().split('/')[-1]
path_to_project_number = f'/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/DAGS/{project_number}'

def write_dag():
    
    # Name of dag
    name_of_dag = 'run.dag'
    
    # Open text file to write the dag instructions
    f_dag = open(name_of_dag, 'w')
    
    # This will contain DAGMAN_USE_DIRECT_SUBMIT = False to avoid obscure bugs of authentication
    f_dag.write(f'CONFIG /lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/dagman.config\n')
    
    # The julia files to run with the given inputs
    file_to_run_finite_temperature = 'run_finite_temperature.jl'
    file_to_run_zero_temperature = 'run_zero_temperature.jl'

    # Path to submission files to run from dag
    path_to_sub_finite_temperature = '/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/run_finite_temperature.sub'
    path_to_sub_zero_temperature = '/lustre/fs24/group/cqta/tangelides/thermal_schwinger_link_model/run_zero_temperature.sub'
    
    # Where to find the inputs
    path_to_inputs_h5 = path_to_project_number + '/inputs.h5'
    f_h5 = h5py.File(path_to_inputs_h5, 'w')
        
    # Create relevant folders if needed
    if not os.path.exists(f'{path_to_project_number}/Plots'):
        os.makedirs(f'{path_to_project_number}/Plots')
    if not os.path.exists(f'{path_to_project_number}/HDF5'):
        os.makedirs(f'{path_to_project_number}/HDF5')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Error'):
        os.makedirs(f'{path_to_project_number}/Logs/Error')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Output'):
        os.makedirs(f'{path_to_project_number}/Logs/Output')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Log'):
        os.makedirs(f'{path_to_project_number}/Logs/Log')
            
    # This will form the job id
    counter_of_jobs = 1
    
    # Finite temperature case
    max_steps_list = [500]
    delta_beta_list = [0.001]
    save_steps_lists = [[50,100,150,200,250,300,350,400,450,500]]
    for S in [1, 2, 3, 4]:
        for x in [10]:
            for mg in [0.5]:
                for N in [20]:
                    for max_steps_idx, max_steps in enumerate(max_steps_list):
                        save_steps_list = save_steps_lists[max_steps_idx]
                        delta_beta = delta_beta_list[max_steps_idx]
                        for cutoff in [1e-9, 1e-10, 1e-11, 1e-12]:
                      
                            # Memory, CPU and maximum number of days to run
                            mem, cpu, days = 32, 8, 5
                            
                            # Job id for the dag job names and path to h5 for results
                            job_id = counter_of_jobs
                            counter_of_jobs += 1 # after assigning the job_id this is incremented for the next job
                                                                        
                            # Write inputs to h5
                            g = f_h5.create_group(f'{job_id}')       
                            g.attrs["N"] = N
                            g.attrs["x"] = x
                            g.attrs["mg"] = mg
                            g.attrs["cutoff"] = cutoff
                            g.attrs["ms"] = max_steps
                            g.attrs["ssl"] = save_steps_list
                            g.attrs["db"] = delta_beta
                            g.attrs["S"] = S
                                                        
                            # Write job to dag
                            job_name = f'{job_id}'
                            f_dag.write(f'JOB ' + job_name + f' {path_to_sub_finite_temperature}\n')
                            f_dag.write(f'VARS ' + job_name + f' job_id="{job_id}" path_to_project_number="{path_to_project_number}" file_to_run="{file_to_run_finite_temperature}" cpu="{cpu}" mem="{mem}" days="{days}"\n')
                            f_dag.write('RETRY ' + job_name + ' 1\n')
    
    # Zero temperature case
    for S in [1, 2, 3, 4]:
        for x in [10]:
            for mg in [0.5]:
                for N in [20]:
                    for nsweeps in [10000]:
                        for cutoff in [1e-9, 1e-10, 1e-11, 1e-12]:
                            for energy_tol in [1e-18]:   
                                for lambd in [100]: # penalty lagrange multiplier for Gauss's law
                      
                                    # Memory, CPU and maximum number of days to run
                                    mem, cpu, days = 16, 8, 1
                                    
                                    # Job id for the dag job names and path to h5 for results
                                    job_id = counter_of_jobs
                                    counter_of_jobs += 1 # after assigning the job_id this is incremented for the next job
                                                                                
                                    # Write inputs to h5
                                    g = f_h5.create_group(f'{job_id}')       
                                    g.attrs["N"] = N
                                    g.attrs["x"] = x
                                    g.attrs["mg"] = mg
                                    g.attrs["cutoff"] = cutoff
                                    g.attrs["et"] = energy_tol
                                    g.attrs["ns"] = nsweeps
                                    g.attrs["S"] = S
                                    g.attrs["lambd"] = lambd
                                                                
                                    # Write job to dag
                                    job_name = f'{job_id}'
                                    f_dag.write(f'JOB ' + job_name + f' {path_to_sub_zero_temperature}\n')
                                    f_dag.write(f'VARS ' + job_name + f' job_id="{job_id}" path_to_project_number="{path_to_project_number}" file_to_run="{file_to_run_zero_temperature}" cpu="{cpu}" mem="{mem}" days="{days}"\n')
                                    f_dag.write('RETRY ' + job_name + ' 1\n')
        
    # Close the dag file and the h5 input file
    f_dag.close() 
    f_h5.close()
    print(f'Total number of jobs in the dag is {counter_of_jobs-1}')                                        

def plot_greens_function_vs_distance():
    
    path_to_HDF5 = f'{path_to_project_number}/HDF5'
    path_to_inputs = f'{path_to_project_number}/inputs.h5'
    f_inputs = h5py.File(path_to_inputs, 'r')
    
    counter = 0
    
    for file in os.listdir(path_to_HDF5):
        
        try:
        
            if int(file[:-3]) <= 15:
                continue
        
            g = f_inputs[file[:-3]]
            attributes_dict = {attr_name: attr_value for attr_name, attr_value in g.attrs.items()}
        
            f = h5py.File(f'{path_to_HDF5}/{file}', 'r')
            
            greens_functions_idxs_list = np.asarray(f['greens_functions_idxs']) # (n, m, n-m), distance, temperature
            greens_functions_list = np.real(np.abs(np.asarray(f['greens_functions']))) # distance, temperature
            penalty_term = np.asarray(f['penalty_term'])
            total_z = np.asarray(f['total_z'])
            
            # print(greens_functions_list)
            # print(greens_functions_idxs_list)
            
            x = greens_functions_idxs_list[2, :]
            y = greens_functions_list[:]
            zipped_lists = list(zip(x, y))
            sorted_zipped_lists = sorted(zipped_lists, key=lambda x: x[0])
            x, y = zip(*sorted_zipped_lists)
            x = x[:50]
            y = y[:50]
            plt.plot(x, y, '-x', label = f'Penalty term = {penalty_term}\nTotal Z = {total_z}')
            plt.yscale('log')
            name = '_'.join(np.ndarray.flatten(np.array([[attr_name, attr_value] for attr_name, attr_value in g.attrs.items()])))
            plt.title(name)
            plt.legend()
            plt.savefig(f'Plots/Greens_fn_vs_d/{file[:-3]}.png')
            plt.close()
                        
            # print(greens_functions_idxs_list.shape)
            # print(greens_functions_list.shape)
                                
        except:
            
            counter += 1
            print(file[:-3], attributes_dict)

# write_dag()

plot_greens_function_vs_distance()