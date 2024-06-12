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

def write_dag():
    
    # Name of dag
    name_of_dag = 'run.dag'
    
    # Open text file to write the dag instructions
    f_dag = open(name_of_dag, 'w')
    
    # This will contain DAGMAN_USE_DIRECT_SUBMIT = False to avoid obscure bugs of authentication
    f_dag.write(f'CONFIG /lustre/fs24/group/cqta/tangelides/OQS/dagman.config\n')
    
    # The julia file to run with the given inputs
    file_to_run = 'run.jl'
    
    # Path to submission files to run from dag
    path_to_sub_finite_temperature = '/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/run_finite_temperature.sub'
    path_to_sub_zero_temperature = '/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/run_zero_temperature.sub'

    # Project number of dag
    project_number = os.getcwd().strip().split('/')[-1]
    
    # Where to find the inputs
    path_to_project_number = f'/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/{project_number}'
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
    max_steps_list = [3]
    delta_beta_list = [0.01]
    save_steps_lists = [[1,2,3]]
    for S in [1]:
        for x in [1]:
            for mg in [1]:
                for N in [4]:
                    for max_steps_idx, max_steps in enumerate(max_steps_list):
                        save_steps_list = save_steps_lists[max_steps_idx]
                        delta_beta = delta_beta_list[max_steps_idx]
                        for cutoff in [1e-6]:
                      
                            # Memory, CPU and maximum number of days to run
                            mem, cpu, days = 1, 1, 1
                            
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
                            f_dag.write(f'VARS ' + job_name + f' job_id="{job_id}" path_to_project_number="{path_to_project_number}" file_to_run="{file_to_run}" cpu="{cpu}" mem="{mem}" days="{days}"\n')
                            f_dag.write('RETRY ' + job_name + ' 1\n')
    
    # Zero temperature case
    for S in [1]:
        for x in [1]:
            for mg in [1]:
                for N in [4]:
                    for nsweeps in [10000]:
                        for cutoff in [1e-6]:
                            for energy_tol in [1e-18]:   
                                for lambd in [100]: # penalty lagrange multiplier for Gauss's law
                      
                                    # Memory, CPU and maximum number of days to run
                                    mem, cpu, days = 1, 1, 1
                                    
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
                                    f_dag.write(f'VARS ' + job_name + f' job_id="{job_id}" path_to_project_number="{path_to_project_number}" file_to_run="{file_to_run}" cpu="{cpu}" mem="{mem}" days="{days}"\n')
                                    f_dag.write('RETRY ' + job_name + ' 1\n')
        
    # Close the dag file and the h5 input file
    f_dag.close() 
    f_h5.close()
    print(f'Total number of jobs in the dag is {counter_of_jobs-1}')                                        

write_dag()
