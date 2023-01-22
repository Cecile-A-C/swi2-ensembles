# Author: Cecile Coulon

import os
working_dir = './'
os.chdir(working_dir)

#%%------------------------- load paths and libraries -------------------------

import pyemu, shutil, glob
import pandas as pd

# Folders with inputs & outputs

no_it = 2
ies_control_file = 'pestpp_ies.pst'
path_ies_in = os.path.join('..', 'ies', 'send_to_ies_v4')
path_ies_out = os.path.join('..', 'ies', 'results_v4')
path_postproc_out = os.path.join('..', 'postproc', 'results_v3')

name_script = 'model_run_steady.py'
model_command_line = 'python3 ' + name_script
name_simu = 'simu_swp_swi'

swp_control_file = 'pestpp_swp.pst'

climatechange = True

if climatechange == False:
    path_swp_in = os.path.join('send_to_swp_v2', 'noclimatechange')
    swp_control_file = swp_control_file.replace('.pst', '_nocc.pst')

elif climatechange == True:
    path_swp_in = os.path.join('send_to_swp_v2', 'climatechange')
    swp_control_file = swp_control_file.replace('.pst', '_cc.pst')


#%%------------------------ Create send_to_swp folder -------------------------

if not os.path.isdir(path_swp_in):
    
    # Copy input files of the IES run
    
    shutil.copytree(path_ies_in, path_swp_in)
    shutil.copy(os.path.join(path_ies_in, ies_control_file), 
                os.path.join(path_swp_in, swp_control_file))
    os.remove(os.path.join(path_swp_in, ies_control_file))
    os.remove(os.path.join(path_swp_in, 'prior_cov_matrix.jcb'))
    os.remove(os.path.join(path_swp_in, name_script))
    shutil.copy(os.path.join('./', name_script), 
                os.path.join(path_swp_in, swp_control_file))

    # Copy parameter ensemble file
    
    if climatechange == False:
        parens_file = ies_control_file.replace('.pst', '.{0}.par.csv'.format(no_it))
        shutil.copy(os.path.join(path_ies_out, parens_file), path_swp_in)
        
    elif climatechange == True:
        parens_file = [os.path.basename(file) for file in glob.glob(
            os.path.join(path_postproc_out, '*.csv')) if '.par' in file][0]
        shutil.copy(os.path.join(path_postproc_out, parens_file), path_swp_in)


#%%--------------------- Prepare .pst for PESTPP-SWP run ----------------------

pst = pyemu.Pst(os.path.join(path_swp_in, swp_control_file)) # Load .pst file
par = pst.parameter_data

#---------- Modify .pst

# qmuni decision variables
p = par.parnme.apply(lambda x: x.startswith('qmuni'))
par.loc[p, 'parval1'] = 0
par.loc[p, 'parlbnd'] = 0 # Mean can't be lower than lower bound

# Set qmuni to 0
parens = pd.read_csv(os.path.join(path_swp_in, parens_file))
parens.loc[:, list(par.loc[p, 'parnme'])] = 0
parens_file_modif = parens_file.replace('.csv', '.q0.csv')
parens.to_csv(os.path.join(path_swp_in, parens_file_modif), index = False)

# Change model command line
pst.model_command[0] = model_command_line


#---------- Modify PESTPP options

pst.control_data.noptmax = 0 # Number of iterations

pst.pestpp_options = {} # Reset
pst.pestpp_options['sweep_parameter_csv_file'] = parens_file_modif
pst.pestpp_options['panther_transfer_on_finish'] = [name_simu + '.hds', name_simu + '.zta']
num_reals = len(pd.read_csv(os.path.join(path_swp_in, parens_file)))
pst.pestpp_options['sweep_chunk'] = num_reals


#%%---------------------------- Save modified .pst 

pst.write(os.path.join(path_swp_in, swp_control_file))

os.remove(os.path.join(path_swp_in, parens_file))
