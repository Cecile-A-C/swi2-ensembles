# Author: Cecile Coulon

import os
working_dir = './'
os.chdir(working_dir)

#%%------------------------- load paths and libraries -------------------------

import pyemu, shutil, glob
import pandas as pd
import numpy as np
import functions_model_run as fmod
import matplotlib.pyplot as plt

# Folders with inputs & outputs

pop_size = 30
base = True # There is a base realization in the parens
swp_control_file = 'pestpp_swp.pst'
swp_out_file = 'sweep_out.csv'

name_script = 'model_run_opt.py'
model_command_line = 'python3 ' + name_script
name_simu = 'simu_opt'

mou_control_file = 'pestpp_mou.pst'

climatechange = True

if climatechange == False:
    path_mou_in = os.path.join('send_to_mou_v3', 'noclimatechange')
    mou_control_file = mou_control_file.replace('.pst', '_nocc.pst')
    path_swp_in = os.path.join('..', 'swp', 'send_to_swp_v2', 'noclimatechange')
    path_swp_out = os.path.join('..', 'swp', 'results_v2', 'noclimatechange')
    swp_control_file = swp_control_file.replace('.pst', '_nocc.pst')

elif climatechange == True:
    path_mou_in = os.path.join('send_to_mou_v3', 'climatechange')
    mou_control_file = mou_control_file.replace('.pst', '_cc.pst')
    path_swp_in = os.path.join('..', 'swp', 'send_to_swp_v2', 'climatechange')
    path_swp_out = os.path.join('..', 'swp', 'results_v2', 'climatechange')
    swp_control_file = swp_control_file.replace('.pst', '_cc.pst')


#%%------------------------ Create send_to_mou folder -------------------------

if not os.path.isdir(path_mou_in):
    
    # Copy input files of the SWP run
    
    shutil.copytree(path_swp_in, path_mou_in)
    shutil.copy(os.path.join(path_swp_in, swp_control_file), 
                os.path.join(path_mou_in, mou_control_file))
    os.remove(os.path.join(path_mou_in, swp_control_file))
    os.remove(os.path.join(path_mou_in, 'model_run_steady.py'))
    shutil.copy(name_script, path_mou_in)
    
    # Copy output files from the SWP run
    
    shutil.copytree(os.path.join(path_swp_out, 'simulation_files'),
                    os.path.join(path_mou_in, 'simulation_files'))
    shutil.copy(os.path.join(path_swp_out, swp_out_file), path_mou_in)
    
    # Modify initial condition filenames from the SWP run
    
    os.mkdir(os.path.join(path_mou_in, 'init'))
    swp_out = pd.read_csv(os.path.join(path_swp_out, swp_out_file), index_col = 0)
    init_df = pd.DataFrame(index = swp_out['input_run_id'].tolist(), columns = ['hds', 'zta'])
    init_df.index.name = 'real_name'
    
    for Str in init_df.columns:
        for filepath in glob.glob(os.path.join(path_mou_in, 'simulation_files', '*.' + Str)):
            
            # Extract the SWP run_id and find the corresponding realization name
            run_id = int(filepath.split('.runid=')[1].split('.simu_swp')[0])
            real_nme = swp_out.loc[swp_out.index == run_id, 'input_run_id'].values[0]
            
            # Save to new filename & add to init_df dataframe
            filename = 'simu_swp_realname_' + real_nme + '.' + Str
            shutil.copy(filepath, os.path.join(path_mou_in, 'init', filename))
            init_df.loc[init_df.index == real_nme, Str] = filename
    
    init_file = 'init_filenames.csv'
    init_df.to_csv(os.path.join(path_mou_in, init_file), index = True)
    os.remove(os.path.join(path_mou_in, swp_out_file))
    shutil.rmtree(os.path.join(path_mou_in, 'simulation_files'))
    
    # Modify the parens filename
    
    parens_file = [os.path.basename(file) for file in glob.glob(
        os.path.join(path_mou_in, '*.csv')) if 'par' in file][0]
    # parens_file_modif = parens_file.replace('q0.csv', 'opt.csv')
    
    shutil.copy(os.path.join(path_mou_in, parens_file), 
                os.path.join(path_mou_in, parens_file.replace('q0.csv', 'opt.csv')))
    os.remove(os.path.join(path_mou_in, parens_file))


#%%--------------------- Prepare .pst for PESTPP-MOU run ----------------------

# Load .pst
pst = pyemu.Pst(os.path.join(path_mou_in, mou_control_file))
par = pst.parameter_data
obs = pst.observation_data

# Read parens file
parens_file = [os.path.basename(file) for file in glob.glob(os.path.join(path_mou_in, '*.csv')) if '.par' in file][0]
parens = pd.read_csv(os.path.join(path_mou_in, parens_file), index_col = 0)
num_reals = len(parens) # Number of realizations


#%%--------------------- Generate alpha_L and M ensembles ---------------------

if not os.path.exists(os.path.join(path_mou_in, '..', 'alphaL_M_ensembles.csv')):
    
    # Retrieve mean and stdev values from pst file 
    # (assuming par bounds represent the 95% CI i.e. the mean +- 2 stdevs)
    
    alpha_L_mean = par.loc['alpha_l', 'parval1']
    alpha_L_stdev = ( par.loc['alpha_l', 'parubnd'] - par.loc['alpha_l', 'parlbnd'] ) / 4
    M_mean = par.loc['m', 'parval1']
    M_stdev = ( par.loc['m', 'parubnd'] - par.loc['m', 'parlbnd'] ) / 4
    
    # Sample the probability distribution
    # (random sampling assuming Gaussian distribution)
    
    if base == True:
        alpha_L = np.random.normal(alpha_L_mean, alpha_L_stdev, num_reals - 1).tolist()
        alpha_L.append(alpha_L_mean) # append base
        M = np.random.normal(M_mean, M_stdev, num_reals - 1).tolist()
        M.append(M_mean) # append base
    
    else:
        alpha_L = np.random.normal(alpha_L_mean, alpha_L_stdev, num_reals).tolist()
        M = np.random.normal(M_mean, M_stdev, num_reals).tolist()
    
    # Export ensembles to csv
    
    pars = parens[['alpha_l', 'm']]
    pars.loc[:, 'alpha_l'] = alpha_L
    pars.loc[:, 'm'] = M
    pars.to_csv(os.path.join(path_mou_in, '..', 'alphaL_M_ensembles.csv'), index = True)


#%%----------------------- Add new parameter 'real_nme' -----------------------

path_model_params = 'model_params'
path_pest_files = 'pest_files'

realnme_init = 9999

if not os.path.exists(os.path.join(path_mou_in, path_pest_files, 'realnme.tpl')):
    
    real_name_df = pd.DataFrame.from_dict({'name': 'realnme', 'mean': str(realnme_init),
                                           'tpl': ["~  {0:^20}  ~".format('realnme')]})
    
    # Create .dat and .tpl files
    fmod.write_df_to_dat(os.path.join(path_mou_in, path_model_params, 'realnme'), 
                         real_name_df[['name', 'mean']])
    fmod.write_df_to_tpl(os.path.join(path_mou_in, path_pest_files, 'realnme'), 
                         real_name_df[['name', 'tpl']])
    
    
    # Add to .pst
    pst.add_parameters(os.path.join(path_mou_in, path_pest_files, 'realnme.tpl'), 
                       os.path.join(path_mou_in, path_model_params, 'realnme.dat'),
                       pst_path = path_pest_files)
    
    pst.model_input_data.pest_file[-1] = os.path.join(
        path_pest_files, 'realnme.tpl').replace('\\', '/')
    pst.model_input_data.model_file[-1] = os.path.join(
        path_model_params, 'realnme.dat').replace('\\', '/')


#%%------------------------ Add new parameter '_risk_' ------------------------

risk_init = 0.001

if not os.path.exists(os.path.join(path_mou_in, path_pest_files, 'risk.tpl')):
    
    risk_df = pd.DataFrame.from_dict({'name': '_risk_', 'mean': str(risk_init),
                                      'tpl': ["~  {0:^20}  ~".format('_risk_')]})
    
    # Create .dat and .tpl files
    fmod.write_df_to_dat(os.path.join(path_mou_in, path_model_params, 'risk'), 
                          risk_df[['name', 'mean']])
    fmod.write_df_to_tpl(os.path.join(path_mou_in, path_pest_files, 'risk'), 
                          risk_df[['name', 'tpl']])
    
    # Add to .pst
    pst.add_parameters(os.path.join(path_mou_in, path_pest_files, 'risk.tpl'), 
                        os.path.join(path_mou_in, path_model_params, 'risk.dat'),
                        pst_path = path_pest_files)
    
    pst.model_input_data.pest_file[-1] = os.path.join(
        path_pest_files, 'risk.tpl').replace('\\', '/')
    pst.model_input_data.model_file[-1] = os.path.join(
        path_model_params, 'risk.dat').replace('\\', '/')


#%%------------------------------- Modify .pst --------------------------------

par = pst.parameter_data

# Parameter transformations

p = par.parnme.apply(lambda x: x in ['alpha_l', 'm']) # Unfix alpha_L and M
par.loc[p, 'partrans'] = 'none'

if climatechange == True:
    par.loc['slvl', 'partrans'] = 'none' # Unfix slvl

# realnme parameter

par.loc['realnme', 'parval1'] = realnme_init
par.loc['realnme', 'partrans'] = 'fixed'

# zmuni_pct constraints

c = obs.loc[obs.obsnme.apply(lambda x: x.endswith('pct')), 'obsnme']
obs.loc[c, 'weight'] = 1

# Change model command line
pst.model_command[0] = model_command_line


#%%---------------------------- Decision variables ----------------------------

#---------- qmuni parameters (m3/day)

qinit = 10
p = par.parnme.apply(lambda x: x.startswith('qmuni'))
par.loc[p, 'pargp'] = 'dv_pars'
par.loc[p, 'partrans'] = 'none'
par.loc[p, 'parval1'] = qinit
par.loc[p, 'parubnd'] = 30
par.loc[p, 'parlbnd'] = 0

pst.rectify_pgroups()
pargp = pst.parameter_groups
dvar_gp = pargp.pargpnme.apply(lambda x: x in ['dv_pars'])
pargp.loc[dvar_gp, 'inctyp'] = 'absolute'
pargp.loc[dvar_gp, 'derinc'] = 1 # (m3/day)

#---------- _risk_ parameter

par.loc['_risk_', 'pargp'] = 'dv_pars'
par.loc['_risk_', 'partrans'] = 'none'
par.loc['_risk_', 'parval1'] = risk_init
par.loc['_risk_', 'parubnd'] = 0.99
par.loc['_risk_', 'parlbnd'] = 0.001


#%%---------------------- Generate initial dv population ----------------------

init_dv_file = 'dv_pop_init_modif.csv'

if not os.path.exists(os.path.join(path_mou_in, '..', init_dv_file)):

    pars_draw = [p for p in parens.columns if p.startswith('qmuni') or p == '_risk_'] # list of dvars
    pars_other = [p for p in parens.columns if p not in pars_draw]  # list of other params
    
    # Draw dvars in pars_draw from uniform distribution using info in the .pst
    
    pe_uniform = pyemu.ParameterEnsemble.from_uniform_draw(pst, num_reals = pop_size)
    dv_pop_draw = pe_uniform._df[pars_draw]
    dv_pop_draw.index.name = 'real_name'
    dv_pop_draw.index = ['gen=0_member=' + str(i) for i in dv_pop_draw.index]
    
    # Create and fill dv_pop_init dataframe
    
    dv_pop_init = pd.DataFrame(index = dv_pop_draw.index, 
                               columns = parens.columns)
    dv_pop_init.index.name = 'real_name'
    
    dv_pop_init.loc[:, pars_draw] = dv_pop_draw.loc[:, pars_draw] # uniform distributions for dvars
    for p in pars_other:
        dv_pop_init.loc[:, p] = par.loc[p, 'parval1'] # parval1 values for other parameters
    
    # Export to .csv
    dv_pop_init.to_csv(os.path.join(path_mou_in, '..', init_dv_file), index = True)

shutil.copy(os.path.join(path_mou_in, '..', init_dv_file), 
            os.path.join(path_mou_in, init_dv_file))

# Modify qmuni bounds AFTER using pyemu.ParameterEnsemble.from_uniform_draw

p = par.parnme.apply(lambda x: x.startswith('qmuni'))
par.loc[p, 'parubnd'] = 200


#%%---------------------------- Modify parens .csv ----------------------------

# Replace alpha_L and M by ensembles in csv file

alphaL_M_ensembles = pd.read_csv(os.path.join(path_mou_in, '..', 'alphaL_M_ensembles.csv'),
                                 index_col = 0)

parens.loc[:, 'alpha_l'] = alphaL_M_ensembles.loc[:, 'alpha_l']
parens.loc[:, 'm'] = alphaL_M_ensembles.loc[:, 'm']

# Add real_nme parameter
parens['realnme'] = parens.index
parens.loc[parens['realnme'] == 'base', 'realnme'] = realnme_init

# Add risk parameter
parens['_risk_'] = risk_init

# Export modified .csv
parens.to_csv(os.path.join(path_mou_in, parens_file), index = True)


#%%-------------------------- Modify PESTPP options ---------------------------

num_gens = 50
recalc_every = 10

pst.control_data.noptmax = num_gens

pst.pestpp_options = {} # Reset
pst.pestpp_options['opt_dec_var_groups'] = ['dv_pars']
pst.pestpp_options['mou_objectives'] = ['_risk_', 'pump_rate']
pst.pestpp_options['mou_risk_objective'] = True
pst.pestpp_options['mou_generator'] = 'de'
pst.pestpp_options['mou_population_size'] = pop_size
pst.pestpp_options['mou_env_selector'] = 'nsga'
pst.pestpp_options['mou_verbose_level'] = 4
pst.pestpp_options['opt_constraint_groups'] = 'l_zmuni_pct'
pst.pestpp_options['opt_chance_points'] = 'all'
pst.pestpp_options['opt_recalc_chance_every'] = recalc_every
pst.pestpp_options['opt_stack_size'] = num_reals
pst.pestpp_options['opt_par_stack'] = parens_file
pst.pestpp_options['mou_dv_population_file'] = init_dv_file


#%%---------------------------- Save modified .pst ----------------------------

pst.write(os.path.join(path_mou_in, mou_control_file))


#%----- Modify objectives

# Objective (maximize = greater_than)
pst.add_pi_equation([x for x in par.parnme if x.startswith('qmuni')], 
                    pilbl = 'pump_rate', obs_group = 'greater_than')

# Objective (maximize = greater_than)
pst.add_pi_equation(['_risk_'], pilbl = '_risk_', obs_group = 'greater_than')

#%----- Save modified .pst

pst.write(os.path.join(path_mou_in, mou_control_file))
