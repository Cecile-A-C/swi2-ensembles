# Author: Cecile Coulon

import os
working_dir = './'
os.chdir(working_dir)
modflow_path = 'mf2005'

#%%------------------------- load paths and libraries -------------------------

import flopy, pyemu, pickle, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functions_model_run as fmod

# Folders with inputs & outputs

path_input_setup_files = 'input_setup_files'
mnw2_data_file = 'mnw2_data.csv'
path_bin_in = 'preproc_IDM_20m_v6.bin'
path_pp_zones_bin_in = 'pp_zones_v1.bin'

path_output_setup_files = 'send_to_ies_v4'
path_model_outputs = 'model_outputs'
path_model_params = 'model_params'
path_pest_files = 'pest_files'

model_run_script = 'model_run_steady.py'
model_command_line = 'python3 ' + model_run_script

pst_control_file = 'pestpp_ies.pst'
name_swi = 'simu_ies_swi'


#%%------------------------ Create send_to_ies folder -------------------------

if not os.path.isdir(path_output_setup_files):
    os.mkdir(path_output_setup_files)
    os.mkdir(os.path.join(path_output_setup_files, path_model_outputs))
    os.mkdir(os.path.join(path_output_setup_files, path_model_params))
    os.mkdir(os.path.join(path_output_setup_files, path_pest_files))
    shutil.copy(os.path.join(path_input_setup_files, path_model_params, mnw2_data_file), path_output_setup_files)
    shutil.copy(os.path.join(path_input_setup_files, path_bin_in), path_output_setup_files)
    shutil.copy(os.path.join('./', model_run_script), path_output_setup_files)
    shutil.copy(os.path.join('./', 'functions_model_run.py'), path_output_setup_files)


#%%----------------------- Load GIS data from bin files -----------------------

with open(os.path.join(path_input_setup_files, path_bin_in), 'rb') as handle:
     objects = pickle.load(handle)
     (nrow, ncol, delr, delc, ibound, sea, dem, geology, thk_valley, cx, cy, 
      muni_wells_row_col, muni_pumping, domestic_wells_row_col, ind_old_wells_row_col,
      obs_wells_row_col, obs_head, obs_zeta, tdem_row_col, ert_row_col) = objects

with open(os.path.join(path_input_setup_files, path_pp_zones_bin_in), 'rb') as handle: # Pilot point zones
     objects_pp = pickle.load(handle)
     (pp_zones) = objects_pp

pp_zones[ (pp_zones == 0) & (ibound == 1) ] = 2 # Zone in which pps are tied (pp zone no.2)


#%%----------------------------- Load parameters ------------------------------

# Prior parameter set (prior values + prior param uncertainty)
# (upper & lower bounds defined by mean +-2 stdev)

param_df = pd.read_csv(
        os.path.join(path_input_setup_files, 'param_prior.csv'), 
        header = 0, names = ['name', 'mean', 'lbnd', 'ubnd'])

# Decision variables
qmuni = pd.read_csv(os.path.join(path_input_setup_files, 'qmuni.csv'))


#%%------------------------- Export .dat & .tpl files -------------------------

# Parameters (except pilot points)

param_df['tpl'] = ["~  {0:^20}  ~".format(parname) for parname in param_df.name]

p = param_df.name.apply(lambda x: x not in ['hk_EDC'])

fmod.write_df_to_dat(os.path.join(path_output_setup_files, 
                                  path_model_params, 'param'), 
                     param_df.loc[p, ['name', 'mean']])              # Write .dat

fmod.write_df_to_tpl(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'param'), 
                     param_df.loc[p, ['name', 'tpl']])               # Write .tpl

param_df.index = param_df.name.str.lower()


# Decision variables

qmuni['tpl'] = ["~  {0:^20}  ~".format(dvarnme) for dvarnme in qmuni.name]

fmod.write_df_to_dat(os.path.join(path_output_setup_files, 
                                  path_model_params, 'qmuni'), 
                     qmuni[['name', 'calibration_m3_s']])                   # Write .dat

fmod.write_df_to_tpl(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'qmuni'), 
                     qmuni[['name', 'tpl']])                    # Write .tpl

qmuni.index = qmuni.name.str.lower()


#%%------------------ Export .dat & .tpl files: pilot points ------------------

p = param_df.name.apply(lambda x: x in ['hk_EDC'])

#---- Load model and its spatial reference

mf = flopy.modflow.Modflow.load(
        os.path.join('./', name_swi + '.nam'), 
        exe_name = modflow_path,  verbose = False, forgive = False) # Load model


xmin, xmax, ymin, ymax = cx[0,0] - delr/2, cx[0,-1] + delr/2, cy[-1,0] - delc/2, cy[0,0] + delc/2

sr = pyemu.helpers.SpatialReference().from_namfile(
        os.path.join(mf.model_ws, mf.namefile),
        delr = mf.modelgrid.delr, delc = mf.modelgrid.delc)
sr.xul = xmin
sr.yul = ymax

#---- Set up pilot point locations (regularly-spaced pp grid), and export to .tpl, .dat, .shp files

prefix_dict= {0: ['hk1']} # dict of pp prefixes {layer: prefix}
every_n_cell = 25 # pp spacing

pp_df = pyemu.pp_utils.setup_pilotpoints_grid(
        ml = None, sr = sr, ibound = pp_zones, use_ibound_zones = True, 
        prefix_dict = prefix_dict, every_n_cell = every_n_cell, 
        pp_dir = os.path.join(path_output_setup_files, path_model_params), 
        tpl_dir = os.path.join(path_output_setup_files, path_pest_files), 
        shapename = os.path.join('./', 'pp.shp')) # dataframe summarizing pp info

#---- Change the default parnames

pp_df['parnme2'] = 'hk1' + pp_df.name.str.split('pp_0', 0, expand = True)[1]
pp_df.tpl = pp_df.parnme2.apply(lambda x: "~    {0}    ~".format(x))

#---- Replace the default information in pp_df with real values

pp_df['parval1'] = param_df.loc[p, 'mean'].values[0]
pp_df['parlbnd'] = param_df.loc[p, 'lbnd'].values[0]#1e-8
pp_df['parubnd'] = param_df.loc[p, 'ubnd'].values[0]#1e-2

#---- Tie the pps where needed

# Define the reference pp for the tied zone
pp_ref_tied = pp_df.loc[pp_df.zone == 2, 'parnme2'].values[0] # 'hk1000' or 'hk1_i:12_j:262_zone:2.0'

# Tie pilot points where pp_df.zone == 2
pp_df.loc[(pp_df.zone == 2) & (pp_df.parnme2 != pp_ref_tied), 'partrans'] = 'tied'
pp_df.loc[(pp_df.zone == 2) & (pp_df.parnme2 != pp_ref_tied), 'partied'] = pp_ref_tied

#---- Update the .dat and .tpl files with the real parval1 value now in pp_df
#(was created with default parval1 value = 1.0)

pp_df['pp_zones'] = pp_df.zone
pp_df.zone = 1

pyemu.utils.pp_utils.write_pp_file(
        os.path.join(path_output_setup_files, path_model_params, 'hk1pp.dat'), pp_df)

pp_df.index = pp_df.name # Change just for .tpl export

pyemu.utils.helpers._write_df_tpl(
    os.path.join(path_output_setup_files, path_pest_files, 'hk1pp.dat.tpl'), 
    pp_df.loc[:, ["name", "x", "y", "zone", "tpl"]], sep=" ", 
    index_label="index", header=False, index=False, quotechar=" ", quoting=2,)

pp_df.index = pp_df.parnme2

#---- Calculate pp variogram

# Create exponential variogram with sill (default) and range
Range = every_n_cell * mf.dis.delr.array[0] * 3.0
variogram = pyemu.geostats.ExpVario(contribution = 1.0, a = Range)

# Create gs object mimicking the behavior of a PEST geostatistical structure
gs = pyemu.geostats.GeoStruct(variograms = variogram, transform = 'log')

#---- Prepare calculation of the prior parameter covariance matrix

gs.to_struct_file(os.path.join('./', 'struct.dat')) # Write gs info to a PEST-style structure file

sd = {'struct.dat': [os.path.join(path_output_setup_files, 
                                  path_pest_files, 'hk1pp.dat.tpl')]} # {structure file: [pp template files]}

#---- DO ONLY ONCE PER GRID:

hk1pp_df = pyemu.pp_utils.pp_file_to_dataframe(
        os.path.join(path_output_setup_files, 
                    path_model_params, 'hk1pp.dat')) # Read the updated .dat file

# Ordinary kriging
ok = pyemu.geostats.OrdinaryKrige(geostruct = gs, point_data = hk1pp_df)

# Calculate kriging factors (structured grid)
ok.calc_factors_grid(sr, zone_array = ibound, minpts_interp = 1,
                    maxpts_interp = 30, search_radius = 1.0e10, 
                    num_threads = 4)

# Write grid-based PEST-style factors file
ok.to_grid_factors_file(os.path.join(path_output_setup_files, 
                                    path_model_params, 'pp_30max.fac'))


#%%---------------------------- Load observations -----------------------------

#---- Load observation files

hwells_obs = pd.read_csv(os.path.join(path_input_setup_files, 'hwells_obs.csv')) # Heads
ztdem_obs = pd.read_csv(os.path.join(path_input_setup_files, 'ztdem_obs.csv')) # TDEM interface obs
zert_obs = pd.read_csv(os.path.join(path_input_setup_files, 'zert_obs.csv')) # ERT interface obs
zwells_all = pd.read_csv(os.path.join(path_input_setup_files, 'zwells_obs.csv')) # Well interface obs + forecasts

zwells_obs = zwells_all.copy()[zwells_all.name.str.startswith('zpz')] # Subgroup: interface observations
zmuni_fcst = zwells_all.copy()[zwells_all.name.str.startswith('zmuni')] # Subgroup: zeta forecasts
zmuni_pct_fcst = zmuni_fcst.copy() # Subgroup: 1% seawater salinity contour forecasts
zmuni_pct_fcst.name = zmuni_pct_fcst.name.astype(str) + '_pct'

vlens_fcst = pd.DataFrame([('vlens', 0)], columns = ['name', 'value']) # Volume of the FW lens


#---- Drop some observations

# Inconsistent TDEM observations
ztdem_obs.drop(ztdem_obs[ (ztdem_obs['name'] == 'ztdem07') | 
        (ztdem_obs['name'] == 'ztdem23') | (ztdem_obs['name'] == 'ztdem67') | 
        (ztdem_obs['name'] == 'ztdem82') ].index, inplace=True)

# ERT observations close (< 100 m) to pumping wells
#zert_obs.drop(zert_obs[ (zert_obs['name'] == 'zert08') | 
#        (zert_obs['name'] == 'zert62') | (zert_obs['name'] == 'zert71') | 
#        (zert_obs['name'] == 'zert72') | (zert_obs['name'] == 'zert73') | 
#        (zert_obs['name'] == 'zert77') | (zert_obs['name'] == 'zert78') | 
#        (zert_obs['name'] == 'zert79') ].index, inplace=True)


#---- Compute weights

# Standard deviation of geophysical observations
sigma_tdem = 15/4 # TDEM (m)
sigma_ert = 20/4 # ERT (m)

# Model observations: weights = inverse of standard deviation
hwells_obs['weight'] = 1./hwells_obs.sigma_tot
zwells_obs['weight'] = 1/zwells_obs.sigma_tot
ztdem_obs['weight'] = 1./sigma_tdem
zert_obs['weight'] = 1./sigma_ert

# Model forecasts
zmuni_fcst['weight'] = 0
zmuni_pct_fcst['weight'] = 0
vlens_fcst['weight'] = 0

#---- Concatenate all observations in a big dataframe to facilitate filling out the .pst

cols = ['name', 'value', 'weight']

all_obs = pd.concat([ hwells_obs[cols], zwells_obs[cols], 
                     ztdem_obs[cols], zert_obs[cols], 
                     zmuni_fcst[cols], zmuni_pct_fcst[cols], 
                     vlens_fcst[cols] ])

all_obs.index = all_obs.name


#%%---------------------------- Export .ins files -----------------------------

# Define width of values (number of characters)
VAL_START = 20
VAL_CHAR_LEN = 45

#---- Create the column containing the instruction line

hwells_obs['ins'] = ['l1 ({0}){1}:{2}'.format(name,VAL_START,VAL_CHAR_LEN) for name in hwells_obs.name]
zwells_obs['ins'] = ['l1 ({0}){1}:{2}'.format(name,VAL_START,VAL_CHAR_LEN) for name in zwells_obs.name]
ztdem_obs['ins'] = ['l1 ({0}){1}:{2}'.format(name,VAL_START,VAL_CHAR_LEN) for name in ztdem_obs.name]
zert_obs['ins'] = ['l1 ({0}){1}:{2}'.format(name,VAL_START,VAL_CHAR_LEN) for name in zert_obs.name]
zmuni_fcst['ins'] = ['l1 ({0}){1}:{2}'.format(name,VAL_START,VAL_CHAR_LEN) for name in zmuni_fcst.name]
zmuni_pct_fcst['ins'] = ['l1 ({0}){1}:{2}'.format(name,VAL_START,VAL_CHAR_LEN) for name in zmuni_pct_fcst.name]
vlens_fcst['ins'] = ['l1 ({0}){1}:{2}'.format(name,VAL_START,VAL_CHAR_LEN) for name in vlens_fcst.name]


#---- Export to .ins files

fmod.write_df_to_ins(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'hwells'), 
                     hwells_obs[['ins']])

fmod.write_df_to_ins(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'zwells'), 
                     zwells_obs[['ins']])

fmod.write_df_to_ins(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'ztdem'), 
                     ztdem_obs[['ins']])

fmod.write_df_to_ins(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'zert'), 
                     zert_obs[['ins']])

fmod.write_df_to_ins(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'zmuni'), 
                     zmuni_fcst[['ins']])

fmod.write_df_to_ins(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'zmuni_pct'), 
                     zmuni_pct_fcst[['ins']])

fmod.write_df_to_ins(os.path.join(path_output_setup_files, 
                                  path_pest_files, 'vlens'), 
                     vlens_fcst[['ins']])


#%%---------------- Create an initial PEST control file (.pst) ----------------

#---- Path to PEST files

in_files = [os.path.join(path_output_setup_files, path_model_params, filename) \
            for filename in os.listdir(os.path.join(path_output_setup_files, path_model_params)) \
            if filename.endswith('.dat')] # Param .dat files

tpl_files = [os.path.join(path_output_setup_files, path_pest_files, filename) \
             for filename in os.listdir(os.path.join(path_output_setup_files, path_pest_files)) \
             if filename.endswith('.tpl')] # Param .tpl files

out_files = [os.path.join(path_output_setup_files, path_model_outputs, filename) \
             for filename in os.listdir(os.path.join(path_output_setup_files, path_model_outputs)) \
             if filename.endswith('.dat')] # Output .dat files

ins_files = [os.path.join(path_output_setup_files, path_pest_files, filename) \
             for filename in os.listdir(os.path.join(path_output_setup_files, path_pest_files)) \
             if filename.endswith('.ins')] # Output .ins files

for List in [tpl_files, in_files, ins_files, out_files]:
    List.sort() # Sort lists so that .tpl and .ins files are assigned to the correct .dat files


#---- Initialize pst control file

pst = pyemu.Pst.from_io_files(tpl_files, in_files, ins_files, out_files) # .pst filled with generic values

pst.model_input_data.model_file = [os.path.join(*(p.split(os.path.sep)[1:])).replace('\\', '/') for p in pst.input_files]
pst.model_input_data.pest_file = [os.path.join(*(p.split(os.path.sep)[1:])).replace('\\', '/') for p in pst.template_files]
pst.model_output_data.model_file = [os.path.join(*(p.split(os.path.sep)[1:])).replace('\\', '/') for p in pst.output_files]
pst.model_output_data.pest_file = [os.path.join(*(p.split(os.path.sep)[1:])).replace('\\', '/') for p in pst.instruction_files]


#%%----------------------- Modify the .pst: parameters ------------------------

par = pst.parameter_data

#---- Method for derivatives calculation

pst.parameter_groups.forcen = 'always_3'

#---- Fixed parameters

p = par.parnme.apply(lambda x: x in ['slvl', 'ss', 'sy', 'ne', 'alpha_l', 'm'])
par.loc[p, 'partrans'] = 'fixed'

#---- Adjustable parameters (except pilot points)

p = par.parnme.apply(lambda x: x in ['rech_mmyr'])
par.loc[p, 'partrans'] = 'none'

p = par.parnme.apply(lambda x: x in param_df.index)
par.loc[p, 'parval1'] = param_df.loc[par.loc[p, 'parnme'], 'mean']

p = par.parnme.apply(lambda x: x in param_df.index and x != 'slvl')
par.loc[p,'parubnd'] = param_df.loc[par.loc[p, 'parnme'], 'ubnd']
par.loc[p,'parlbnd'] = param_df.loc[par.loc[p, 'parnme'], 'lbnd']

p = par.parnme.apply(lambda x: x in ['rech_mmyr']) #?
par.loc[p, 'partrans'] = 'none'

#---- Pilot point parameters

p = par.parnme.apply(lambda x: x.startswith('hk1'))
par.loc[p, 'parval1'] = pp_df.loc[p, 'parval1']
par.loc[p,'parubnd'] = pp_df.loc[p, 'parubnd']
par.loc[p,'parlbnd'] = pp_df.loc[p, 'parlbnd']
par.loc[p,'partrans'] = pp_df.loc[p, 'partrans']
par.loc[p,'partied'] = pp_df.loc[p, 'partied']
par.loc[p, 'pargp'] = 'hk_pp'

pst.rectify_pgroups()


#%%------------------- Modify the .pst: decision variables --------------------

p = par.parnme.apply(lambda x: x.startswith('qmuni')) # qmuni decision variables

par.loc[p, 'parval1'] = qmuni.loc[p, 'calibration_m3_s']
par.loc[p, 'partrans'] = 'fixed' # fixed for IES
par.loc[p, 'pargp'] = 'qmuni'

pst.rectify_pgroups()

#---- Change perturbation increment

pargp = pst.parameter_groups

dvar_gp = pargp.pargpnme.apply(lambda x: x in ['qmuni'])
pargp.loc[dvar_gp, 'inctyp'] = 'absolute'
pargp.loc[dvar_gp, 'derinc'] = 999

#---- Compute the prior parameter covariance matrix

cov = pyemu.helpers.geostatistical_prior_builder(
        pst, struct_dict = sd, sigma_range = 4)

cov.to_binary(os.path.join(path_output_setup_files, 
                           'prior_cov_matrix.jcb'))


#%%----------------------- Modify the .pst: model outputs ------------------------

obs = pst.observation_data

#---- Define observation & constraint groups

obgnmes = ['hmuni', 'hdeep', 'hshallow', 'zwells', 'ztdem', 'zert', 
           'zmuni', 'zmuni_pct', 'v_lens'] # Obs group names

obgconditions = [ obs.obsnme.apply(lambda x:x.startswith('hmuni')),
          obs.obsnme.isin(hwells_obs.loc[hwells_obs.well_type == 'hdeep', 'name'].tolist()),
          obs.obsnme.isin(hwells_obs.loc[hwells_obs.well_type == 'hshallow', 'name'].tolist()),
          obs.obsnme.apply(lambda x:x.startswith('zpz')),
          obs.obsnme.apply(lambda x:x.startswith('ztdem')),
          obs.obsnme.apply(lambda x:x.startswith('zert')),
          obs.obsnme.apply(lambda x:x.startswith('zmuni') and not x.endswith('pct')),
          obs.obsnme.apply(lambda x:x.endswith('pct')), 
          obs.obsnme.apply(lambda x:x.startswith('vlens')) ]

dict_obg = dict(zip(obgnmes, obgconditions))

for obgname, obgcondition in dict_obg.items():
    groups = obs.groupby(obgcondition).groups
    obs.loc[groups[True], 'obgnme'] = obgname


#---- zmuni_pct constraints

c = obs.loc[obs.obsnme.apply(lambda x: x.endswith('pct')), 'obsnme']

# Load right-hand side of the constraint
muni_data = pd.read_csv(os.path.join(path_input_setup_files, mnw2_data_file))
muni_data.index = 'z' + muni_data.wellid.astype(str) + '_pct'

obs.loc[c, 'obgnme'] = 'l_' + obs.loc[c, 'obgnme'] # 'less than' constraint
obs.loc[c, 'obsval'] = muni_data.loc[c, 'zbotm'] # Set RHS of the constraint

obs.loc[c, 'weight'] = 0

#---- Observations

notc = obs.obsnme.apply(lambda x: x not in c)
obs.loc[notc, 'obsval'] = all_obs.loc[notc, 'value']
obs.loc[notc, 'weight'] = all_obs.loc[notc, 'weight']

                                                                
#%%--------------------------- Final .pst settings ----------------------------

pst.control_data.noptmax = 2 # Number of iterations
pst.model_command[0] = model_command_line

#---- IES options

pst.pestpp_options = {}
pst.pestpp_options['ies_num_reals'] = 200
pst.pestpp_options['parcov'] = 'prior_cov_matrix.jcb'
pst.pestpp_options['ies_add_base'] = True
pst.pestpp_options['ies_enforce_bounds'] = True
pst.pestpp_options['enforce_tied_bounds'] = True
pst.pestpp_options['ies_use_approx'] = False
pst.pestpp_options['ies_save_rescov'] = True
pst.pestpp_options['ies_no_noise'] = False
pst.pestpp_options['ies_drop_conflicts'] = True

#---------- Write pst to PEST control file

pst.write(os.path.join(path_output_setup_files, pst_control_file))