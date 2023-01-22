# Author: Cecile Coulon

import os
working_dir = './'
os.chdir(working_dir)
modflow_path = 'mf2005'


#%%------------------------- load paths and libraries -------------------------

# Libraries and function
import pickle, copy, pyemu, flopy
import pandas as pd
import numpy as np
import functions_model_run as fmod

# Input & output files
path_model_outputs = 'model_outputs' # Folder with model outputs
path_model_params = 'model_params' # Folder with model inputs
mnw2_data_file = 'mnw2_data.csv' # MNW2 information
krig_factors_file = 'pp_30max.fac' # Kriging factors
hk_pp_file = 'hk1pp.dat' # hk values at pilot points
param_file = 'param.dat' # Other parameter values
qmuni_file = 'qmuni.dat' # qmuni values (= decision variables)

# Bin files generated with QGridder
path_bin_in = 'preproc_IDM_20m_v6.bin' #path_pp_zones_bin_in = 'pp_zones_v1.bin'

name_simu = 'simu_ies' # Model name
model_ws = './' # Path to save simulation files


#%%----------------------- Load GIS data from bin file ------------------------

# Retrieve data from the bin file generated with QGridder
with open(path_bin_in, 'rb') as handle:
     objects = pickle.load(handle)
     (nrow, ncol, delr, delc, ibound, sea, dem, geology, thk_valley, cx, cy, 
      muni_wells_row_col, muni_pumping, domestic_wells_row_col, ind_old_wells_row_col,
      obs_wells_row_col, obs_head, obs_zeta, tdem_row_col, ert_row_col) = objects

# Create row_col dictionaries for wells with head and zeta observations or forecasts
Dict_z = {key: value for key, value in obs_zeta.items() if value == 1}
obs_z_row_col = {key: value for key, value in obs_wells_row_col.items() if key in Dict_z.keys()}
obs_z_row_col.update(muni_wells_row_col)

Dict_h = {key: value for key, value in obs_head.items() if value == 1}
obs_h_row_col = {key: value for key, value in obs_wells_row_col.items() if key in Dict_h.keys()}
obs_h_row_col.update(muni_wells_row_col)

# Spatial references
xmin, xmax, ymin, ymax = cx[0,0] - delr/2, cx[0,-1] + delr/2, cy[-1,0] - delc/2, cy[0,0] + delc/2

#%%----------------------------- Model parameters -----------------------------

# Read parameters
param_df = pd.read_csv(os.path.join(path_model_params, param_file), 
                       index_col = 0, delim_whitespace = True, 
                       header = None, names = ['name', 'value'])

hk_dunes = param_df.loc['hk_dunes', 'value']
hk_valley = param_df.loc['hk_valley', 'value']
hk_seabed = param_df.loc['hk_seabed', 'value']
rech_mmyr = param_df.loc['rech_mmyr', 'value']
alpha_T = param_df.loc['alpha_T', 'value']
sy = param_df.loc['sy', 'value']
ss = param_df.loc['ss', 'value']
ne = param_df.loc['ne', 'value']
slvl = param_df.loc['slvl', 'value']
alpha_L = param_df.loc['alpha_L', 'value']
M = param_df.loc['M', 'value']

# Create hk grid from hk values determined at pilot points ('hk1pp.dat' pp_file) & kriging factors ('pp.fac' factors_file)
hk_EDC_array = pyemu.utils.fac2real(pp_file = os.path.join(path_model_params, hk_pp_file), 
                                    factors_file = os.path.join(path_model_params, krig_factors_file), 
                                    out_file = None, fill_value = np.nan)

# Decision variables: Q at municipal wells
qmuni_df = pd.read_csv(os.path.join(path_model_params, qmuni_file), 
                       index_col = 0, delim_whitespace = True, 
                       header = None, names = ['name', 'value'])


#%%================================ 1) No SWI =================================

# Discretization package
nlay = 1 # Number of layers (along height)
top = 100 # Top elevation of layer
botm = -300 # Bottom elevation
itmuni = 1 # Time units (seconds)
lenuni = 2 # Length units (meters)

nper = 1 # Number of stress periods
perlen = 2*24*3600 # Length of stress periods (2 days, converted to seconds)
nstp = 1 # Number of time steps in each stress period
steady = True # Type of simulation (steady or transient)
tsmult = 1 # Time step multiplier

# Basic package
h0 = 0 # Initial heads (m)

# Layer-Property Flow package
laytyp = 1 # Convertible layer (T = f(head) throughout simulation)
thk_dunes = 10 # m Thickness of sand dunes

# General-Head Boundary package
seabed_thk = 150 # Seabed thickness (m) = aquifer thickness/2
rho_f = 1000 # Freshwater density (kg/m3)
rho_s = 1025 # Saltwater density (kg/m3)

# Well package WEL
domestic_use_m3d = 80 # Total amount of pumped water from individual wells (m3/day)

# Revised Multi-Node Well package MNW2
nnodes = 0 # Number of cells associated with the well
losstype = 'thiem' # Model for well loss
ppflag = 0 # Partial penetration flag


#%%--------------------------- Parameter processing ---------------------------

# Basic package
h_init = h0 * np.ones((nlay, nrow, ncol)) # Grid of initial heads

# Layer-Property Flow package
hk_array = copy.deepcopy(hk_EDC_array) # Initialize
hk_valleyH = ( hk_EDC_array * (dem-botm - thk_valley) + hk_valley * thk_valley ) / (dem-botm)
hk_dunesH = ( hk_EDC_array * (dem-botm - thk_dunes) + hk_dunes * thk_dunes ) / (dem-botm)
hk_array[geology == 3] = hk_valleyH[geology == 3] # Paleoglacial valley
hk_array[geology == 4] = hk_dunesH[geology == 4] # Sand dunes

# Recharge package
recharge_ms = rech_mmyr / (1000*365.25*24*3600) # Convert recharge values from mm/year to m/s
rech = np.zeros((nrow, ncol)) # Recharge grid
rech[ sea == 0 ] = recharge_ms # For onshore cells, assign recharge

# General-Head Boundary package
epsilon = (rho_s - rho_f) / rho_f # Density ratio
corr_factor = 1 - ( alpha_T / (dem.mean() - botm) ) ** (1/4) # Density ratio correction factor (Lu & Werner 2013)
epsilon_corr = epsilon * corr_factor # Corrected density ratio

ghbc_EDC = (hk_seabed / seabed_thk) * delr * delc # Boundary conductance for EDC geology (m2/s)
ghb_cond = np.ones((nrow, ncol)) * ghbc_EDC # Default geology = EDC
ghb_stage = (epsilon + 1) * slvl + (-epsilon) * dem # Stage (FW head at the ocean bottom) (m)
ghb_stage[sea != 1] = np.nan # GHB stage = nan for onshore cells

idx = np.logical_and(sea == 1, ibound == 1) # Boolean (indexes of active & ocean cells)
lay = 0 # Layer in which to implement GHB

ghb_sp_data = []
for row_col in np.argwhere(idx): # For indices where idx is True
    row, col = row_col
    ghb_sp_data.append( [lay, row , col, ghb_stage[row,col], ghb_cond[row,col]] )
    # Each ghb cell defined by: [lay, row, col, stage, conductance]

# Well package WEL
domestic_use_m3s = domestic_use_m3d / (24 * 3600) # Total domestic use: conversion from m3/day to m3/s
domestic_use_m3s_pc = domestic_use_m3s / len(domestic_wells_row_col) # Use per well = total use/number of wells

wel_sp_data = {}
lay = 0 # Layer in which wells are pumping

for i in range(nper): # For each stress period
    
    well_data = []
    for well in domestic_wells_row_col.keys():
        row = domestic_wells_row_col[well][0][0]
        col = domestic_wells_row_col[well][0][1]
        well_data.append([lay, row , col, -domestic_use_m3s_pc])
        
    wel_sp_data[i] = well_data
    # WEL dictionary = {nper: [lay, row, col, pumping rate]}

# Revised Multi-Node Well package MNW2
node_data_df = pd.read_csv(os.path.join(path_model_params, mnw2_data_file)) # Table with MNW2 data by node
node_data_df['nnodes'] = nnodes
node_data_df['ppflag'] = ppflag
node_data_df['losstype'] = losstype
node_data = node_data_df.to_records() # Convert dataframe to rec array

# Stress period information for MNW2
per = [0] # List of stress periods
active_mnw2_wells = len(node_data_df) # Number of active MNW2 wells
wel_sp_data_mnw2_df = pd.DataFrame(list(zip(per * active_mnw2_wells, 
                                            node_data_df.wellid, 
                                            -qmuni_df.value)), 
               columns =['per', 'wellid', 'qdes']) # qdes = actual volumetric pumping rate at the well (m3/s, <0 for withdrawal)

pers = wel_sp_data_mnw2_df.groupby('per') # Group wel_sp_data_mnw2_df by periods
wel_sp_data_mnw2 = {i: pers.get_group(i).to_records() for i in range(nper)} # Convert df to dictionary

# Multi-Node Well Information Package MNWI
unit_mnwi, qndflag, qbhflag = 35, 0, 0 # Unit number for output files, flags for writing additional flow info in the output files
mnwi_list = [list(x) for x in zip(node_data_df.wellid, 
             [unit_mnwi] * active_mnw2_wells, 
             [qndflag] * active_mnw2_wells, 
             [qbhflag] * active_mnw2_wells)]

# Output Control package
oc_list = ['save head']#, 'save budget']
spd = {(i,0) : oc_list for i in range(nper)} # Save head & budget for (stress period i, timestep 1)


#%%------------------------------- Build model --------------------------------

name_noswi = name_simu + '_noswi'

# Create modflow model object
mf = flopy.modflow.Modflow(modelname = name_noswi, 
                           model_ws = model_ws, 
                           exe_name = modflow_path)


# Discretization package
dis = flopy.modflow.ModflowDis(mf, nlay = nlay, nrow = nrow, ncol = ncol, 
                               nper = nper, delr = delr, delc = delc, 
                               top = top, botm = botm, 
                               perlen = perlen, nstp = nstp, 
                               tsmult = tsmult, steady = steady, 
                               itmuni = itmuni, lenuni = lenuni)

# Set spatial reference info (after defining dis package)
mf.modelgrid.set_coord_info(xoff = xmin, yoff = ymin, angrot = 0, epsg = 2946)

# Basic package
bas = flopy.modflow.ModflowBas(mf, ibound = ibound, strt = h_init)
# Recharge package
rch = flopy.modflow.ModflowRch(mf, rech = rech)
# Layer-Property Flow package
lpf = flopy.modflow.ModflowLpf(mf, laytyp = laytyp, hk = hk_array, 
                               sy = sy, ss = ss)
# General-Head Boundary package
ghb = flopy.modflow.ModflowGhb(mf, stress_period_data = ghb_sp_data)
# Well package
wel = flopy.modflow.ModflowWel(mf, stress_period_data = wel_sp_data)
# Revised Multi-Node Well package
mnw2 = flopy.modflow.ModflowMnw2(mf, mnwmax = active_mnw2_wells, 
                                 node_data = node_data,
                                 stress_period_data = wel_sp_data_mnw2,
                                 itmp = [active_mnw2_wells] * nper)
# Multi-Node Well Information Package
mnwi = flopy.modflow.ModflowMnwi(mf, byndflag = 36, 
                                 mnwobs = active_mnw2_wells, # bynd = output file for MNW2 information
                                 wellid_unit_qndflag_qhbflag_concflag = mnwi_list)
# Output Control package
oc = flopy.modflow.ModflowOc(mf, stress_period_data = spd)
# Preconditioned Conjugate-Gradient package
pcg = flopy.modflow.ModflowPcg(mf)

# Write input files
mf.write_input()

#%%-------------------------------- Run model ---------------------------------

mf.run_model(silent=False)

#%%------------------------------- Read outputs -------------------------------

hds_noswi = flopy.utils.binaryfile.HeadFile(name_noswi + '.hds')
h_noswi = hds_noswi.get_alldata()[0,0,:,:] # Retrieve all data
h_noswi[ h_noswi < -888 ] = 0 # Assign nan to no data values


#%%================================== 2) SWI ==================================

# Discretization package
nper = 1 # Number of stress periods
perlen_yr = 500 # Length of stress periods (years)
stplen_obj_day = 200 # Desired length of time steps (days)
steady = False # Type of simulation (transient)
tsmult = 1 # Time step multiplier

perlen = int(perlen_yr*365.25*24*3600) # Convert perlen from years to seconds
stplen_obj = stplen_obj_day*24*3600 # Convert stplen_obj from days to seconds
nstp = round(perlen/stplen_obj) # Number of time steps = length of stress period / length of time steps
stplen = round( (perlen/nstp)/(24*3600) , 2) # Real length of time steps (days)

# Saltwater Intrusion package
nsolver = 2 # (PCG solver)
toeslope = 0.16
tipslope = toeslope
nu = [0, epsilon_corr] # [0, density ratio]

zeta_init = -40 * h_noswi # Ghyben-Herzberg of h_no_SWI
zeta_init[ zeta_init < botm ] = botm # Avoid SWI to go below aquifer bottom
zeta_init = [zeta_init] # List of initial FW-SW interface elevations
isource = np.zeros((nrow, ncol), np.int) # Sources/sinks have same fluid density as active zone at the top of the aquifer 
isource[ sea == 1 ] = -2 # Ocean bottom: infiltrating water = SW, exfiltrating water = same type as water at the top of the aquifer

# Output Control package
spd = {(0, nstp-1): ['save head']}#, 'save budget']} # Save head and budget at the end of the simulation only
#output_nstp = 30 # Output frequency for zeta 
#spd = {}
#for kstp in range(0, nstp, output_nstp): # For istp ranging from 0 to nstp, with step = output_nstp
#    spd[ (0, kstp) ] = ['save head']#, 'save budget']  # Creat dictionary = {(stress period i, timestep i): save head}

#%%------------------------------- Build model --------------------------------

name_swi = name_simu + '_swi'
ipakcb = None # Don't save cell-by-cell budget data

# New model
mf = flopy.modflow.Modflow(modelname = name_swi, 
                           model_ws = model_ws, 
                           exe_name = modflow_path)

# Discretization package
dis = flopy.modflow.ModflowDis(mf, nlay = nlay, nrow = nrow, ncol = ncol, 
                               nper = nper, delr = delr, delc = delc, 
                               top = top, botm = botm, 
                               perlen = perlen, nstp = nstp, 
                               tsmult = tsmult, steady = steady, 
                               itmuni = itmuni, lenuni = lenuni)

# Set spatial reference info (after defining dis package)
mf.modelgrid.set_coord_info(xoff = xmin, yoff = ymin, angrot = 0, epsg = 2946)

# Basic package
bas = flopy.modflow.ModflowBas(mf, ibound = ibound, strt = h_noswi) 
# Recharge package
rch = flopy.modflow.ModflowRch(mf, rech = rech, ipakcb = ipakcb)
# Layer-Property Flow package
lpf = flopy.modflow.ModflowLpf(mf, laytyp = laytyp, hk = hk_array, 
                               sy = sy, ss = ss)
# General-Head Boundary package
ghb = flopy.modflow.ModflowGhb(mf, stress_period_data = ghb_sp_data, 
                               ipakcb = ipakcb)
# Well package
wel = flopy.modflow.ModflowWel(mf, stress_period_data = wel_sp_data, 
                               ipakcb = ipakcb)
# Revised Multi-Node Well package
mnw2 = flopy.modflow.ModflowMnw2(mf, mnwmax = active_mnw2_wells, 
                                 node_data = node_data, 
                                 stress_period_data = wel_sp_data_mnw2, 
                                 itmp = [active_mnw2_wells] * nper, 
                                 ipakcb = ipakcb)
# Multi-Node Well Information Package
mnwi = flopy.modflow.ModflowMnwi(mf, byndflag = 36, 
                                 mnwobs = active_mnw2_wells,
                                 wellid_unit_qndflag_qhbflag_concflag = mnwi_list)
# Output control package
oc = flopy.modflow.ModflowOc(mf, stress_period_data = spd)
# Preconditioned Conjugate-Gradient package
pcg = flopy.modflow.ModflowPcg(mf)
# Saltwater-Intrusion package
swi = flopy.modflow.ModflowSwi2(mf, nsrf = 1, istrat = 1, zeta = zeta_init, 
                                toeslope = toeslope, tipslope = tipslope, 
                                nu = nu, alpha = 0.1, beta = 0.1, 
                                ssz = ne, isource = isource, nsolver = nsolver, 
                                iswizt = 55, ipakcb = ipakcb)
# Write input files
mf.write_input()

#%%-------------------------------- Run model ---------------------------------

mf.run_model(silent=False)

#%%------------------------------- Read outputs -------------------------------

#mf = flopy.modflow.Modflow.load(os.path.join(model_ws, name_simu + '.nam'), verbose=True, forgive=True)

#---- Freshwater heads

hds = flopy.utils.binaryfile.HeadFile(
        os.path.join(model_ws, name_swi + '.hds'))

times = hds.get_times() #[i/(365.25*24*3600) for i in times]
h = hds.get_alldata()[:,:,:,:]

#---- MNW2 data

bynd_calib_swi = pd.read_csv(os.path.join(model_ws, name_swi + '.bynd'), 
                             header = 0, delim_whitespace=True,  
                             index_col = 0, usecols = [0,5,7,8], 
                             names = ['name', 'totim', 'hwell', 'hcell'])

bynd_calib_swi.index = bynd_calib_swi.index.str.lower()
hwell_calib_last = bynd_calib_swi['hwell'].tail(9).to_dict()

#---- Interface data

zta = flopy.utils.binaryfile.CellBudgetFile(
        os.path.join(model_ws, name_swi + '.zta'))

kstpkper = zta.get_kstpkper() # List of unique stress periods & timesteps (nper,nstp)
zeta = []
for kk in kstpkper:
    zeta.append(zta.get_data(kstpkper = kk, text = 'ZETASRF  1')[0])
zeta = np.array(zeta)


#%%--------------------------- Post-process outputs ---------------------------

#---- Post-process head and zeta outputs

h_P = copy.deepcopy(h) # Freshwater heads
zeta_P = copy.deepcopy(zeta) # Interface elevations

# Assign Nan value to inactive model cells
for i in range(len(spd)):
    h_P[i, 0, :, :][ ibound == 0 ] = np.nan
    zeta_P[i, 0, :, :][ ibound == 0 ] = np.nan

#---- Compute thickness of the freshwater lens (m)

thk_lens = (h_P - zeta_P)

#---- Compute volume of the lens (m3)

vlens = []
for i in range(len(thk_lens)):
    fw = np.nansum(thk_lens[i, 0, :, :]) * delr * delc # sum(hswi - zeta) * delr * delc
    vlens.append(fw)


#%%------------------------- Prepare PESTPP-OPT file --------------------------

tstep_save = -1 # We're interested in the new steady state (i.e. last timestep)


#---- Create dataframes with the data simulated at observation/forecast locations

# Head observations

hwells_df = fmod.sim_to_df(obs_h_row_col, h, tstep_save) # contains hcell values
for well in muni_wells_row_col.keys(): # Municipal wells: save hwell instead of hcell
    hwells_df.loc[hwells_df['name'] == well, 'value'] = hwell_calib_last[well]


# Interface observations

ztdem_df = fmod.sim_to_df(tdem_row_col, zeta, tstep_save) # TDEM
zert_df = fmod.sim_to_df(ert_row_col, zeta, tstep_save) # ERT

zwells_all = fmod.sim_to_df(obs_z_row_col, zeta, tstep_save) # Wells: contains zcell values

zwells_df = zwells_all[zwells_all['name'].str.startswith('pz')] # Observations
zmuni_df = zwells_all[zwells_all['name'].str.startswith('muni')] # zcell forecasts

zmuni_pct_df = copy.deepcopy(zmuni_df) # z_1pct forecasts: save z_1pct instead of zcell
#for well in zmuni_pct_df.name:
#    zmuni_pct_df.loc[zmuni_pct_df.name == well, 'value'] = z_1pct[well]


# Volume of the lens (m3)

vlens_df = pd.DataFrame({'name': 'vlens', 'value': vlens[tstep_save]}, index = [0])


#---- Drop several observations

# Inconsistent TDEM observations
ztdem_df.drop(ztdem_df[ (ztdem_df['name'] == 'tdem07') | 
        (ztdem_df['name'] == 'tdem23') | (ztdem_df['name'] == 'tdem67') | 
        (ztdem_df['name'] == 'tdem82') ].index, inplace=True)

# ERT observations close (< 100 m) to pumping wells
#zert_df.drop(zert_df[ (zert_df['name'] == 'zert08') | 
#        (zert_df['name'] == 'zert62') | (zert_df['name'] == 'zert71') | 
#        (zert_df['name'] == 'zert72') | (zert_df['name'] == 'zert73') | 
#        (zert_df['name'] == 'zert77') | (zert_df['name'] == 'zert78') | 
#        (zert_df['name'] == 'zert79') ].index, inplace=True)


#---- Assign observation name for the model

# Model observations
hwells_df.name = 'h' + hwells_df.name.astype(str)
zwells_df.name = 'z' + zwells_df.name.astype(str)
zmuni_df.name = 'z' + zmuni_df.name.astype(str)
ztdem_df.name = 'z' + ztdem_df.name.astype(str)

# Model forecasts
zmuni_pct_df.name = 'z' + zmuni_pct_df.name.astype(str) + '_pct'


#---- Export simulated data to .dat files

# Observations

fmod.write_df_to_dat(os.path.join(path_model_outputs, 'hwells'), 
                     hwells_df[['name', 'value']])

fmod.write_df_to_dat(os.path.join(path_model_outputs, 'zwells'), 
                     zwells_df[['name', 'value']])

fmod.write_df_to_dat(os.path.join(path_model_outputs, 'ztdem'), 
                     ztdem_df[['name', 'value']])

fmod.write_df_to_dat(os.path.join(path_model_outputs, 'zert'), 
                     zert_df[['name', 'value']])

# Forecasts

fmod.write_df_to_dat(os.path.join(path_model_outputs, 'zmuni'), 
                     zmuni_df[['name', 'value']])

fmod.write_df_to_dat(os.path.join(path_model_outputs, 'zmuni_pct'), 
                     zmuni_pct_df[['name', 'value']])

fmod.write_df_to_dat(os.path.join(path_model_outputs, 'vlens'), 
                     vlens_df[['name', 'value']])

