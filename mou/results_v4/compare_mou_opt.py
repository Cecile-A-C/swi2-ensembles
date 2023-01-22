# Author: Cecile Coulon

import os
working_dir = './'
os.chdir(working_dir)
modflow_path = 'mf2005'

#%%------------------------- load paths and libraries -------------------------

import pyemu
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# Folders with inputs & outputs

name_mou_nocc = 'pestpp_mou_nocc'
path_mou_nocc = 'noclimatechange'

name_mou_cc = 'pestpp_mou_cc'
path_mou_cc = 'climatechange'

# SPL-FOSM run with no climate change
name_opt_nocc_fosm = 'pestpp_opt'
path_opt_nocc_fosm = os.path.join('pestpp-opt run')

# Water demand projections
water_demand_2050 = list(pd.read_csv('water_demand_GE.csv', header = None)[0]) # m3/day
water_demand_2020 = 93 # m3/day


#%%----------------------- Define MOU reading functions -----------------------

#------ Read pareto archive summary csv to dataframe

def read_pareto_archive(path_mou, name_mou):
    
    pareto_archive_df =  pd.read_csv(
        os.path.join(path_mou, '{0}.pareto.archive.summary.csv'.format(name_mou)))
    pareto_archive_df = pareto_archive_df[pareto_archive_df['nsga2_front'] == 1] # 1st nondominated front    
    pareto_archive_df.loc[:, 'memgen'] = pareto_archive_df.member.apply(lambda x: int(x.split('_')[0].split('=')[1]))
    
    print('read', '{0}.pareto.archive.summary.csv'.format(name_mou))
    return pareto_archive_df

#------ Read archive decision variable population files to dictionary

def read_archive_pop(path_mou, name_mou, gens_list):

    archive_pop_dict = {}
    for gen in gens_list:
        try:
            archive_pop_dict[gen] = pd.read_csv(
                os.path.join(path_mou, '{0}.{1}.archive.dv_pop.csv'.\
                             format(name_mou, gen)))[dvars_list + ['real_name']]
        except:
            pass
            # print('{0}.{1}.archive.dv_pop.csv'.format(name_mou, gen), 'does not exist')

    return archive_pop_dict

#------ Merge objective function & dvar information contained in archive files

def merge_archive_info(pareto_archive_df, archive_pop_dict):
    
    pareto_archive_dict = {}
    for gen in pareto_archive_df.generation.unique():
        try:
            pareto_archive_dict[gen] = pareto_archive_df.loc[pareto_archive_df.generation == gen]
            pareto_archive_dict[gen] = pareto_archive_dict[gen].merge(archive_pop_dict[gen], 
                how = 'left', left_on = ['member'], right_on = ['real_name']) # merge with dv values
            pareto_archive_dict[gen].pop('real_name')
        except:
            pass

    return pareto_archive_dict


#------ Append the restart MOU run to the initial MOU run: dictionaries

def merge_dictionnaries(init_dict, restart_dict):

    restart_dict_cp = restart_dict.copy()
    restart_dict_cp.pop(0, None)
    restart_dict_cp = dict(zip([key + list(init_dict.keys())[-1] \
                                  for key in restart_dict_cp.keys()], 
                            list(restart_dict_cp.values())))
    
    dict_all = init_dict.copy()
    dict_all.update(restart_dict_cp)
    
    return dict_all

#------ Read nested parameter stack csv files to dictionary

def read_nest_par_stack(path_mou, name_mou, gens_list):

    nested_par_stack_dict = {}
    for gen in gens_list:
        try:
            nested_par_stack_dict[gen] = pd.read_csv(
                os.path.join(path_mou, '{0}.{1}.nested.par_stack.csv'.\
                             format(name_mou, gen)), index_col = 0)
            
            nested_par_stack_dict[gen].loc[:, 'member'] = \
                nested_par_stack_dict[gen].index.map(lambda x: x.split('||')[1])
        except:
            pass
        
    return nested_par_stack_dict

#------ Read nested observation stack csv files to dictionary

def read_nest_obs_stack(path_mou, name_mou, gens_list):

    nested_obs_stack_dict = {}
    for gen in gens_list:
        try:
            nested_obs_stack_dict[gen] = pd.read_csv(
                os.path.join(path_mou, '{0}.{1}.nested.obs_stack.csv'.\
                             format(name_mou, gen)), index_col = 0)
                
            nested_obs_stack_dict[gen].loc[:, 'member'] = \
                nested_obs_stack_dict[gen].index.map(lambda x: x.split('||')[1])
        except:
            pass
    
    return nested_obs_stack_dict

#------ Read observation population csv files to dictionary

def read_obs_pop(path_mou, name_mou, gens_list):

    obs_pop_dict = {}
    for gen in gens_list:
        try:
            obs_pop_dict[gen] = pd.read_csv(
                os.path.join(path_mou, '{0}.{1}.obs_pop.csv'.\
                             format(name_mou, gen)), index_col = 0)
        except:
            pass
            # print('{0}.{1}.obs_pop.csv'.format(name_mou, gen), 'does not exist')
    
    return obs_pop_dict

#------ Read chance observation population csv files to dictionary

def read_chance_obs_pop(path_mou, name_mou, gens_list):

    chance_obs_pop_dict = {}
    for gen in gens_list:
        try:
            chance_obs_pop_dict[gen] = pd.read_csv(
                os.path.join(path_mou, '{0}.{1}.chance.obs_pop.csv'.\
                             format(name_mou, gen)), index_col = 0)
        except:
            pass
            # print('{0}.{1}.chance.obs_pop.csv'.format(name_mou, gen), 'does not exist')

    return chance_obs_pop_dict

#------ Read MOU output files

def read_mou_outputs(path_mou, name_mou, gens_list):
    
    pareto_archive_df = read_pareto_archive(path_mou, name_mou) # pareto archive summary
    archive_pop_dict = read_archive_pop(path_mou, name_mou, gens_list) # archive dvar pop files
    pareto_archive_dict = merge_archive_info(pareto_archive_df, archive_pop_dict) # merge info
    
    nested_par_stack_dict = read_nest_par_stack(path_mou, name_mou, gens_list) # nested par stacks
    nested_obs_stack_dict = read_nest_obs_stack(path_mou, name_mou, gens_list) # nested obs stacks
    
    obs_pop_dict = read_obs_pop(path_mou, name_mou, gens_list) # obs pop files
    chance_obs_pop_dict = read_chance_obs_pop(path_mou, name_mou, gens_list) # chance obs pop files

    return pareto_archive_dict, nested_par_stack_dict, nested_obs_stack_dict,\
        obs_pop_dict, chance_obs_pop_dict

#------ Print total number of gens vs gens recorded in the archive files

def print_gens_in_archive(pareto_archive_dict, gens_list):
    
    print('total number of gens:', gens_list[-1],
            '\ngens in archive files:', 
            [gen for gen in gens_list if gen in pareto_archive_dict.keys()], 
            ' - total:', 
            len([gen for gen in gens_list if gen in pareto_archive_dict.keys()]))

#------ List the members of the final generation

def list_final_members(pareto_archive_dict):
    
    lastgen = list(pareto_archive_dict.keys())[-1]
    list_mem = pareto_archive_dict[lastgen]['member'].tolist()
    list_mem.sort()
    
    print('no. of members at last gen:', len(list_mem))
    return list_mem

#------ Get pestpp-opt output dvars and risks

def get_opt_dvars_risk(name_pst, folder, dvar_nme, list_risk):
    
    ''' Create a dataframe containing optimal decision variables & decision 
    objective function for each risk stance '''
    
    par_dic = {}
    
    # Retrieve decision vars from all .par files
    
    for risk in list_risk:
        
        parfile = os.path.join(folder, risk.replace('.', '_') + '_' + name_pst + '.par')
        par_df = pyemu.pst_utils.read_parfile(parfile)
        par_df = par_df.loc[dvar_nme, :]
        
        # When the solution is infeasible, pestpp-opt writes extreme negative values to the par file:
        if par_df.parval1.sum() < 0: 
            print('infeasible at risk ' + risk)
            continue
        
        par_dic[risk] = par_df.parval1
    
    # Process the dvar dictionary for plotting
    
    par_df = pd.concat(par_dic, axis = 1).T
    par_df.index = [int(float(idx) * 100) for idx in par_df.index]
    par_df['phi'] = par_df.sum(axis=1) # Compute decision objective function
    
    return par_df


#%%---------------------- Read PESTPP-MOU control files -----------------------

pst_cc = pyemu.Pst(os.path.join(path_mou_cc, name_mou_cc + '.pst'))

obs = pst_cc.observation_data
par = pst_cc.parameter_data

dvars_list = [x for x in pst_cc.par_names if x.startswith('qmuni')]
obj_list = pst_cc.pestpp_options['mou_objectives'].split(',')

fcst_list = [s for s in obs.obsnme if s.startswith('zmuni') and not s.endswith('pct')]
cons_list = [s for s in obs.obsnme if s.endswith('pct')]
cons_rhs_df = obs.loc[cons_list, 'obsval'].to_frame()


#%%----------------- Read MOU output files, no climate change -----------------

pst_nocc = pyemu.Pst(os.path.join(path_mou_nocc, name_mou_nocc + '.pst'))
print('read', '{0}.pst'.format(name_mou_nocc))
gens_list_nocc = list(np.arange(0, pst_nocc.control_data.noptmax + 1))

pareto_archive_dict_nocc, nested_par_stack_dict_nocc, nested_obs_stack_dict_nocc,\
    obs_pop_dict_nocc, chance_obs_pop_dict_nocc \
        = read_mou_outputs(path_mou_nocc, name_mou_nocc, gens_list_nocc)

# Restart files

name_mou_r =  name_mou_nocc + '_restart'
pst_r = pyemu.Pst(os.path.join(path_mou_nocc, '{0}.pst'.format(name_mou_r)))
print('read', '{0}.pst'.format(name_mou_r))
gens_list_r = list(np.arange(0, pst_r.control_data.noptmax + 1))

pareto_archive_dict_r, nested_par_stack_dict_r, nested_obs_stack_dict_r,\
    obs_pop_dict_r, chance_obs_pop_dict_r \
        = read_mou_outputs(path_mou_nocc, name_mou_r, gens_list_r)

# Merge

gens_list_nocc = list(np.arange(0, pst_nocc.control_data.noptmax + pst_r.control_data.noptmax + 1))

pareto_archive_dict_nocc = merge_dictionnaries(pareto_archive_dict_nocc, pareto_archive_dict_r) # pareto archive files
nested_par_stack_dict_nocc = merge_dictionnaries(nested_par_stack_dict_nocc, nested_par_stack_dict_r) # nested par stacks
nested_obs_stack_dict_nocc = merge_dictionnaries(nested_obs_stack_dict_nocc, nested_obs_stack_dict_r) # nested obs stacks
obs_pop_dict_nocc = merge_dictionnaries(obs_pop_dict_nocc, obs_pop_dict_r) # obs pop files
chance_obs_pop_dict_nocc = merge_dictionnaries(chance_obs_pop_dict_nocc, chance_obs_pop_dict_r) # chance obs pop files

print_gens_in_archive(pareto_archive_dict_nocc, gens_list_nocc) # total number of gens vs gens recorded in the archive files
list_mem_nocc = list_final_members(pareto_archive_dict_nocc) # list of final members

# highest-reliability scenario
gen = list(pareto_archive_dict_nocc.keys())[-1]
highestreliab_nocc = pareto_archive_dict_nocc[gen].loc[pareto_archive_dict_nocc[gen]['_risk_'].idxmax()]
print(highestreliab_nocc)

# Lastgen run

pareto_archive_dict_nocc_lastgen, nested_par_stack_dict_nocc_lastgen, nested_obs_stack_dict_nocc_lastgen,\
    obs_pop_dict_nocc_lastgen, chance_obs_pop_dict_nocc_lastgen \
        = read_mou_outputs(os.path.join(path_mou_nocc, 'mou_restart_lastgen'), 
                           '{0}_restart'.format(name_mou_nocc), [0])


#%%------------------ Read MOU output files, climate change -------------------

pst_cc = pyemu.Pst(os.path.join(path_mou_cc, name_mou_cc + '.pst'))
print('read', '{0}.pst'.format(name_mou_cc))
gens_list_cc = list(np.arange(0, pst_cc.control_data.noptmax + 1))

pareto_archive_dict_cc, nested_par_stack_dict_cc, nested_obs_stack_dict_cc,\
    obs_pop_dict_cc_cc, chance_obs_pop_dict_cc \
        = read_mou_outputs(path_mou_cc, name_mou_cc, gens_list_cc)

print_gens_in_archive(pareto_archive_dict_cc, gens_list_cc) # total number of gens vs gens recorded in the archive files
list_mem_cc = list_final_members(pareto_archive_dict_cc) # list of final members

# highest-reliability scenario
gen = list(pareto_archive_dict_cc.keys())[-1]
highestreliab_cc = pareto_archive_dict_cc[gen].loc[pareto_archive_dict_cc[gen]['_risk_'].idxmax()]
print(highestreliab_cc)

# Lastgen run

pareto_archive_dict_cc_lastgen, nested_par_stack_dict_cc_lastgen, nested_obs_stack_dict_cc_lastgen,\
    obs_pop_dict_cc_lastgen, chance_obs_pop_dict_cc_lastgen \
        = read_mou_outputs(os.path.join(path_mou_cc, 'mou_restart_lastgen'), 
                           '{0}_restart'.format(name_mou_cc), [0])


#%%---------------------- Read PESTPP-OPT control files -----------------------

# Post calibration risk values
list_risk_fosm = sorted([0.01] + list(np.arange(0.05, 1.00, 0.05)) + [0.96] + [0.97] + [0.98])
list_risk_fosm = ['%.2f' % x for x in list_risk_fosm]

par_df_fosm = get_opt_dvars_risk(name_opt_nocc_fosm, path_opt_nocc_fosm, dvars_list, list_risk_fosm) # post calibration
par_df_fosm_m3day = par_df_fosm * 24 * 3600 # Convert m3/s to m3/day


#%%------------------------ Define plotting functions -------------------------

fs = 11
plt.rcParams['font.size'] = fs
plt.rc('font', family = 'serif', size = fs)

xaxis = [1] + list(np.arange(5, 100, 5)) + [99] # x-axis for all figures
xaxis_bis = [1] + list(np.arange(10, 100, 10)) + [99]

color = {dvars_list[i]: 
    plt.rcParams['axes.prop_cycle'].by_key()['color'][i] for i in range(len(dvars_list))}


def cm2inch(*tupl):
    ''' From a tuple containing centimeters, return a tuple containing inches '''
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


#%%-------------------------------FIGURE PAPER --------------------------------

# Pareto front

def plot_pumping_vs_probability_demand(pareto_archive_dict_nocc, pareto_archive_dict_cc):

    fig, ax = plt.subplots(1, 1, figsize = cm2inch(14, 12))
    fig.subplots_adjust(wspace=0.2, hspace=0.2, left=0.125, right=1, bottom=0, top=1)
    
    #------- X-AXIS 1
    
    # MOU run, no climate change
    gen = list(pareto_archive_dict_nocc.keys())[-1]
    x = pareto_archive_dict_nocc[gen]['_risk_']
    y = pareto_archive_dict_nocc[gen]['pump_rate']
    ax.scatter((1 - x) * 100, y, color = 'blue', alpha = 0.6, 
                label = 'No climate change')
    
    # MOU run, with climate change
    gen = list(pareto_archive_dict_cc.keys())[-1]
    x = pareto_archive_dict_cc[gen]['_risk_']
    y = pareto_archive_dict_cc[gen]['pump_rate']
    ax.scatter((1 - x) * 100, y, color = 'red', alpha = 0.6, 
                label = 'Climate change')
    
    # FOSM-based OPT run, no climate change
    ax.plot(100 - par_df_fosm_m3day.index.values, par_df_fosm_m3day.phi, marker = '.', ms = 7, 
            lw = 0.3, ls = '-', c = 'blue', 
            label = 'No climate change \n(FOSM-based opt.)')
    
    # Legend
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [1,2,0]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
               loc = 'upper left')
    
    #------- X-AXIS 2
    
    ax2 = ax.twiny()
    
    # Projected water demand
    ax2.hist(water_demand_2050, bins = 30, orientation = 'horizontal', 
             density = True, color = 'green', alpha = 0.3, label = 'Projected water demand')
    
    # Current water demand
    ax2.axhline(water_demand_2020, color = 'green', alpha = 0.8, label = 'Current water demand')
    
    # Legend
    ax2.legend(bbox_to_anchor = (0.582, 0.74))
    
    #------- PARAMETERS
    
    # x-axis 1
    ax.set_xticks(xaxis)
    ax.set_xticklabels([str(x) for x in xaxis], fontsize = fs - 1)
    ax.set_xlim([0, 100])
    ax.set_xlabel('Probability of well salinization (%)', fontsize = fs)
    
    # x-axis 2
    ax2.set_xticks([0, 0.01, 0.02, 0.03, 0.04, 0.05])
    ax2.set_xticklabels([str(x) for x in [0, 0.01, 0.02, 0.03, 0.04, 0.05]], 
                        fontsize = fs - 1, color = 'green')
    ax2.set_xlabel('Probability density', fontsize = fs, color = 'green')
    
    # Others
    ax.set_axisbelow(True)
    ax.grid(axis = 'both', lw = 0.3)
    ax.set_ylabel('Total pumping rate (m$^{3}$/day)', fontsize = fs)
    ax.set_yticks(list(np.arange(0, 1100, 100)))
    ax.set_ylim([0, 700])
    plt.tight_layout()

plot_pumping_vs_probability_demand(pareto_archive_dict_nocc, pareto_archive_dict_cc)
plt.savefig(os.path.join('./', '00paper_pumping_vs_probability_demand.png'), dpi = 300)


#%%-------------------------------FIGURE PAPER --------------------------------

# Allocation of pumping

def plot_dvars_vs_probability(pareto_archive_dict_nocc, pareto_archive_dict_cc):
    
    xaxis2 = [1] + list(np.arange(10, 100, 10)) + [99]
    
    fig, ax = plt.subplots(5, 2, figsize = cm2inch(19, 24), sharex = True, sharey = True)
    fig.subplots_adjust(wspace = 0, hspace = 0)
    ax = ax.flatten()
    fig.delaxes(ax[9])
    
    for i, dv in enumerate(dvars_list):
        
        # FOSM-based OPT run, no climate change
        ax[i].plot(100 - par_df_fosm_m3day.index.values, par_df_fosm_m3day[dv], 
                marker = '.', ms = 5, lw = 0.3, ls = '-', c = 'blue', 
                )#label = 'No climate change\n(FOSM-based)')
        
        # MOU run, no climate change
        gen = list(pareto_archive_dict_nocc.keys())[-1]
        x = pareto_archive_dict_nocc[gen]['_risk_']
        y = pareto_archive_dict_nocc[gen][dv]
        ax[i].scatter((1 - x) * 100, y, color = 'blue', alpha = 0.5, 
                    )#label = 'No climate change\n(ensemble-based)')
        
        # MOU run, with climate change
        gen = list(pareto_archive_dict_cc.keys())[-1]
        x = pareto_archive_dict_cc[gen]['_risk_']
        y = pareto_archive_dict_cc[gen][dv]
        ax[i].scatter((1 - x) * 100, y, color = 'red', alpha = 0.5,
                    )#label = 'Climate change\n(ensemble-based)')
        
        ax[i].scatter(0, 0, s = 0.001, color = 'None', label = 'Well no. ' + dv.split('qmuni')[-1])
        ax[i].legend(loc = 'upper left', handletextpad = 0)
        
        ax[i].grid(axis = 'both', lw = 0.3)
    
    ax[0].set_xlim([0, 100])
    ax[0].set_xticks(xaxis2)
    ax[0].set_yticks(list(np.arange(0, 160, 20)))
    ax[0].set_ylim(-5, 150)
    ax[8].set_xlabel( 'Probability of well salinization (%)', fontsize = fs)
    
    fig.text(0.01, 0.5, 'Pumping rate (m$^{3}$/day)', va = 'center', rotation = 'vertical', fontsize = fs)
    plt.tight_layout(pad = 1.08, h_pad = 0.2, w_pad = 0, rect = (0.02, 0.0, 1, 1))


plot_dvars_vs_probability(pareto_archive_dict_nocc, pareto_archive_dict_cc)
plt.savefig(os.path.join('./', '00paper_dvars_vs_probability.png'), dpi = 300)
