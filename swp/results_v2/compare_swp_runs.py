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

# Run with no climate change
path_swp_control_file_nocc = os.path.join('noclimatechange', 'pestpp_swp_nocc.pst')
path_swp_out_file_nocc = os.path.join('noclimatechange', 'sweep_out.csv')

# Run with climate change
path_swp_control_file_cc = os.path.join('climatechange', 'pestpp_swp_cc.pst')
path_swp_out_file_cc = os.path.join('climatechange', 'sweep_out.csv')

# IES run
path_ies = os.path.join('..', '..', 'ies', 'results_v4', 'num_reals_200')
name_pst_ies = 'pestpp_ies'


#%%---------------------------- Read control files ----------------------------

# Run with no climate change
pst = pyemu.Pst(path_swp_control_file_nocc)
fcstnames = [s for s in pst.observation_data.obsnme if s.startswith('zmuni') and not s.endswith('pct')]
obsens_nocc = pd.read_csv(path_swp_out_file_nocc, index_col = 1)[fcstnames]


# Run with climate change
pst = pyemu.Pst(path_swp_control_file_cc)
obsens_cc = pd.read_csv(path_swp_out_file_cc, index_col = 1)[fcstnames]


# IES run
no_iterations = 2
obsens_post_ies = pd.read_csv(os.path.join(path_ies, name_pst_ies + '.{0}.obs.csv'.format(no_iterations)),
                              index_col = 0) # Posterior observation ensemble

# Right-hand side of the constraint
obs = pst.observation_data
c = obs.loc[obs.obsnme.apply(lambda x: x.endswith('pct')), 'obsnme']
cons_rhs = obs.loc[c, 'obsval'].to_frame()


#%%---------------------------- Plotting functions ----------------------------

fs = 11
plt.rc('font', family = 'serif', size = fs)

def cm2inch(*tupl):
    ''' From a tuple containing centimeters, return a tuple containing inches '''
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


#%%------------------------ Prepare violin plot figure ------------------------

fcst_nocc = []
fcst_cc = []
for forecast in fcstnames:
    fcst_nocc.append(obsens_nocc.loc[:, forecast])
    fcst_cc.append(obsens_cc.loc[:, forecast])
fcst_nocc = tuple(fcst_nocc)
fcst_cc = tuple(fcst_cc)


fcst_nocc_stats = pd.DataFrame(columns = ['mean', 'median', 'stdev'])
fcst_cc_stats = pd.DataFrame(columns = ['mean', 'median', 'stdev'])

for forecast in fcstnames:
    data_nocc = [np.mean(obsens_nocc.loc[:, forecast]), 
                 np.median(obsens_nocc.loc[:, forecast]), 
                 np.std(obsens_nocc.loc[:, forecast], ddof = 1)]
    fcst_nocc_stats.loc[forecast,:] = data_nocc
    
    data_cc = [np.mean(obsens_cc.loc[:, forecast]), 
                 np.median(obsens_cc.loc[:, forecast]), 
                 np.std(obsens_cc.loc[:, forecast], ddof = 1)]
    fcst_cc_stats.loc[forecast,:] = data_cc


#%%------------------------------- FIGURE PAPER -------------------------------

# Violin plots with/without climate change

dic_props = dict(color = 'k', lw =  0.5)
meanprops = dict(marker = 'o', mfc = 'none', ms = 3, lw = 0.5, mec = 'k')
flierprops = dict(marker = 'o', mfc = 'none', ms = 3, lw = 0.1, mec = 'k')
showpoints = [True, True]
color_list = ['tab:blue', 'tab:red']
xlabels = ['No climate change', 'Climate change']


fig, ax = plt.subplots(figsize = cm2inch(19, 8))
xticks_list = list(range(1, len(fcstnames)+1))
labels = [s.replace('zmuni', 'Well ') for s in fcstnames]

# White violin plots no CC

v1_white = ax.violinplot(fcst_nocc, widths = 0.8, points = 100, 
                         positions = np.arange(1, len(fcst_nocc) + 1),
                         showmeans = False, showextrema = False, showmedians = False)
for b in v1_white['bodies']:
    m = np.mean(b.get_paths()[0].vertices[:, 0]) # get the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m) # modify the paths to not go further right than the center
    b.set_color('white')
    b.set_alpha(1)

# White violin plots CC

v2_white = ax.violinplot(fcst_cc, widths = 0.8, points = 100, 
                         positions = np.arange(1, len(fcst_nocc) + 1),
                         showmeans = False, showextrema = False, showmedians = False)
for b in v2_white['bodies']:
    m = np.mean(b.get_paths()[0].vertices[:, 0]) # get the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m) # modify the paths to not go further right than the center
    b.set_color('white')
    b.set_alpha(1)

# Color violin plots no CC

v1 = ax.violinplot(fcst_nocc, widths = 0.8, points = 100, 
                   positions = np.arange(1, len(fcst_nocc) + 1),
                   showmeans = False, showextrema = False, showmedians = True)
for b in v1['bodies']:
    m = np.mean(b.get_paths()[0].vertices[:, 0]) # get the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m) # modify the paths
    b.set_color('blue')
    b.set_alpha(0.4)

# Plot median no CC

cm = v1['cmedians']
for i in range(len(cm.get_paths())):
    p = cm.get_paths()[i].vertices[:, 0]
    m = np.mean(p) # get the center
    p[0] -= 0.1
    cm.get_paths()[i].vertices[:, 0] = np.clip(p, -np.inf, m) # modify the paths
cm.set_color('blue')

# Color violin plots CC

v2 = ax.violinplot(fcst_cc, widths = 0.8, points = 100, 
                   positions = np.arange(1, len(fcst_cc) + 1), 
                   showmeans = False, showextrema = False, showmedians = True)
for b in v2['bodies']:
    m = np.mean(b.get_paths()[0].vertices[:, 0]) # get the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf) # modify the paths
    b.set_color('red')
    b.set_alpha(0.4)

# Plot median CC

cm = v2['cmedians']
for i in range(len(cm.get_paths())):
    p = cm.get_paths()[i].vertices[:, 0]
    m = np.mean(p) # get the center
    p[1] += 0.1
    cm.get_paths()[i].vertices[:, 0] = np.clip(p, m, np.inf) # modify the paths
cm.set_color('red')

# Well botm

for ii in xticks_list:
    cons = 'zmuni' + str(ii) + '_pct'
    plt.hlines(cons_rhs.loc[cons, 'obsval'], xmin = ii - 0.25, xmax = ii + 0.25, 
               color = 'k', linestyles = '--', linewidths = 1)

# Show points

for i, tick in enumerate(xticks_list):
    y = fcst_nocc[i]
    x = np.random.normal(tick, 0.04, size = len(y))
    plt.plot(x, y, 'b.', alpha = 0.1)
    y = fcst_cc[i]
    x = np.random.normal(tick, 0.04, size = len(y))
    plt.plot(x, y, 'r.', alpha = 0.1)


ax.legend([v1['bodies'][0], v2['bodies'][0]], 
          ['No climate change', 'Climate change'],
          loc = 'lower right')

ax.text(1, -16, r'$\zeta_{50\%}$'' Well bottom')

plt.xticks(xticks_list, labels, rotation = 0)
plt.yticks(np.arange(-90, 0, 10))
ax.yaxis.grid()
plt.ylabel('Elevation (masl)')
plt.tight_layout()
plt.savefig(os.path.join('./', '00paper_violin_all.png'), dpi = 300)

