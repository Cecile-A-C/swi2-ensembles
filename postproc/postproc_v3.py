# Author: Cecile Coulon

import os
working_dir = './'
os.chdir(working_dir)

#%%------------------------- load paths and libraries -------------------------

import itertools, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import scipy.stats as stats

# Folders with inputs & outputs

num_it = 2 # number of IES iteration to use
base = True # There is a base realization in the parens
ies_control_file = 'pestpp_ies'
path_ies_out = os.path.join('..', 'ies', 'results_v4')

path_postproc_out = 'results_v3'
if not os.path.isdir(path_postproc_out):
    os.mkdir(path_postproc_out)


#%% Plotting functions

fs = 11
plt.rc('font', family = 'serif', size = fs)

def cm2inch(*tupl):
    ''' From a tuple containing centimeters, return a tuple containing inches '''
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def plot_hist(ens, mu, stdev, fmt, color, vline=None, vline_label=None):
    ''' Function to plot histogram of an ensemble + mean + standard deviation (+ optional) vertical lines
        ens = ensemble (list)
        mu = mean of the ensemble (float)
        stdev = standard deviation of the ensemble (float)
        fmt = format when printing float to string (string: ex '1.2f', '1.0f')
        color = color of the histogram (string)
        vline = vertical line to plot (float) (default: None)
        vline_label = legend of vertical line (string) (default: None) '''
    
    fig = plt.figure(figsize = cm2inch(15, 8))
    ax = fig.add_subplot(111)
    
    # Histogram
    ax.hist(ens, bins = int(np.sqrt(len(ens))), alpha = 0.6, color = color, label = 'Ensemble')
    
    # Mean + standard deviation vertical lines
    label = r'$\mu$' + ' = {:' + fmt + '} m'
    ax.axvline(mu, lw = 1.2, color = 'tab:red', label = label.format(mu))
    label = r'$\mu$ $\pm$ 2$\sigma$''\n'r'($\sigma$ = {:' + fmt + '} m)'
    ax.axvline(mu - stdev * 2, lw = 1.2, color = 'grey', label =  label.format(stdev))
    ax.axvline(mu + stdev * 2, lw = 1.2, color = 'grey', label = '_')
    
    # Additional vertical line (optional)
    if vline is not None:
        label = ' = {:' + fmt + '} m'
        ax.axvline(vline, lw = 1.2, color = 'tab:blue', label = vline_label + label.format(vline))
    
    plt.legend()
    plt.ylabel('Frequency', fontsize = fs)

dic_props = dict(color = 'k', lw =  0.5)
meanprops = dict(marker = 'o', mfc = 'none', ms = 3, lw = 0.5, mec = 'k')
flierprops = dict(marker = 'o', mfc = 'none', ms = 3, lw = 0.1, mec = 'k')

def plot_violin(ens, color, xlabel, hline=None, hline_label=None, showpoints=False):
    ''' Function to plot violin plot of an ensemble + optional horizontal line
        ens = ensemble (list)
        color = color of the histogram (string)
        hline = horizontal line to plot (float) (default: None)
        hline_label = legend of horizontal line (string) (default: None)
        showpoints = Flag determining whether to superimpose data points (Boolean)'''
    
    xticks_list = [1]
    fig, ax = plt.subplots(figsize = cm2inch(7, 8))
    
    # Boxplot
    plt.boxplot(ens, widths = 0.4, showbox = True, showcaps = False, 
                showmeans = True, showfliers = False, boxprops = dic_props, 
                whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
                meanprops = meanprops, flierprops = flierprops)
    
    # Violin plot (white)
    vp = ax.violinplot(ens, widths = 0.8, points = 100, showmeans = False, 
                       showmedians = False, showextrema = False)
    for pc in vp['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('white')
        pc.set_alpha(1)
    
    # Violin plot (color)
    vp = ax.violinplot(ens, widths = 0.8, points = 100, showmeans = False, 
                       showmedians = False, showextrema = False)
    for pc in vp['bodies']:
        pc.set_facecolor(color)
        pc.set_edgecolor(color)
        pc.set_alpha(0.4)
    
    # Additional horizontal line (optional)
    if hline is not None:
        plt.hlines(hline, xmin = 1 - 0.15, xmax = 1 + 0.15, color = 'k', 
                   linestyles = '-', linewidths = 1.5, label = hline_label)
        plt.legend(loc = 'upper right')
    
    if showpoints == True:
        x = np.random.normal(xticks_list[0], 0.04, size = len(ens))
        plt.plot(x, ens, 'k.', alpha = 0.1)
    
    ax.set_xticks(xticks_list)
    ax.set_xticklabels([xlabel])#, color = color)
    # ax.xaxis.label.set_color('red')
    ax.set_axisbelow(True)
    ax.yaxis.grid()

def plot_violin_multiple(ens_list, color_list, xlabels, showpoints=False):
    ''' Function to plot several violin plots side by side
        ens = list of lists (= the ensembles)
        color_list = list of colors (strings))
        showpoints = list of flags determining whether to superimpose data points (Booleans)'''
    
    xticks_list = list(range(1, len(ens_list) + 1))
    color_dic = dict(zip(list(range(len(ens_list))), color_list))
    
    if len(ens_list) == 2:
        figsize = cm2inch(12, 8)
    elif len(ens_list) == 4:
        figsize = cm2inch(16, 8)
    
    fig, ax = plt.subplots(figsize = figsize)
    
    # Boxplot
    plt.boxplot(ens_list, widths = 0.4, showbox = True, showcaps = False, 
                showmeans = True, showfliers = False, boxprops = dic_props, 
                whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
                meanprops = meanprops, flierprops = flierprops)
    
    # Violin plot (white)
    vp = ax.violinplot(ens_list, widths = 0.8, points = 100, showmeans = False, 
                       showmedians = False, showextrema = False)
    for pc in vp['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('white')
        pc.set_alpha(1)
    
    # Violin plot (color)
    vp = ax.violinplot(ens_list, widths = 0.8, points = 100, showmeans = False,
                       showmedians = False, showextrema = False)
    
    for i, pc in enumerate(vp['bodies']):
        pc.set_facecolor(color_dic[i])
        pc.set_edgecolor(color_dic[i])
        pc.set_alpha(0.4)
    
    if showpoints != False: # Points
        for i, tick in enumerate(xticks_list):
            if showpoints[i] == True:
                y = ens_list[i]
                x = np.random.normal(tick, 0.04, size = len(y))
                plt.plot(x, y, 'k.', alpha = 0.1)
    
    plt.xticks(xticks_list, xlabels, rotation = 0)
    ax.set_axisbelow(True)
    ax.yaxis.grid()


#%%------------------------ READ POSTERIOR IES PARAMS -------------------------

# Read last parens file from IES

parens_file = ies_control_file + '.{0}.par.csv'.format(num_it) # get last parens file
parens = pd.read_csv(os.path.join(path_ies_out, parens_file), index_col = 0)
num_reals = len(parens) # includes 'base' realization

# Historical recharge ensemble (mm/yr)

rech_ies = list(parens['rech_mmyr'])
rech_ies_mean = np.mean(rech_ies)
rech_ies_stdev = np.std(rech_ies, ddof = 1)

# Historical sea level value (m)

slvl_2020_mean = round(parens.slvl.unique()[0], 3) # 0.014 m


#%%------------------- READ 2050 SLVL & RECHARGE ENSEMBLES --------------------
# (if they've already been created)

save_ens_file = 'rech_slvl_ensembles_details.csv'
rech_slvl_ensembles = pd.read_csv(os.path.join(path_postproc_out, save_ens_file),
                                 index_col = 0)
slvl_2050 = rech_slvl_ensembles.loc[:, 'slvl']
# rech_2050 = rech_slvl_ensembles.loc[:, 'rech_mmyr']
rech_2050_v1 = rech_slvl_ensembles.loc[:, 'rech_mmyr_v1']
rech_2050_v2 = rech_slvl_ensembles.loc[:, 'rech_mmyr_v2']
delta_rech_gauss = rech_slvl_ensembles.loc[:, 'delta_rech_gauss']


#%%------------------------------ 2050 SEA LEVEL ------------------------------

#------ 2050 sea level predictions (Barnett et al 2017)

slvl_2050_mean = slvl_2020_mean + 0.19
slvl_2050_stdev = 0.11

#------ Generate the 2050 sea level ensemble
# (sampling assuming Gaussian distribution, and truncating so that no slvl_2050 < 0

lbnd, ubnd = 0, 100

if base == True:
    slvl_2050 = stats.truncnorm.rvs((lbnd - slvl_2050_mean) / slvl_2050_stdev, 
                                    (ubnd - slvl_2050_mean) / slvl_2050_stdev, 
                                    loc = slvl_2050_mean, scale = slvl_2050_stdev, 
                                    size = num_reals - 1).tolist()
    slvl_2050.append(slvl_2050_mean) # append base

else:
    slvl_2050 = stats.truncnorm.rvs((lbnd - slvl_2050_mean) / slvl_2050_stdev, 
                                    (ubnd - slvl_2050_mean) / slvl_2050_stdev, 
                                    loc = slvl_2050_mean, scale = slvl_2050_stdev, 
                                    size = num_reals).tolist()


#%%---------------------------------- FIGURE ----------------------------------

color = 'red'
fig, ax = plt.subplots(figsize = cm2inch(7.5, 7))

# Boxplot
plt.boxplot(slvl_2050, widths = 0.35, showbox = True, showcaps = False, 
            showmeans = True, showfliers = False, boxprops = dic_props, 
            whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
            meanprops = meanprops, flierprops = flierprops)

# Violin plot (white)
vp = ax.violinplot(slvl_2050, widths = 0.7, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor('white')
    pc.set_edgecolor('white')
    pc.set_alpha(1)

# Violin plot (color)
vp = ax.violinplot(slvl_2050, widths = 0.7, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor(color)
    pc.set_edgecolor(color)
    pc.set_alpha(0.4)

# Plot historical-to-current
plt.hlines(slvl_2020_mean, xmin = 1 - 0.20, xmax = 1 + 0.20, color = 'blue', 
           linestyles = '-', linewidths = 1.3)

# Show points
x = np.random.normal(1, 0.04, size = len(slvl_2050))
plt.plot(x, slvl_2050, 'k.', alpha = 0.1)

ax.set_xticks([1])
ax.set_xticklabels(['Sea level'])#, color = color)
ax.set_axisbelow(True)
ax.yaxis.grid()
plt.ylabel('Elevation (masl)')
plt.ylim(bottom = 0)
plt.tight_layout()
plt.savefig(os.path.join(path_postproc_out, '00paper_slvl2050_vp.png'), dpi = 300)


#%%------------------------------ SWB2 RECHARGE -------------------------------

# SWB2 historical recharge value (1989-2019 average) (mm/year)

rech_swb2_2019 = 524

# Read the SWB2 predicted 2021-2050 recharge ensembles (mm/year)

rech_swb2_allyears = pd.read_csv('SWB2_recharge_scenarios.csv', index_col = 0)
rech_swb2_allyears.columns = rech_swb2_allyears.columns.astype(int)

# Read the SWB2 predicted 2050 recharge ensemble (mm/year)

rech_swb2_2050 = list(rech_swb2_allyears[2050])
rech_swb2_2050_mean = np.mean(rech_swb2_2050)
rech_swb2_2050_stdev = np.std(rech_swb2_2050, ddof = 1)

rech_swb2_allyearstats = pd.DataFrame(
    {'Mean': rech_swb2_allyears.mean(), 'Std': rech_swb2_allyears.std(ddof = 1), 
     'Max': rech_swb2_allyears.max(), 'Min': rech_swb2_allyears.min()})


#%%---------------------------------- FIGURE ----------------------------------

color = 'orange'
fig, ax = plt.subplots(1, 1, figsize = cm2inch(19.5, 6))

# Plot historical-to-current
line = ax.plot(2020, rech_swb2_2019, color = 'k', marker = '.', ls = 'None', 
          label = 'Current')

# Plot projected range
ax.fill_between(rech_swb2_allyearstats.index, rech_swb2_allyearstats['Min'], 
                   rech_swb2_allyearstats['Max'], color = color, alpha = 0.4, 
                   label = 'Projected (range)')

# Plot projected min, max, mean
ax.plot(rech_swb2_allyearstats['Min'], lw = 0.5, marker = 'None', color = color)
ax.plot(rech_swb2_allyearstats['Max'], lw = 0.5, marker = 'None', color = color)
ax.plot(rech_swb2_allyearstats['Mean'], lw = 0.8, marker = '.', color = color, 
        label = 'Projected (mean)')


ax.set_xlim(2019.95, 2050.05)
ax.xaxis.set_minor_locator(MultipleLocator(1))

ax.legend(loc = 'upper left', ncol = 3, bbox_to_anchor = (-0.02, 1.27))
ax.set_ylabel('Recharge (mm/year)', fontsize = fs)
plt.tight_layout()
plt.savefig(os.path.join(path_postproc_out, '00paper_rech_swb2_allyears.png'), dpi = 300)


#%%---------------------------------- FIGURE ----------------------------------

color = 'orange'
fig, ax = plt.subplots(figsize = cm2inch(6.7, 7))

# Boxplot
plt.boxplot(rech_swb2_2050, widths = 0.4, showbox = True, showcaps = False, 
            showmeans = True, showfliers = False, boxprops = dic_props, 
            whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
            meanprops = meanprops, flierprops = flierprops)

# Violin plot (white)
vp = ax.violinplot(rech_swb2_2050, widths = 0.8, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor('white')
    pc.set_edgecolor('white')
    pc.set_alpha(1)

# Violin plot (color)
vp = ax.violinplot(rech_swb2_2050, widths = 0.8, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor(color)
    pc.set_edgecolor(color)
    pc.set_alpha(0.4)

# Plot historical-to-current
plt.hlines(rech_swb2_2019, xmin = 1 - 0.20, xmax = 1 + 0.20, color = 'k', 
           linestyles = '-', linewidths = 1.3)

# Show points
x = np.random.normal(1, 0.04, size = len(rech_swb2_2050))
plt.plot(x, rech_swb2_2050, 'k.', alpha = 0.1)

ax.set_xticks([1])
ax.set_xticklabels([r'$R_{\rm 2050,SWB2}$'])#, color = color)
ax.set_axisbelow(True)
ax.yaxis.grid()
plt.yticks(np.arange(300, 900, 100))
plt.ylabel('Recharge (mm/year)')
plt.tight_layout()
plt.savefig(os.path.join(path_postproc_out, '00paper_rech_swb2_2050_vp.png'), dpi = 300)


#%%------------------ GENERATE SWB2 delta_recharge ENSEMBLE -------------------
# (delta_recharge = R2050 - R2019)

delta_rech = [x - rech_swb2_2019 for x in rech_swb2_2050]
delta_rech_mean = np.mean(delta_rech)
delta_rech_stdev = np.std(delta_rech, ddof = 1)


#%%------------------------------- FIGURE PAPER -------------------------------

color = 'grey'
fig, ax = plt.subplots(figsize = cm2inch(6.7, 7))

# Boxplot
plt.boxplot(delta_rech, widths = 0.4, showbox = True, showcaps = False, 
            showmeans = True, showfliers = False, boxprops = dic_props, 
            whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
            meanprops = meanprops, flierprops = flierprops)

# Violin plot (white)
vp = ax.violinplot(delta_rech, widths = 0.8, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor('white')
    pc.set_edgecolor('white')
    pc.set_alpha(1)

# Violin plot (color)
vp = ax.violinplot(delta_rech, widths = 0.8, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor(color)
    pc.set_edgecolor(color)
    pc.set_alpha(0.4)

# Show points
x = np.random.normal(1, 0.04, size = len(delta_rech))
plt.plot(x, delta_rech, 'k.', alpha = 0.1)

ax.set_xticks([1])
ax.set_xticklabels([r'$\Delta R$'])#, color = color)
ax.set_axisbelow(True)
ax.yaxis.grid()
plt.yticks(np.arange(-300, 400, 100))
ax.set_ylim(-320, 360)
plt.ylabel(r'$\Delta R$ (mm/year)')
plt.tight_layout()
plt.savefig(os.path.join(path_postproc_out, '00paper_delta_rech_vp.png'), dpi = 300)


#%%------------------------- 2050 RECHARGE, OPTION 1 --------------------------

#------ Resample delta_recharge --> delta_recharge_gauss 
# (assuming Gaussian distribution)

if base == True:
    delta_rech_gauss = np.random.normal(delta_rech_mean, delta_rech_stdev, num_reals - 1).tolist()
    delta_rech_gauss.append(delta_rech_mean) # append base
else:
    delta_rech_gauss = np.random.normal(delta_rech_mean, delta_rech_stdev, num_reals).tolist()


#%%------------------------------- FIGURE PAPER -------------------------------

color = 'darkgrey'
fig, ax = plt.subplots(figsize = cm2inch(6.7, 7))

# Boxplot
plt.boxplot(delta_rech_gauss, widths = 0.4, showbox = True, showcaps = False, 
            showmeans = True, showfliers = False, boxprops = dic_props, 
            whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
            meanprops = meanprops, flierprops = flierprops)

# Violin plot (white)
vp = ax.violinplot(delta_rech_gauss, widths = 0.8, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor('white')
    pc.set_edgecolor('white')
    pc.set_alpha(1)

# Violin plot (color)
vp = ax.violinplot(delta_rech_gauss, widths = 0.8, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor(color)
    pc.set_edgecolor(color)
    pc.set_alpha(0.4)

# Show points
x = np.random.normal(1, 0.04, size = len(delta_rech_gauss))
plt.plot(x, delta_rech_gauss, 'k.', alpha = 0.1)

ax.set_xticks([1])
ax.set_xticklabels([r'$\Delta R$ (Gaussian)'])#, color = color)
ax.set_axisbelow(True)
ax.yaxis.grid()
plt.yticks(np.arange(-300, 400, 100))
ax.set_ylim(-320, 360)
plt.ylabel(r'$\Delta R$ (mm/year)')
plt.tight_layout()
plt.savefig(os.path.join(path_postproc_out, '00paper_delta_rech_gauss_vp.png'), dpi = 300)



#%%----- GENERATE rech_2050_v1 = rech_modflow_ies + delta_recharge_gauss ------

rech_2050_v1 = [x + y for x, y in zip(rech_ies, delta_rech_gauss)]

rech_2050_v1_mean = np.mean(rech_2050_v1)
rech_2050_v1_stev = np.std(rech_2050_v1, ddof = 1)


#%%------------------------------- FIGURE PAPER -------------------------------

ens_list = [rech_ies, rech_2050_v1]
color_dic = dict(zip(list(range(len(ens_list))), ['blue', 'red']))

fig, ax = plt.subplots(figsize = cm2inch(11.5, 7))

# Boxplot
plt.boxplot(ens_list, widths = 0.4, showbox = True, showcaps = False, 
            showmeans = True, showfliers = False, boxprops = dic_props, 
            whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
            meanprops = meanprops, flierprops = flierprops)

# Violin plot (white)
vp = ax.violinplot(ens_list, widths = 0.8, points = 100, showmeans = False, 
                   showmedians = False, showextrema = False)
for pc in vp['bodies']:
    pc.set_facecolor('white')
    pc.set_edgecolor('white')
    pc.set_alpha(1)

# Violin plot (color)
vp = ax.violinplot(ens_list, widths = 0.8, points = 100, showmeans = False,
                   showmedians = False, showextrema = False)
for i, pc in enumerate(vp['bodies']):
    pc.set_facecolor(color_dic[i])
    pc.set_edgecolor(color_dic[i])
    pc.set_alpha(0.4)

# Show points
for i, tick in enumerate([1, 2]):
    y = ens_list[i]
    x = np.random.normal(tick, 0.04, size = len(y))
    plt.plot(x, y, 'k.', alpha = 0.1)

plt.xticks([1, 2], ['Current', '2050'], rotation = 0)
ax.set_axisbelow(True)
ax.yaxis.grid()
plt.yticks(np.arange(100, 1000, 100))
ax.set_ylim(150, 1000)
plt.ylabel('Recharge (mm/year)')
plt.tight_layout()
plt.savefig(os.path.join(path_postproc_out, '00paper_rech_2050_vs_hist.png'), dpi = 300)


#%%------------------ EXPORT 2050 SLVL & RECHARGE ENSEMBLES -------------------

# Export ensembles to csv

save_ens_file = 'rech_slvl_ensembles.csv'
pars = parens[['slvl', 'rech_mmyr']].copy()

pars.loc[:, 'slvl'] = slvl_2050
pars.drop(columns = ['rech_mmyr'], inplace = True)
pars['rech_mmyr_v1'] = rech_2050_v1
pars['rech_mmyr_v2'] = rech_2050_v2

pars.to_csv(os.path.join(path_postproc_out, save_ens_file), index = True)

pars['delta_rech_gauss'] = delta_rech_gauss
pars.to_csv(os.path.join(path_postproc_out, save_ens_file.replace('.csv', '_details.csv')), index = True)


#%%---------------------------- MODIFY THE PARENS -----------------------------

rech_slvl_ensembles = pd.read_csv(os.path.join(path_postproc_out, save_ens_file),
                                 index_col = 0)
parens.loc[:, 'rech_mmyr'] = rech_slvl_ensembles.loc[:, 'rech_mmyr_v1']
parens.loc[:, 'slvl'] = rech_slvl_ensembles.loc[:, 'slvl']

# Export modified parens

parens_file_modif = parens_file.replace('par.csv', 'par.cc.csv')
parens.to_csv(os.path.join(path_postproc_out, parens_file_modif), index = True)
