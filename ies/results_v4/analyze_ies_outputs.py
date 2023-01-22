# Author: Cecile Coulon

import os
working_dir = './'
os.chdir(working_dir)
modflow_path = 'mf2005'

if not os.path.isdir('figures'):
        os.mkdir('figures')

# Load python libraries

import pyemu, flopy, pickle, copy, collections
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

name_pst = 'pestpp_ies'
krig_factors_file = os.path.join('..', '..', 'send_to_ies_v4', 'model_params', 'pp_30max.fac') # Kriging factors
pp_df = pyemu.pp_utils.pp_file_to_dataframe(os.path.join('..', '..', 'send_to_ies_v4', 'model_params', 'hk1pp.dat')) # Pilot point file
pp_df.index = pp_df.name


#%%-------------------------- Read PEST control file --------------------------

pst = pyemu.Pst(os.path.join('./', name_pst + '.pst'))

obs = pst.observation_data
par = pst.parameter_data

fcstnames = [s for s in pst.obs_names if s.startswith('zmuni') and not s.endswith('pct')]
obsnmes = obs.loc[obs.weight != 0, 'obsnme'].tolist()

pp_names_adj = par[(par.parnme.str.startswith('hk1')) & (par.partrans != 'tied') ].index.tolist()#[s for s in pst.adj_par_names if s.startswith('hk1')]
pp_names_adj2 = [s.replace('hk1', 'pp_0') for s in pp_names_adj]

pp_names_tied = par[(par.parnme.str.startswith('hk1')) & (par.partrans == 'tied') ].index.tolist()
pp_names_tied2 = [s.replace('hk1', 'pp_0') for s in pp_names_tied]

pp_names_subset = [s for s in pst.adj_par_names if s.startswith('hk1') and s != 'hk1000']
par_names_others = [s for s in pst.adj_par_names if s not in pp_names_subset]

# Right-hand side of the constraint
c = obs.loc[obs.obsnme.apply(lambda x: x.endswith('pct')), 'obsnme']
cons_rhs = obs.loc[c, 'obsval'].to_frame()


#--- Read phi

phi = pd.read_csv(os.path.join('./', name_pst + '.phi.actual.csv'), index_col = 0)
iterations = list(range(phi.index[-1] + 1))


#%%------------------------------ Read ensembles ------------------------------

#--- Read parameter ensembles

parens_dic = {}
parens_reals_dic = collections.defaultdict(dict)
for i in iterations:
    parens_dic[i] = pd.read_csv(os.path.join('./', name_pst + '.{0}.par.csv'.format(i)), index_col = 0)
    for real in parens_dic[i].index:
        parens_reals_dic[real][i] = parens_dic[i].loc[real, :]

for real in parens_dic[i].index:
    df = pd.DataFrame(parens_reals_dic[real])
    df['parnme'] = df.index
    parens_reals_dic[real] = df


#--- Read observation ensembles

obsens_dic = {}
for i in iterations:
    obsens_dic[i] = pd.read_csv(os.path.join('./', name_pst + '.{0}.obs.csv'.format(i)), index_col = 0)

obsplusnoise_ens = pd.read_csv(os.path.join('./', 'pestpp_ies.obs+noise.csv'), index_col = 0)

#--- Compute parameter ensemble statistics

parstats_dic = {}
parstats_pp_dic = {}
parstats_condensed_dic = {}

for i in iterations:
    
    df = parens_dic[i].loc[:, pst.adj_par_names]
    parstats_dic[i] = pd.DataFrame(
        {'mean': df.mean(), 'stdev': df.std(ddof = 1), 'pct5': df.quantile(0.05), 'pct95': df.quantile(0.95)})
    
    df = parens_dic[i].loc[:, par_names_others]
    parstats_condensed_dic[i] = pd.DataFrame(
        {'mean': df.mean(), 'stdev': df.std(ddof = 1), 'pct5': df.quantile(0.05), 'pct95': df.quantile(0.95)})
    
    ll = parens_dic[i].loc[:, pp_names_subset].stack().tolist()
    parstats_pp_dic[i] = pd.DataFrame(
        {'mean': [np.mean(ll)], 'stdev': [np.std(ll, ddof = 1)], 
         'pct5': [np.quantile(ll, 0.05)],  'pct95': [np.quantile(ll, 0.95)]},
        index = ['hk_pp'])
    
    parstats_condensed_dic[i] = parstats_condensed_dic[i].append(parstats_pp_dic[i])

#--- Compute forecast ensemble statistics

fcststats_dic = {}
fcststats_dic2 = {}

mean_dic, stdev_dic, pct1_dic, pct5_dic, pct95_dic, pct99_dic = {}, {}, {}, {}, {}, {}
for i in iterations:
    fcststats_dic[i] = pd.DataFrame({'mean': obsens_dic[i].loc[:, fcstnames].mean(),
                                 'stdev': obsens_dic[i].loc[:, fcstnames].std(),
                                 'pct1': obsens_dic[i].loc[:, fcstnames].quantile(0.01),
                                 'pct5': obsens_dic[i].loc[:, fcstnames].quantile(0.05), 
                                 'pct95': obsens_dic[i].loc[:, fcstnames].quantile(0.95),
                                 'pct99': obsens_dic[i].loc[:, fcstnames].quantile(0.99)})
    
    mean_dic[i] = fcststats_dic[i]['mean']
    stdev_dic[i] = fcststats_dic[i]['stdev']
    pct1_dic[i] = fcststats_dic[i]['pct1']
    pct5_dic[i] = fcststats_dic[i]['pct5']
    pct95_dic[i] = fcststats_dic[i]['pct95']
    pct99_dic[i] = fcststats_dic[i]['pct99']

fcststats_dic2 = dict({'mean': pd.DataFrame(mean_dic).T,
                       'stdev': pd.DataFrame(stdev_dic).T,
                       'pct1': pd.DataFrame(pct1_dic).T,
                       'pct5': pd.DataFrame(pct5_dic).T, 
                       'pct95': pd.DataFrame(pct95_dic).T,
                       'pct99': pd.DataFrame(pct99_dic).T})


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

def fmt_label(x, pos):
    ''' Function later used to label the colorbar using scientic notation '''
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


#%%------------------------------- FIGURE PAPER -------------------------------

# Plot simulated vs observed

i = iterations[-1]
ensemble = obsens_dic[i]

for obs_type in ['h', 'z']:
    
    #--- Load observations
    
    obs = pst.observation_data
    grouper = obs.groupby(obs.obgnme).groups # Group observations
    skip_groups = ['zmuni', 'l_zmuni_pct', 'v_lens']
    
    #--- Define legends and colors
    
    if obs_type == 'h':
        skip_groups = skip_groups + ['zwells', 'ztdem', 'zert']
        grouper_lgd = {'hdeep': r'$h_{\rm f\ deep\ wells}$', 
                       'hshallow': r'$h_{\rm f\ shallow\ wells}$', 
                       'hmuni': r'$h_{\rm f\ pumping\ wells}$'}    
        grouper_clr = {'hdeep': 'red', 'hshallow': 'lime', 'hmuni': 'black'}
        title = 'Freshwater heads ' + r'$h_{\rm f}$'
        ylabel = 'Simulated (m)'
    
    elif obs_type == 'z':
        skip_groups = skip_groups + ['hdeep', 'hshallow', 'hmuni']
        grouper_lgd = {'zwells': r'$\zeta_{\rm\ deep\ wells}$', 
                       'ztdem': r'$\zeta_{\rm\ TDEM}$', 
                       'zert': r'$\zeta_{\rm\ ERT}$'}
        grouper_clr = {'zwells': 'cyan', 'ztdem': 'cornflowerblue','zert': 'navy'}
        title = 'Interface elevations ' + r'$\zeta$'
        ylabel = 'Residual (m)'
    
    for skip_group in skip_groups: # Drop the obs groups that are not to be plotted
        grouper.pop(skip_group)
    
    #--- Plot
    
    fig, ax = plt.subplots(figsize = cm2inch(9.5, 9.5))
    
    for groupname, obsnames in grouper.items(): # For each obs group
        
        obs_g = obs.loc[obsnames, :].sort_values(by = "obsval") # observation group
        obs_g = obs_g.loc[obs_g.weight > 0, :] # Exclude observations with weight = 0
        
        en_g = ensemble.loc[:, obs_g.obsnme] # simulated ensembles
        base_g = obsplusnoise_ens.loc[:, obs_g.obsnme] # observation + noise ensemble
        
        
        # Plot simulated ensembles
        
        linemin = en_g.min().values
        linemax = en_g.max().values
        
        [ax.plot([obsval, obsval], [emin, emax], label = grouper_lgd[groupname], 
                 color = grouper_clr[groupname], zorder = 1, alpha = 0.8)
         for obsval, emin, emax in zip(obs_g.obsval.values, linemin, linemax)]
        
        # Plot simulated mean
        ax.scatter(obs_g.obsval.values, en_g.mean().values, s = 5,
                   color = grouper_clr[groupname], alpha = 0.8)
        
    #--- Plot limits
    
    if obs_type == 'h':
        Max = 3
        Min = 0
    elif obs_type == 'z':
        Max = 0
        Min = -90
        
    ax.plot([Min,Max], [Min,Max], 'k--', lw = 1.0, zorder = 3) # Plot the 1:1 diagonal line
    ax.set_xlim(Min,Max)
    ax.set_ylim(Min,Max)
    loc = 'upper left'
    
    #--- Legend
    
    handles, labels = plt.gca().get_legend_handles_labels()
    labels_bis = list(np.unique(labels))
    idx = [labels.index(labels_bis[i]) for i in range(len(labels_bis))]
    handles_bis = [handles[i] for i in idx]
    order = [2,1,0]
    ax.legend([handles_bis[i] for i in order], [labels_bis[i] for i in order], loc = loc)
    
    #--- Other figure parameters
    
    ax.axis('square')
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel('Observed (m)', labelpad = 0.1)
    ax.set_ylabel(ylabel, labelpad = 0.1)
    ax.set_title(title, loc = 'left', fontsize = fs)
    
    plt.tight_layout()
    plt.savefig(os.path.join('figures', '00paper_1to1_' + obs_type + '_it' + str(i) + '.png'), dpi = 300)


#%%------------------------ Export ensemble statistics ------------------------

writer = pd.ExcelWriter(os.path.join('figures', 'par_stats.xlsx'), engine = 'xlsxwriter')
for i in iterations:
    parstats_dic[i].to_excel(writer, sheet_name = 'iteration' + str(i))
writer.save()

writer = pd.ExcelWriter(os.path.join('figures', 'par_stats_condensed.xlsx'), engine = 'xlsxwriter')
for i in iterations:
    parstats_condensed_dic[i].to_excel(writer, sheet_name = 'iteration' + str(i))
writer.save()


#%%------------------------------- FIGURE PAPER -------------------------------

# Plot K maps

# Retrieve data from the bin file generated with QGridder
path_bin_in = os.path.join('..', '..', 'send_to_ies_v3', 'preproc_IDM_20m_v6.bin')

with open(path_bin_in, 'rb') as handle:
     objects = pickle.load(handle)
     (nrow, ncol, delr, delc, ibound, sea, dem, geology, thk_valley, cx, cy, 
      muni_wells_row_col, muni_pumping, domestic_wells_row_col, ind_old_wells_row_col,
      obs_wells_row_col, obs_head, obs_zeta, tdem_row_col, ert_row_col) = objects

mf = flopy.modflow.Modflow.load(os.path.join('..', '..', 'simu_ies_swi.nam'))
botm = -300
thk_dunes = 10

# Plot K fields at iteration 2

i = 2 

for real in ['60', '162', '180', 'base']:
    
    pp_file = os.path.join('figures', 'hk1pp_' + real + '_it' + str(i) + '.dat')    
    
    #--- Read pp values
    
    # Adjustable pps
    parens_reals_dic2 = parens_reals_dic.copy()
    parens_reals_dic2[real]['parnme2'] = parens_reals_dic2[real].index.str.replace('hk1', 'pp_0')
    parens_reals_dic2[real].index = parens_reals_dic2[real].parnme2
    pp_df.loc[pp_names_adj2, 'parval1'] = parens_reals_dic2[real].loc[pp_names_adj2, i]
    
    # Tied pps
    pp_df.loc[pp_names_tied2, 'parval1'] = parens_reals_dic2[real].loc[pp_names_tied2, i]
    
    pyemu.utils.pp_utils.write_pp_file(pp_file, pp_df)
    
    # Create hk grid from hk values determined at pilot points ('hk1pp.dat' pp_file) & kriging factors ('pp.fac' factors_file)
    hk_EDC_array = pyemu.utils.fac2real(pp_file = pp_file, factors_file = krig_factors_file, 
                                        out_file = None, fill_value = np.nan)
    
    #--- Read hk_valley, hk_dunes values
    
    hk_valley = parens_reals_dic2[real].loc['hk_valley', i]
    hk_dunes = parens_reals_dic2[real].loc['hk_dunes', i]
    
    #--- Compute hk array
    hk_array = copy.deepcopy(hk_EDC_array) # Initialize
    hk_valleyH = ( hk_EDC_array * (dem - botm - thk_valley) + hk_valley * thk_valley ) / (dem - botm)
    hk_dunesH = ( hk_EDC_array * (dem - botm - thk_dunes) + hk_dunes * thk_dunes ) / (dem - botm)
    hk_array[geology == 3] = hk_valleyH[geology == 3] # Paleoglacial valley
    hk_array[geology == 4] = hk_dunesH[geology == 4] # Sand dunes
    
    #--- Plot hk_array
    
    fig = plt.figure(figsize = cm2inch(9, 6))
    ax = fig.add_subplot(1, 1, 1, aspect = 'equal')
    fig.subplots_adjust(wspace = 0.2, hspace = 0.2, left = 0.125, right = 1, 
                        bottom = 0, top = 1)
    
    mapview = flopy.plot.PlotMapView(model = mf)
    
    # Hydraulic conductivity grid
    im = mapview.plot_array(hk_array)
    # Coastline
    mapview.contour_array(sea, colors = 'k', levels = np.array([0.5]), linewidths = 0.5)
    
    #--- Figure parameters
    
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))  
    plt.tight_layout()
    plt.savefig(os.path.join('figures', '00paper_real' + real + '_Kmap.png'), dpi = 300)


parens_reals_dic['base'].loc['rech_mmyr', 2]
parens_reals_dic['60'].loc['rech_mmyr', 2]
parens_reals_dic['162'].loc['rech_mmyr', 2]
parens_reals_dic['180'].loc['rech_mmyr', 2]
