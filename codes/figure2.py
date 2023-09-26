"""----------------------------------------------
Antoine Hermant, Sep 2023
----------------------------------------------"""

import xarray as xr
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import ticker
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
from sklearn.linear_model import LinearRegression

from functions import individual_emissions, dataset_dir, figure_dir

def plot_figure():
    dataset = xr.open_dataset(f'{dataset_dir}/historical_simple-plumes_fldmean_yearmean.nc').squeeze()
    MACSP = xr.open_dataset(f'{dataset_dir}/MAC-SP.nc')
    individualEmissions = individual_emissions()
    fig = plt.figure()
    gs = GridSpec(2,2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax3 = fig.add_subplot(gs[1, 0], sharex=ax1)
 
    subplots = {'direct': ax1,
                 'indirect': ax2,
                 'clear-sky':ax3}

    effects = ['direct', 'indirect','clear-sky']

    for ax in [ax1, ax2, ax3]:
        ax.spines['bottom'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_position('zero')
        for k, spine in ax.spines.items():
            spine.set_zorder(20)
        ax.set_xlim([-0.5, 51])
        ax.set_xlabel(r"Emissions [Tg SO$_2$]")
        ax.tick_params(axis='x', direction='out', which='major', pad=-15)


    coefs = pd.DataFrame({'direct': [0,0,0,0,0,0,0,0,0],
                  'clear-sky': [0,0,0,0,0,0,0,0,0]},
                   index=['Europe', 'North America', 'East Asia', 'South Asia', 'North Africa', 'South America', 'Maritime Continent', 'South Central Africa', 'Australia'])
    ms=4
    zorder = [1, 4, 3, 2, 9, 8, 7, 6, 5]
    for effect in effects:
        for plume_number in range(1, 10, 1):
            label = MACSP.attrs[f'plume{plume_number}_region']
            emissions = MACSP.year_weight.sel(plume_number=plume_number).isel(years=slice(0,164))*individualEmissions.iloc[:, plume_number-1].loc[2005]

            if effect == 'total':
                forcing = dataset[f'dR_sp{plume_number}_srad0'] + dataset[f'dR_sp{plume_number}_trad0']
                subplots[effect].scatter(emissions, forcing, label=label, s=ms, zorder=zorder[plume_number-1])
            elif effect == 'direct':
                forcing = dataset[f'dR_spd{plume_number}_srad0'] + dataset[f'dR_spd{plume_number}_trad0']
                subplots[effect].scatter(emissions, forcing, label=label, s=ms, zorder=zorder[plume_number-1], marker='x', linewidths=0.8)
                X = emissions.values.reshape(-1, 1)
                reg = LinearRegression().fit(X, forcing)
                subplots[effect].plot(X, X*reg.coef_ + reg.intercept_, linewidth=1, linestyle='-')
                coefs['direct'][label] = reg.coef_

            elif effect == 'indirect':
                forcing = dataset[f'dR_spi{plume_number}_srad0'] + dataset[f'dR_spi{plume_number}_trad0']
                subplots[effect].scatter(emissions, forcing, label=label, s=ms*0.6, zorder=zorder[plume_number-1])
            
            elif effect == 'clear-sky':
                forcing = dataset[f'dR_spd{plume_number}_sraf0'] + dataset[f'dR_spd{plume_number}_traf0']
                subplots[effect].scatter(emissions, forcing, label=label, s=ms, zorder=zorder[plume_number-1], marker='x', linewidths=0.8)
                X = emissions.values.reshape(-1, 1)
                reg = LinearRegression().fit(X, forcing)
                subplots[effect].plot(X, X*reg.coef_ + reg.intercept_, linewidth=1, linestyle='-')
                coefs['clear-sky'][label] = reg.coef_

    ax1.annotate(str(round(coefs['direct']['Europe']*1000,2))+r'$\cdot10^{-3}$', [0.86, 0.82],xycoords="axes fraction", fontsize=8)
    ax1.annotate(str(round(coefs['direct']['North America']*1000,2))+r'$\cdot10^{-3}$', [0.66, 0.72],xycoords="axes fraction", fontsize=8)
    ax1.annotate(str(round(coefs['direct']['East Asia']*1000,2))+r'$\cdot10^{-3}$', [0.78, 0.5],xycoords="axes fraction", fontsize=8)
    ax1.annotate(str(round(coefs['direct']['South Asia']*1000,2))+r'$\cdot10^{-3}$', [0.45, 0.02],xycoords="axes fraction", fontsize=8)

    ax3.annotate(str(round(coefs['clear-sky']['Europe']*1000,2))+r'$\cdot10^{-3}$', [0.97, -0.02],xycoords="axes fraction", fontsize=8)
    ax3.annotate(str(round(coefs['clear-sky']['North America']*1000,2))+r'$\cdot10^{-3}$', [0.66, 0.55],xycoords="axes fraction", fontsize=8)
    ax3.annotate(str(round(coefs['clear-sky']['East Asia']*1000,2))+r'$\cdot10^{-3}$', [0.78, 0.3],xycoords="axes fraction", fontsize=8)
    ax3.annotate(str(round(coefs['clear-sky']['South Asia']*1000,2))+r'$\cdot10^{-3}$', [0.44, 0.02],xycoords="axes fraction", fontsize=8)
    
    ax1.set_title('Direct effect', pad=20)
    ax1.annotate('a)', [-0.15, 1.1],xycoords="axes fraction", fontsize=13)
    ax1.xaxis.set_label_coords(.82, .97)
    ax1.set_ylim([-0.205, 0.003])
    ax1.set_yticks([-0.2, -.15, -0.1, -0.05, 0])
    ax1.set_ylabel("Radiative Forcing [Wm$^{-2}$]")  
    ax1.tick_params(axis='y', direction='in', which='both', pad=5)
   
    ax2.annotate('b)', [-0.05, 1.1],xycoords="axes fraction", fontsize=13)
    ax2.set_title('Indirect effect', pad=20)
    ax2.xaxis.set_label_coords(.82, .95)
    ax2.set_ylim([-0.205, 0.003])
    ax2.tick_params(axis='y', labelleft=False, which='both', direction='in')

    ax3.annotate('c)', [-0.15, 1.1],xycoords="axes fraction", fontsize=13)
    ax3.set_title('Clear-sky', pad=20)
    ax3.xaxis.set_label_coords(.82, .95)
    ax3.set_ylim([-0.36, 0.005])
    ax3.set_yticks([-0.3, -.2, -0.1, 0])
    ax3.set_ylabel("Radiative Forcing [Wm$^{-2}$]")  
    ax3.tick_params(axis='y', direction='in', which='both', pad=5)

    leg = ax2.legend(frameon=True, ncol=1,
                handletextpad=0.5, handlelength=1.0, loc=(.2,-1.2), markerscale=2)
    leg.set_zorder(1)
    fig = plt.gcf()
    fig.set_size_inches(9, 7.5)
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.1)

    plt.savefig(f'{figure_dir}/figure2.pdf', bbox_inches='tight')

    return 0

def main():
    plot_figure()

if __name__ == '__main__':
    main()