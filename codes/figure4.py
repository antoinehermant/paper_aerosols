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

from functions import WeightedMean, individual_emissions, dataset_dir, figure_dir

def plot_figure():

    dataset = xr.open_dataset(f'{dataset_dir}/historical_simple-plumes_yearmean.nc').squeeze()
    dataset_ssa = xr.open_dataset(f'{dataset_dir}/historical_simple-plumes_ssa093_yearmean.nc').squeeze()
    MACSP = xr.open_dataset(f'{dataset_dir}/MAC-SP.nc')
    individualEmissions = individual_emissions()

    fig = plt.figure()
    gs = GridSpec(2,2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    #ax3 = fig.add_subplot(gs[0, 1], sharey=ax2, sharex=ax1)

    effects = ['clear-sky', 'clear-sky-alb','ssa']

    for ax in [ax1, ax2]:
        ax.spines['bottom'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_position('zero')
        for k, spine in ax.spines.items():
            spine.set_zorder(20)
        ax.set_xlim([-0.0001, 0.0082])
        ax.xaxis.set_label_coords(.55, 1.13)
        ax.set_ylim([-0.38, 0.008])
        ax.tick_params(axis='x', direction='out', which='major', pad=-16)

    coefs = pd.DataFrame({'clear-sky': [0,0,0,0,0,0,0,0,0],
                  'clear-sky-alb': [0,0,0,0,0,0,0,0,0]},
                   index=['Europe', 'North America', 'East Asia', 'South Asia', 'North Africa', 'South America', 'Maritime Continent', 'South Central Africa', 'Australia'])
    ms=4
    zorder = [1, 4, 3, 2, 7, 9, 8, 6, 5]
    for effect in effects:
        for plume_number in range(1, 10, 1):
            label = MACSP.attrs[f'plume{plume_number}_region']
            clearsky = dataset[f'dR_spd{plume_number}_sraf0'] + dataset[f'dR_spd{plume_number}_traf0']
            aod = dataset[f'aod_sp{plume_number}']
            #albedo = 1 - (dataset.srafs)/(dataset.srad0d)  #clear-sky albedo
            #albedo = 1 / (1 - (- dataset.sradsu /(dataset.srads + dataset.sradsu)))
            #albedo_ssa = 1 - (dataset_ssa.srafs)/(dataset_ssa.srad0d)
            albedo_ssa = - dataset_ssa.sradsu /(dataset_ssa.srads - dataset_ssa.sradsu)
            clearsky_ssa = dataset_ssa[f'dR_spd{plume_number}_sraf0'] + dataset_ssa[f'dR_spd{plume_number}_traf0']

            if effect == 'clear-sky':
                x = WeightedMean(aod)
                X = x.values.reshape(-1, 1)
                y = WeightedMean(clearsky)
                reg = LinearRegression().fit(X, y)
                ax1.scatter(x, y, label=label, s=ms, zorder=zorder[plume_number-1], marker='x', linewidths=0.8)
                ax1.plot(X, x*reg.coef_ + reg.intercept_, linewidth=1, linestyle='-')
                coefs['clear-sky'][label] = reg.coef_

            #if effect == 'clear-sky-alb':
            #    x = WeightedMean(aod)
            #    X = x.values.reshape(-1, 1)
            #    y = WeightedMean(clearsky*(albedo))
            #    reg = LinearRegression().fit(X, y)
            #    ax2.scatter(x, y, label=label, s=ms, zorder=zorder[plume_number-1], marker='x', linewidths=0.8)
            #    ax2.plot(X, X*reg.coef_ + reg.intercept_, linewidth=1, linestyle='-')
            #    coefs['clear-sky-alb'][label] = reg.coef_

            if effect == 'ssa':
                ax2.scatter(WeightedMean(dataset_ssa[f'aod_sp{plume_number}']), WeightedMean(clearsky_ssa), label=label, s=ms*0.6, zorder=zorder[plume_number-1])

    ax1.annotate(str(round(coefs['clear-sky']['Europe'],3)), [0.85, 0.02],xycoords="axes fraction", fontsize=10)
    ax1.annotate(str(round(coefs['clear-sky']['North America'],3)), [0.12, 0.55],xycoords="axes fraction", fontsize=10)
    ax1.annotate(str(round(coefs['clear-sky']['East Asia'],3)), [0.37, 0.32],xycoords="axes fraction", fontsize=10)
    ax1.annotate(str(round(coefs['clear-sky']['South Asia'],3)), [0.58, 0.06],xycoords="axes fraction", fontsize=10)

    #ax2.annotate(str(round(coefs['clear-sky-alb']['Europe'],3)), [0.88, 0.0],xycoords="axes fraction", fontsize=10)
    #ax2.annotate(str(round(coefs['clear-sky-alb']['North America'],3)), [0.12, 0.6],xycoords="axes fraction", fontsize=10)
    #ax2.annotate(str(round(coefs['clear-sky-alb']['East Asia'],3)), [0.41, 0.39],xycoords="axes fraction", fontsize=10)
    #ax2.annotate(str(round(coefs['clear-sky-alb']['South Asia'],3)), [0.6, 0.2],xycoords="axes fraction", fontsize=10)

    ax1.annotate('a) Clear-sky Aerosol Forcing', [-0.05, 1.2],xycoords="axes fraction", fontsize=13)
    #ax1.set_title('Clear-sky', pad=20)
    ax1.set_xlabel(r"Global Mean AOD")
    ax1.set_yticks([-0.3, -0.2, -0.1, 0])
    ax1.set_ylabel("Radiative Forcing [Wm$^{-2}$]")  
    ax1.tick_params(axis='y', direction='in', which='both', pad=5)

    formatter_x = ticker.FormatStrFormatter('%.3f')
    formatter_y = ticker.FormatStrFormatter('%.2f')
    
    def custom_formatter(x, pos):
        if x == 0.000:
            return '0.0%'
        return formatter_x(x, pos)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))

    def custom_formatter(y, pos):
        if y == 0.00:
            return '0.0%'
        return formatter_y(y, pos)
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))
    
    ax2.annotate('b) Clear-sky with Equal SSA', [-0.05, 1.2],xycoords="axes fraction", fontsize=13)
    #ax2.set_title('Clear-sky with Surface Albedo', pad=20)
    ax2.set_xlabel(r"Global Mean AOD")
    ax2.tick_params(axis='y', direction='in', which='both', labelleft=False)

    leg = ax2.legend(frameon=True, ncol=1,
                handletextpad=0.5, handlelength=1.0, loc=(-0.35,-0.8), markerscale=1.8)
    leg.set_zorder(1)
    fig = plt.gcf()
    fig.set_size_inches(9, 7.5)
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.1)

    plt.savefig(f'{figure_dir}/figure4.pdf', bbox_inches='tight')

    return 0

def main():
    plot_figure()

if __name__ == '__main__':
    main()