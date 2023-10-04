"""----------------------------------------------
Antoine Hermant, Sep 2023
----------------------------------------------"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

from functions import global_emissions, colors, dataset_dir, figure_dir

def plot():
    dataset = xr.open_dataset(f'{dataset_dir}/historical_simple-plumes_fldmean_yearmean.nc').squeeze()
    globalEmissions = global_emissions()
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    time = np.arange(1850, 1850 + dataset.time.size, 1)
    lw = 1.6
    ax2.plot(time, (dataset.dR_spd_srad0.values + dataset.dR_spd_trad0.values), color=colors['direct'], linestyle='-', linewidth=lw, label='direct effect')
    ax2.plot(time, (dataset.dR_spi_srad0.values + dataset.dR_spi_trad0.values), color=colors['indirect'], linestyle='-', linewidth=lw, label='indirect effect')
    ax2.plot(time, (dataset.dR_sp_srad0.values + dataset.dR_sp_trad0), 'darkorange', linestyle='-', linewidth=lw, label='total effect')
    ax2.plot(time, (dataset.dR_sp_sraf0.values + dataset.dR_sp_traf0), 'darkblue', linestyle=':', linewidth=0.8, alpha=0.9, label='aerosol effect clear-sky')
    ax1.plot(time, globalEmissions, 'black', linewidth=lw, label='Global emissions')
    plt.axhline(0, color='black', linewidth=0.8)
    
    ax2.set_ylim([-.9, .9])    
    ax1.set_ylim([-max(globalEmissions), max(globalEmissions)+1])   

    fs=10.5
    ax1.set_xlabel('Year', fontsize=fs)
    ax1.set_xlim([1850, 2013])
    plt.xticks(np.arange(1850, 2001, 50))
    ax1.set_ylabel(r'Emissions [Tg SO$_2$-eq]', fontsize=fs, labelpad=30)
    ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.yaxis.set_label_coords(-.11,0.71)

    ax2.set_yticks([-0.8, -0.6, -0.4, -0.2, 0])
    ax2.set_ylabel(r'Radiative Forcing [Wm$^{-2}$]', rotation=-90, labelpad=70, fontsize=fs)
    ax2.yaxis.set_label_coords(1.2,0.26)

    ax1.tick_params(axis='both', which='major', direction='in', labelsize=fs,  pad=7)
    ax2.tick_params(axis='both', which='major', direction='in', labelsize=fs)
    
    legend_elements = [
    Line2D([0],[2],color="black",markerfacecolor="black",markersize=10,linewidth=1.1)]
    labels = ["Global emissions"]
    ax1.legend(handles=legend_elements, frameon=True, ncol=1, labels=labels, fontsize=fs,
    handletextpad=0.5, handlelength=1.0, loc=(0.05,0.85))

    legend_elements = [
    Line2D([0],[2],color=colors['direct'],markerfacecolor="black",markersize=10,linewidth=1.2),
    Line2D([0],[2],color=colors['indirect'],markerfacecolor="black",markersize=10,linewidth=1.2),
    Line2D([0],[2],color="darkorange",markerfacecolor="black",markersize=10,linewidth=1.2),
    Line2D([0],[2],color="darkblue",markerfacecolor="black",markersize=10,linewidth=0.8, alpha=0.9, linestyle=':'),
    ]
    labels = ["Direct effect", "Indirect effect", "Total effect", "Clear-sky"]
    leg = ax2.legend(handles=legend_elements, frameon=True, ncol=1, labels=labels, fontsize=fs,
    handletextpad=0.5, handlelength=1.0, loc=(0.05,0.05))

    ax1.axvline(x=1972, color='black', linewidth=0.3, linestyle='--')
    ax1.axvline(x=2005, color='black', linewidth=0.3, linestyle='--')
    plt.subplots_adjust(wspace=0.4)

    leg.set_zorder(1)
    fig = plt.gcf()
    fig.set_size_inches(4.5, 5)
    
    plt.tight_layout()
    plt.savefig(f'{figure_dir}/figure1.pdf',bbox_inches='tight')
    
    return 0

def main():
    plot()

if __name__ == '__main__':
    main()


