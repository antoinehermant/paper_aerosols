"""----------------------------------------------
Antoine Hermant, Sep 2023
----------------------------------------------"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy.util import add_cyclic_point
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import ticker
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

from functions import WeightedMean, plume_locs, dataset_dir, figure_dir

def effect_map(ax, effect, type, year, sub):
    dataset = xr.open_dataset(f'{dataset_dir}/historical_simple-plumes_yearmean.nc').squeeze().isel(time=year-1850)

    if type == 'all-sky':
        r = 'rad0'
    elif type == 'clear-sky':
        r = 'raf0'
    
    if effect == 'total':
        forcing = dataset[f'dR_sp_s{r}'] + dataset[f'dR_sp_t{r}']
        title = 'Total Aerosol Effect (all-sky)'
    elif effect == 'direct':
        forcing = dataset[f'dR_spd_s{r}'] + dataset[f'dR_spd_t{r}']
        title = 'Aerosol Direct Effect (all-sky)'
    elif effect == 'indirect':
        forcing = dataset[f'dR_spi_s{r}'] + dataset[f'dR_spi_t{r}']
        title = 'Aerosol Indirect Effect (all-sky)'
    
    if type == 'clear-sky':
        r = 'raf0'
        title = 'Clear-sky Aerosol Effect'
    
    var = forcing
    lon = var.coords['lon']
    lon_idx = var.dims.index('lon')
    wrap_var, wrap_lon = add_cyclic_point(var.values, coord=lon, axis=lon_idx)

    pc = ax.contourf(wrap_lon,
        var.lat,
        wrap_var,
        transform=ccrs.PlateCarree(),
        cmap='RdBu_r',
        levels= np.linspace(-3, 3, 16),
        extend="both",
    )
    axins = inset_axes(ax, width='100%', height='10%', loc="lower center", borderpad=-2)
    cb = plt.colorbar(pc, cax=axins, ax=ax, orientation='horizontal', pad=0.04)
    cb.set_label(label=r'Radiative Forcing [Wm$^{-1}$]', size=10)

    ax.coastlines()
    ax.annotate(r'$\overline{dR}$ = '+'${:.3f}$'.format(WeightedMean(var)),
                [0.05, 0.3], xycoords="axes fraction", fontsize=10)
    ax.annotate(sub, [-0.05, 1.05], xycoords="axes fraction", fontsize=13)
    ax.set_title(title, fontsize=13)   

    return 0 

def cldcov(ax, year, sub):

    dataset = xr.open_dataset(f'{dataset_dir}/historical_simple-plumes_yearmean.nc').squeeze().isel(time=year-1850)
    var = dataset.aclcov
    lon = var.coords['lon']
    lon_idx = var.dims.index('lon')
    wrap_var, wrap_lon = add_cyclic_point(var.values, coord=lon, axis=lon_idx)

    pc = ax.contourf(wrap_lon,
        var.lat,
        wrap_var,
        transform=ccrs.PlateCarree(),
        cmap='YlGnBu',
        #cmap='PuBu',
        levels= np.linspace(0, 1, 11)
    )
    ax.coastlines()
    axins = inset_axes(ax, width='100%', height='10%', loc="lower center", borderpad=-2)
    cb = plt.colorbar(pc, cax=axins, ax=ax, orientation='horizontal', pad=0.04)
    cb.set_label(label='Fraction', size=10)
    
    ax.set_title('Cloud Cover Fraction', fontsize=13)
    ax.annotate(sub, [-0.05, 1.05], xycoords="axes fraction", fontsize=13)

    return 0

def aod(ax, year, sub):
    dataset = xr.open_dataset(f'{dataset_dir}/historical_simple-plumes_yearmean.nc').squeeze().isel(time=year-1850)
    var = dataset.aod_sp
    lon = var.coords['lon']
    lon_idx = var.dims.index('lon')
    wrap_var, wrap_lon = add_cyclic_point(var.values, coord=lon, axis=lon_idx)
    formatter = ticker.FormatStrFormatter('%.2f')

    level = [0, 0.01, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.38]
    pc = ax.contourf(wrap_lon,
        var.lat,
        wrap_var,
        transform=ccrs.PlateCarree(),
        cmap='Purples',
        levels= level,
        extend='max',
    )
    C = ax.contour(wrap_lon,
        var.lat,
        wrap_var,
        transform=ccrs.PlateCarree(),
        colors='black',
        levels= level,
        linewidths=.4,
        extend='max'
    )
    cont = ax.contour(wrap_lon,
        var.lat,
        wrap_var,
        transform=ccrs.PlateCarree(),
        levels=[0.0025], colors='black', linewidths=0.4, linestyles='--')
    ax.coastlines()

    axins = inset_axes(ax, width='100%', height='10%', loc='lower center', borderpad=-2)
    cb = plt.colorbar(pc, cax=axins, ax=ax, orientation='horizontal', pad=0.04, ticks=level)
    cb2 = plt.colorbar(C, cax=axins, ax=ax, orientation='horizontal', pad=0.04, ticks=level)
    cb.add_lines(cont)
    cb.lines[-1].set_linestyles(cont.linestyles)

    def custom_formatter(x, pos):
        if x == 0.002:
            return '0.002%'
        elif x == 0.00:
            return '0%'
        return formatter(x, pos)
    cb.ax.xaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))

    ax.annotate(r'$\overline{AOD}$ = '+'${:.3f}$'.format(WeightedMean(var)),
            [0.05, 0.3], xycoords="axes fraction", fontsize=10)
    
    cb.set_label(label='Aerosol Optical Depth', size=10)
    ax.set_title('Column Simple-Plumes AOD at 550nm', fontsize=13)
    ax.annotate(sub, [-0.05, 1.05], xycoords="axes fraction", fontsize=13)

    return 0

def plot_figure():
    year = 2005
    central_longitude = 0
    projectionStyle = ccrs.Robinson(central_longitude=central_longitude)
    fig, ax = plt.subplots(3,2, subplot_kw=dict(projection=projectionStyle))

    effect_map(ax[0,0], 'total', 'all-sky', year, 'a)')
    effect_map(ax[0,1], 'indirect', 'all-sky', year, 'b)')
    effect_map(ax[1,0], 'direct', 'all-sky', year, 'c)')
    effect_map(ax[1,1], 'direct', 'clear-sky', year, 'd)')
    cldcov(ax[2,0], year, 'e)')
    aod(ax[2,1], year, 'f)')

    for ax in [
        ax[0,0],
        ax[0,1],
        ax[1,0],
        ax[1,1],
        ax[2,1],
        ]:
        plume_locs(ax)

    fig = plt.gcf()
    fig.set_size_inches(10, 11)
    plt.subplots_adjust(wspace=0.1)

    plt.savefig(f'{figure_dir}/figure3.pdf', bbox_inches='tight')

    return 0

def main():
     plot_figure()
    
if __name__ == "__main__":
    main()

