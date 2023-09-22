"""----------------------------------------------
Antoine Hermant, Sep 2023
----------------------------------------------"""

import xarray as xr
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.colors as mcolors

dataset_dir = '../datasets/'
figure_dir = '../figures/'

colors={'direct':'#3d7f1e',
        'indirect':'#b1116d',
        'total':'darkorange'
        }

def WeightedMean(var):
    weights = np.cos(np.deg2rad(var.lat))
    if len(var.shape) >= 2:
        return var.weighted(weights).mean(dim=('lat', 'lon'))
    else:
        return var.weighted(weights).mean("lat")

def individual_emissions():
    return pd.DataFrame({
    'Europe': [8.95, 18.26, 56.80, 16.41, 11.35],
    'North America': [7.65, 24.09, 29.32, 17.45, 7.39],
    'East Asia': [0.17, 1.69, 14.72, 37.36, 34.89],
    'South Asia': [0.18, 1.57, 9.18, 17.17, 22.89],
    'North Africa': [0.08, 0.20, 1.02, 1.70, 1.94],
    'South America': [0.12, 1.24, 4.81, 4.88, 5.26],
    'Maritime Continent': [0.03, 0.22, 2.13, 4.15, 4.43],
    'South Central Africa': [0.04, 1.42, 3.83, 3.35, 4.29],
    'Australia': [0.23, 0.63, 1.61, 1.57, 1.39],
    #'Global': [17.54, 49.32, 123.42, 104.04, 70.94]
    },
    index=[1900, 1950, 1980, 2005, 2013])

def global_emissions():
    individualEmissions = individual_emissions()
    MACSP = xr.open_dataset(f'{dataset_dir}/MAC-SP.nc')
    globalEmissions = MACSP.year_weight.sel(plume_number=1).isel(years=slice(0,164)).values * 0
    for year in range(0, 164, 1):
        for plume_number in range(1, 10 ,1):
            globalEmissions[year] = globalEmissions[year] + MACSP.year_weight.sel(plume_number=plume_number).isel(years=year)*individualEmissions.iloc[:, plume_number-1].loc[2005]

    return globalEmissions

def plume_locs(ax):
        macsp = xr.open_dataset(f'{dataset_dir}/MAC-SP.nc')
        transform = ccrs.PlateCarree()._as_mpl_transform(ax)
        
        for iplume in macsp.plume_number:
            ax.plot(
                macsp.plume_lon.sel(plume_number=iplume),
                macsp.plume_lat.sel(plume_number=iplume),
                marker="^",
                color="darkorange",
                markersize=3,
                transform=transform,
            )

def norm():
    class MidpointNormalize(mcolors.Normalize):
        def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
            self.midpoint = midpoint
            mcolors.Normalize.__init__(self, vmin, vmax, clip)

        def __call__(self, value, clip=None):
            v_ext = np.max( [ np.abs(self.vmin), np.abs(self.vmax) ] )
            x, y = [-v_ext, self.midpoint, v_ext], [0, 0.5, 1]
            return np.ma.masked_array(np.interp(value, x, y))
        
    return MidpointNormalize( midpoint = 0 )