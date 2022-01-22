# -*- coding: utf-8 -*-
"""
Spyder Editor

This script is used to compare the results of two simulations;
Date: Jan 21, 2022
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import cm
import os

# load in the two datasets for comparison
datadir = '/work/noaa/ome/xychen/SAM_UW_RICO_test/OUT_2D'
netCDF_FN=dict()
netCDF_FN['MPDATA']='RICO_512x512x126_dx250m_largedomain6km_ADV_MPDATA_RRTM4PBL_perpetualsun_32x32gc_256.2Dbin_1.nc'
netCDF_FN['SELPPM']='RICO_512x512x126_dx250m_largedomain6km_test_RRTM4PBL_radon_perpetualsun_32x32gc_256.2Dbin_1.nc'

abs_svpath='/home/xychen/GitHub_repo/SAM_LES/Figs/local/RadOn_Perpetual_Top6km_2days'
figname='CWP_spatial_evolution_comparison_MPDATA_SELPPM.jpg'

# read both:
ds_2D=dict()
for key in netCDF_FN:
    ds_2D[key] = xr.open_dataset(os.path.join(datadir, netCDF_FN[key]))
    

simThr = [8, 16, 24, 32]
# ready to make some plots:
var = 'CWP'
nrow = 2
ncol = 4

colormap=cm.Blues


fig, axes = plt.subplots(nrow, ncol, figsize=(15,6))
# row: different datasets;column: different time in the simulations;
for i, key in enumerate(ds_2D):
    val = np.log10(ds_2D[key][var])
    time = (ds_2D[key]['time']-ds_2D[key]['time'][0]) * 24.  # hours
    x = ds_2D[key]['x']
    y = ds_2D[key]['y']
    [xx, yy] = np.meshgrid(x,y)
    nc = ds_2D[key]
    
    for t, ax in zip (simThr, axes[i].flatten()):
        # find the index for the simulation hours:
        #it = np.where(time==t)[0][0]
        timedif = np.abs(time-t)
        it = np.where(timedif == np.min(timedif))[0][0]
        
        labelstr = '{0:.0f}th hour'.format(t)
        hc = ax.pcolor(xx, yy, np.squeeze(val[it,:,:]), cmap=colormap, 
                         vmin=-9, vmax=1.5
                        )
        # add contours
        #cs = ax.contour(xx, yy, np.squeeze(val[it,:,:]), colors = 'k', linewidths=1)
    
        #ax.clabel(cs, cs.levels)
        
        
        if i==nrow-1:
            ax.set_xlabel(nc['x'].name + '(' + nc['x'].units + ')', fontsize=12)
        if t==simThr[0]:
            ax.set_ylabel(nc['y'].name + '(' + nc['y'].units + ')', fontsize=12)
        if i==0:
            ax.set_title(labelstr, fontsize=14)
        
        hb = plt.colorbar(hc, ax=ax)
        #plt.clim(0,12)
        hb.ax.set_title(val.units, fontsize=12)
        #if i%4==0:
        hb.ax.set_ylabel('log10(CWP)',fontsize=12)
        
        ax.tick_params(axis='x',labelsize=12)
        ax.tick_params(axis='y',labelsize=12)
        #ax.axis('square')
fig.tight_layout()

abs_figname = os.path.join(abs_svpath, figname)
fig.savefig(abs_figname)

