#! /home/xychen/env/default/bin/python3

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#from mpl_toolkits import mplot3d
#from matplotlib.ticker import LinearLocator
from matplotlib import cm
import os


netCDF_FN=dict()
netCDF_FN['contrl']='tmp_dataRICO_128x128x120_dx100m_standard_test_RRTM4PBL_radON_higher_sounding_quartile_mean_ds_for_testing_updated_3days.nc'
netCDF_FN['large']='tmp_dataRICO_512x512x120_dx250m_largedomain_test_RRTM4PBL_radON_realistic_higher_sounding_quartile_mean_ds.nc'

# read both:
ds_rad=dict()
ds_rad['contrl']= xr.open_dataset(netCDF_FN['contrl'])
ds_rad['large'] = xr.open_dataset(netCDF_FN['large'])


# In[7]:



# plot variable names direclty. 
dt = (ds_rad['contrl'].time[1]-ds_rad['contrl'].time[0])*24
nt = int(2/dt)
tsm_window=nt
varname='PW'
def movmean(x, w):
   return np.convolve(x, np.ones(w), 'valid')/w

colors=plt.cm.RdBu(np.linspace(0,1, 4))
fig, ax = plt.subplots(1,2, figsize = (12,6))
for k, key in enumerate(ds_rad):
   t = ds_rad[key].time - ds_rad[key].time[0]
   quartiles = ds_rad[key].quartile.values


   for i, qrt in enumerate(quartiles):
       var = ds_rad[key].sel(quartile=qrt)[varname]

       #if smoothing is requested for the time evolution:

       t = movmean(ds_rad[key].time - ds_rad[key].time[0], tsm_window)
       var = movmean(var, tsm_window)

       ax[k].plot(t, var, '-', linewidth=1.2, color=colors[i], label=qrt)


   ax[k].set_xlabel(ds_rad[key].time.long_name + ' (' + ds_rad[key].time.units +')', fontsize=14 )
#   ax[k].set_ylabel(ds_2D[varname].long_name + ' ('+ ds_2D[varname].units +')' +'\n quartile mean', fontsize=14)
   ax[k].set_ylim(35.5, 39.5)
   ax[k].set_xlim(0, 1.5)
   ax[k].legend()
   ax[k].grid(True)
   ax[k].set_title('rad ' + key,fontsize=18)
   ax[k].tick_params(axis='x',labelsize=14)
   ax[k].tick_params(axis='y',labelsize=14)

fig.tight_layout()
# save the two figures:
abs_svpath='/home/xychen/GitHub_repo/SAM_LES/Figs/demo'
figname='PW_quartile_1day_evolution_contrl_vs_largedomain_updated.jpg'
abs_figname = os.path.join(abs_svpath, figname)
fig.savefig(abs_figname)


