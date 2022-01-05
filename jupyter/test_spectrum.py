#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 10:45:35 2022

@author: xchen

This script is used to test wnspectrum.py functions.
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import metpy.calc as metcalc
from metpy.units import units

import os
import sys
from tqdm import tqdm


datFN =

ds_2D = xr.open_dataset(dataFN)

def compute_divergence(U,V, dx, dy):
    # divergence: du/dx + dv/dy
    ndim = len(U.dims)
    print('U dimension length = {0:d}'.format(ndim))
    gradU = np.gradient(U, axis = ndim-1)  # along the x axis;
    gradV = np.gradient(V, axis = ndim-2)  # along the y axis;

    div = gradU/dx + gradV/dy
    return div


# test with 2D data, but can apply to 3D wind as well.
U_Surf = ds_2D.USFC
V_Surf = ds_2D.VSFC

dx = 100
dy = 100

div_Surf =compute_divergence(ds_2D.USFC, ds_2D.VSFC, dx, dy)
div_P850 = compute_divergence(ds_2D.U850, ds_2D.V850, dx, dy)



dimlist = list(ds_2D.USFC.dims)
ds_2D['div_SFC'] = (dimlist, div_Surf)
ds_2D['div_SFC'] = ds_2D.div_SFC.assign_attrs(units='s-1', long_name='surface wind divergence')

ds_2D['div_850'] = (dimlist, div_P850)
ds_2D['div_850'] = ds_2D.div_SFC.assign_attrs(units='s-1', long_name='wind divergence at 850mb')

# 3. compute cross-spectrum of these two fields:
from spectrum_analysis.wnspectrum import cross_spectrum
#from spectrum_analysis.wnspectrum import cross_coh2pha

#XX = ds_2D.sel(time=slice(338,338.01)).PW #- ds_2D.CWP.mean(axis=(1,2))
XX = ds_2D.PW[-2:-1,:,:]
#YY = ds_3D.sel(time=slice(338,338.01)).W[:,10,:,:]
YY = ds_2D.div_SFC[-2:-1,:,:]
XX = XX - XX.mean(axis=(1,2))


XX2 = XX.isel(x=slice(0,100))
YY2 = YY.isel(x=slice(0,100))
print(np.shape(XX2), np.shape(YY2))
# the following function will return power spectrum of the two input map, co-spectrum, quadrature spectrum,
# coherence (mag squared), phase and the two components that describe the phase (for plotting). 
results = cross_spectrum(XX2, YY2, dx, dy)


STC = results['STC']
kx = results['kx']
ky = results['ky']
PX = STC[0, :,:,:]
PY = STC[1,:,:,:]
CXY = STC[2,:,:,:]
QXY = STC[3,:,:,:]
COH2 = STC[4,:,:,:]
PHAS = STC[5,:,:,:]


# ready to test the omni-directional spectrum:
from spectrum_analysis.wnspectrum import omnidirectional_wavenumber_spectrum
results2 = omnidirectional_wavenumber_spectrum(PX, kx, ky, dx, dy)

PK = results2['PK']
wn = results2['wvnum']
wlen = results2['wvlen']


fig, ax = plt.subplots(1,1, figsize=(12, 8))
ax.plot(wlen/1e3, PK[1,:], '-r')
plt.gca().invert_xaxis()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(10**2, 10**10)
ax.set_xlim( 10**(-1), 10**2)
ax.set_xlabel('wavelength (km)')
ax.set_ylabel('P(k) ')