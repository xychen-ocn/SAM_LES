#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits import mplot3d
from matplotlib.ticker import LinearLocator
from matplotlib import cm

import os
from tqdm import tqdm


# In[31]:


class SAM_DataProcess():
    def __init__(self):
        self.var = ['PW','CWP','TWP','ZC','TB']
        self.nc = None
        self.caseID = None
        self.path = '/data/xchen/SAM_LES_Orion'
             
    def block_averaging(self, nblock=64, domain_shape='square'):
        """ Purpose: this method is used to compute block average. 
            note: I think I need to connect the first point with the end point as well 
            (repeat the first column of points after the last column) for periodic boundary condition.
            May not have significant influence on the ultimate result.
        """
        nc = self.nc
        
        # pad extra grid point column at the end of x and y (to represent periodic boundary condition)
        ds_firstcolumn = nc.isel(x=0)
        ds_firstcolumn.coords['x']=nc.x[-1].values+nc.x[1].values-nc.x[0].values
        ds_firstcolumn.expand_dims(x=1)
        nc_extx = xr.concat([nc, ds_firstcolumn],"x")
        
        ds_firstrow = nc_extx.isel(y=0)
        ds_firstrow.coords['y']=nc_extx.y[-1].values+nc_extx.y[1].values-nc_extx.y[0].values
        ds_firstrow.expand_dims(y=1)
        ds_2D = xr.concat([nc_extx, ds_firstrow],"y")
        
        
        # add one extra column of grid point at the end due to periodic boundary condition (?)
        domain_size = np.max(ds_2D.x) - np.min(ds_2D.x) 
        
        # we know a priori that the domain is a square;
        nblock_per_side = np.sqrt(nblock)
        block_size = domain_size.values/nblock_per_side

        # we also know that the domain resolution is horizontally uniform
        dx = ds_2D.x[1] - ds_2D.x[0]
        dy = ds_2D.y[1] - ds_2D.y[0]

        # construct coodinate for the 64 blocks:
        xincs = np.int32(np.round(block_size/dx.values))
        n = np.arange(0, nblock_per_side, dtype=int)
        
        xidx_block = np.int8(np.zeros(np.shape(n), dtype=int) + n*xincs + xincs/2*np.ones(np.shape(n), dtype=int))
        yidx_block = xidx_block
        
        ds_sub= ds_2D.isel(x=xidx_block, y=yidx_block)
        
        #xidx_vtx = np.int8(np.zeros(np.shape(n), dtype=int) + n*xincs)
        #yidx_vtx = xidx_vtx

        #ds_vtx = ds_2D.isel(x=xidx_vtx, y=yidx_vtx)

        #x2d, y2d = np.meshgrid(ds_2D.x, ds_2D.y)
        #xb2d, yb2d = np.meshgrid(ds_sub.x, ds_sub.y)
        
        # average data point within each block;
        r = block_size/2
        ds_blockave = [None]*nblock
        cnt = 0
        for i in np.arange(nblock_per_side, dtype=int):
            for j in np.arange(nblock_per_side, dtype=int):        
                critx = (ds_2D.x>=ds_sub.x[i]-r) * (ds_2D.x<=ds_sub.x[i]+r)   # use math operator instead of and/or
                crity = (ds_2D.y>=ds_sub.y[j]-r) * (ds_2D.y<=ds_sub.y[j]+r)
                crit_block = critx * crity
                ds_tmp = ds_2D.where(crit_block, drop=True)

                # find block average for all the variables:
                # how to set up a new dimenion named block
                ds_blockave_tmp = ds_tmp.mean(dim=['x','y'])
                ds_blockave_tmp = ds_blockave_tmp.assign_coords(block=cnt)
                ds_blockave[cnt] = ds_blockave_tmp.expand_dims('block')
                cnt = cnt +1
                
        # combined all the blocks into one xarray dataset.
        ds_blkave = xr.combine_by_coords(ds_blockave)
        
        return ds_blkave
    
    
    def TWPsorted_quartile_statistics(self,TWP_bins=None, TWP_bin_labels=['Q1','Q2','Q3','Q4'] ):
        
        ds_blkave = self.nc
        
        # sort the blocks by TWP 
        #TWP_bins=np.percentile(ds_blkave.TWP, np.linspace(0,100,5))
        #TWP_bin_labels = ['Q1','Q2','Q3','Q4']
        
        ds_quartile_mean = [None]*len(ds_blkave.time)
        for ip in tqdm(range(10)):
            for it in np.arange(len(ds_blkave.time)):
                # get the time instance at current time:
                ds_blkave_tmp = ds_blkave.isel(time=it)

                # group blocks into four quartiles:
                TWP_bins=np.percentile(ds_blkave_tmp.TWP, np.linspace(0,100,5))
                ds_blkave_PW_grped = ds_blkave_tmp.groupby_bins('TWP', TWP_bins, 
                                                       labels=TWP_bin_labels).groups

                # compute the quartile averaged value for all variables:
                cnt=0
                ds_QN_mean = [None]*4
                for label, grpID in ds_blkave_PW_grped.items():
                    #print(label)
                    ds_QN = ds_blkave_tmp.isel(block=grpID)
                    # average over all the selected blocks:
                    ds_tmp =ds_QN.mean(dim='block')
                    ds_tmp = ds_tmp.assign_coords(quartile=label)
                    ds_QN_mean[cnt] = ds_tmp.expand_dims('quartile')
                    cnt=cnt+1

                ds_quartile_mean[it] = xr.combine_by_coords(ds_QN_mean)

        # concatenate data:
        ds_qrts = xr.concat(ds_quartile_mean, dim='time')
        return ds_qrts

