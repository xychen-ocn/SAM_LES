##!/usr/bin/env python
#!/home/xychen/env/default/bin/python3
# coding: utf-8

# In[ ]:


# This is an object file that will help to make plotting different variables with the same format easier;
# Four types of plots can be generated.
# Date: 10/13/2021
# X.C.


# In[3]:


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import cm
import os 


# In[1]:


class SAM_DataVisualize():
    def __init__(self):
        self.var = ['CLD','PRECIP', 'RADQR', 'U','V', 'TKE']
        self.nc = None
        self.caseID = None
        self.path = '/data/xchen/SAM_LES_Orion'
        
        
    def domain_mean_evolution(self, varname = None, colormap= 'viridis', zmax= 4000, svfig = False, 
                              svdir = None, figsize = (14,10),figname_suff=None):
        # Purpose: this function is for plotting data in ds_stat to show how it varies in height and time;
        # plot as color shading
        nc = self.nc
        z = self.nc['z']
        time_abs = self.nc['time']
        time = (time_abs - time_abs[0])* 24.  # hours
        var_list = self.var;
        TT, ZZ = np.meshgrid(time, z)
      
    
        
        if (varname is None) and (var_list):
            print('plotting standard variables')
            # plot the default variables:
            nrow = np.int32(np.ceil(len(var_list)/3))
            fig, axes = plt.subplots(nrow,3, figsize=figsize)
            
            for varname, ax in zip(var_list, axes.flatten()):
                #print(varname, ax)
                hm = ax.pcolor(TT, ZZ, nc[varname].T, shading='auto', cmap=colormap)
                # add the inversion height:
                ax.plot(time, nc['ZINV']*1000, '--w',linewidth=1.2)
                ax.set_xlabel('Simulation Time (hrs)')
                ax.set_ylabel('Height (m)')
                ax.set_title(nc[varname].long_name)
                
                # set up the colorbar:
                hb = plt.colorbar(hm, ax=ax)
                hb.ax.set_title(nc[varname].units)
                ax.set_ylim(0, zmax)
            
            fig.tight_layout()
            if (figname_suff is None):
                figname = self.caseID + '_domain_mean_evolution_of_standard_variables.jpg'
            else:
                figname = self.caseID + '_domain_mean_evolution_of_'+ figname_suff +'.jpg'
                
            
        else:
                               
            # plot variable names direclty. 
            fig, ax = plt.subplots(1,1, figsize = (10,8))
            
            hm = ax.pcolor(TT, ZZ, nc[varname].T, cmap=colormap, shading='auto')
            ax.set_xlabel('Simulation Time (hrs)')
            ax.set_ylabel('Height (m)')
            ax.set_title(nc[varname].long_name)

            # set up the colorbar:
            hb = plt.colorbar(hm, ax=ax)
            hb.ax.set_title(nc[varname].units)
            ax.set_ylim(0, zmax)
            
            figname = self.caseID + '_domain_mean_evolution_of_' + varname +'.jpg'
            
        
        # save figure:
        if svfig:
            print('saving figure as ' + figname)
            abs_svpath = os.path.join(self.path, svdir)
            print('to ' + abs_svpath)
            
            # make directory if doesn't exist already:
            if not os.path.exists(abs_svpath):
                os.makedirs(abs_svpath)
            
            abs_figname = os.path.join(abs_svpath, figname)
            fig.savefig(abs_figname)
            
            
    
    def domain_mean_profiles(self, var_list = None, simThr = [0, 24], zmax=4000, ncol = 2, figsize=(12,10), 
                             svfig = False, svdir =None, figname_suff=None):
        # Purpose: this function will plot the domain mean vertical profile 
        # update: change the input argument "tidx" to simThr (simulation hour is easier to specify than the index.)
        nc = self.nc
        z = self.nc['z']      
        ZINV = nc['ZINV']*1000
        time = (self.nc['time']-self.nc['time'][0]) * 24.  # hours
        
         
                
        # plot the input list of variables:
        nrow = np.int32(np.ceil(len(var_list)/ncol))
        
        fig, axes = plt.subplots(nrow, ncol, figsize=figsize)

        for varname, ax in zip(var_list, axes.flatten()):
            #print(varname, ax)
            val = nc[varname]
            
            hline = None
            for t in simThr:
                
                # find the index for the simulation hours:
                #it = np.where(time==t)[0][0]
                # update to make it find the closest time stamp.
                timedif = np.abs(time-t)
                it = np.where(timedif == np.min(timedif))[0][0]
                
                ZINV_it = ZINV[it]

                labelstr = 't={0:.1f}th hr'.format(time[it].values)
                ax.plot(val[it,:],z, linestyle='-',linewidth=1.2, label=labelstr)
                # add the inversion height at t = it
                hzinv = ax.axhline(ZINV_it, linestyle='--',linewidth=1.0, color='gray', label='_nolegend_')
          
            ax.set_xlabel(nc[varname].long_name + '(' + nc[varname].units + ')')
            ax.set_ylabel('Height (m)')
            ax.set_ylim(0, zmax)
            ax.legend()
            


        fig.tight_layout()
        if (figname_suff is None):
            figname = self.caseID + '_domain_mean_vertical_profiles_of_key_parameters.jpg'
        else:
            figname = self.caseID + '_domain_mean_vertical_profiles_of_' + figname_suff +'.jpg'

       
        
        # save figure:
        if svfig:
            print('saving figure as ' + figname)
            abs_svpath = os.path.join(self.path, svdir)
            print('to ' + abs_svpath)
            
            # make directory if doesn't exist already:
            if not os.path.exists(abs_svpath):
                os.makedirs(abs_svpath)
            
            abs_figname = os.path.join(abs_svpath, figname)
            fig.savefig(abs_figname)
            
            
    
    def spatial_map(self, var_list = ['PW','CWP', 'ZC'], zlev = 500, ColIntFlag=False, simThr=[8,16], colormap= cm.Blues, 
                    figsize=(12,10), svfig = False, svdir =None, figname_suff=None):
        # Purpose: this function will plot the spatial distribution of a quantity 
        nc = self.nc
        x = nc['x'] 
        y = nc['y']
        xx, yy = np.meshgrid(x, y)
        time = (nc['time']-nc['time'][0]) * 24.  # hours
        
        # plot the input list of variables:
        ncol = len(simThr)
        nrow = len(var_list)
        
        width_ratio = np.ones(ncol)
        width_ratio[-1] = 1.25
        
        
        fig, axes = plt.subplots(nrow, ncol, figsize=figsize, 
                                 gridspec_kw={'width_ratios':width_ratio})

        # time in coloumns
        # variable in rows:          
        for iv, varname in enumerate(var_list):
            #print(varname, ax)
            # use xarray .sel to select at the desired levels:
            if not ColIntFlag:
                nc_sliced = nc.sel(z=zlev, method = 'nearest')
                val = nc_sliced[varname]
            else:
                val = nc[varname]
            
            for t, ax in zip(simThr, axes[iv].flatten()):
                
                # find the index for the simulation hours:
                #it = np.where(time==t)[0][0]
                timedif = np.abs(time-t)
                it = np.where(timedif == np.min(timedif))[0][0]
                
                labelstr = '{0:.0f}th hour'.format(t)
                hc = ax.contourf(xx, yy, np.squeeze(val[it,:,:]), cmap=colormap
                                )
                # add contours
                cs = ax.contour(xx, yy, np.squeeze(val[it,:,:]), colors = 'k', linewidths=1)
            
                ax.clabel(cs, cs.levels)
                
                if iv==nrow-1:
                    ax.set_xlabel(nc['x'].name + '(' + nc['x'].units + ')')
                if t==simThr[0]:
                    ax.set_ylabel(nc['y'].name + '(' + nc['y'].units + ')')
                if iv==0:
                    ax.set_title(labelstr)
                    
#                 yticks = ax.get_yticks()
#                 ax.set_xticks(yticks)
                    
                #plt.axis('square')
            
            # set colorbar
            # add an axis for colorbar:
            
#             pos1 = ax.get_position()
#             pos2 = [pos1.x0 + 1.15*pos1.width, pos1.y0, pos1.width/10, pos1.height]
    
#             cbar_ax = fig.add_axes(pos2)
            hb = plt.colorbar(hc, ax=ax)
            hb.ax.set_title(val.units)
            hb.ax.set_ylabel(val.long_name)

           

        fig.tight_layout()
        if (figname_suff is None):
            figname = self.caseID + '_spatial_evolution_of_moisture_variables.jpg'
        else:
            figname = self.caseID + '_spatial_evolution_of_' + figname_suff + '.jpg'

       
        
        # save figure:
        if svfig:
            print('saving figure as ' + figname)
            abs_svpath = os.path.join(self.path, svdir)
            print('to ' + abs_svpath)
            
            # make directory if doesn't exist already:
            if not os.path.exists(abs_svpath):
                os.makedirs(abs_svpath)
            
            abs_figname = os.path.join(abs_svpath, figname)
            fig.savefig(abs_figname)
            
            
    
    def quartile_evolution(self, varname=None, tsm=False, tsm_window=None,
                           figsize=(12,10), lw=1.2, svfig = False, svdir = None, figname_suff=None):
     
        # Purpose: this function will plot how the PW sorted profiles evolution in time;
        # this is a 2D line polt.
        
        # define a nested function here:
        def movmean(x, w):
            return np.convolve(x, np.ones(w), 'valid')/w
        
        nc = self.nc
        var_list = self.var
        t = nc.time - nc.time[0]
        quartiles = nc.quartile.values
        
        colors=plt.cm.RdBu(np.linspace(0,1, len(quartiles)))
        
        # plot standard variables defined in the obj self
        if (varname is None) and (var_list):
            print('plotting standard variables')
            # plot the default variables:
            nrow = np.int32(np.ceil(len(var_list)/3))
            fig, axes = plt.subplots(nrow,3, figsize=figsize)
            
            for varname, ax in zip(var_list, axes.flatten()):
                
                for i, qrt in enumerate(quartiles):
                    var = nc.sel(quartile=qrt)[varname]
                    
                    #if smoothing is requested for the time evolution:
                    if (tsm) and (tsm_window is not None):
                        t = movmean(nc.time - nc.time[0], tsm_window)
                        var = movmean(var, tsm_window)

                    ax.plot(t, var, '-', linewidth=lw, color=colors[i], label=qrt)
    
                plt.legend()
                ax.set_xlabel(nc.time.long_name + ' (' + nc.time.units +')' )
                ax.set_ylabel( nc[varname].long_name + ' ('+ nc[varname].units +')' +'\n quartile mean')
            
            
            
            fig.tight_layout()
            if (figname_suff is None):
                figname = self.caseID + '_TWP_quartile_mean_evolution_of_standard_variables.jpg'
            else:
                figname = self.caseID + '_TWP_quartile_mean_evolution_of_'+ figname_suff +'.jpg'
                
            
        else:
            # plot variable names direclty. 
            fig, ax = plt.subplots(1,1, figsize = (10,8))
            
            for i, qrt in enumerate(quartiles):
                var = nc.sel(quartile=qrt)[varname]
                
                #if smoothing is requested for the time evolution:
                if (tsm) and (tsm_window is not None):
                    t = movmean(nc.time - nc.time[0], tsm_window)
                    var = movmean(var, tsm_window)

                ax.plot(t, var, '-', linewidth=lw, color=colors[i], label=qrt)
            
            plt.legend()
            ax.set_xlabel(nc.time.long_name + ' (' + nc.time.units +')' )
            ax.set_ylabel( nc[varname].long_name + ' ('+ nc[varname].units +')' +'\n quartile mean')

            
            figname = self.caseID + '_TWP_quartile_domain_mean_evolution_of_' + varname +'.jpg'
            
        

        
        ### after finish making plots, decide whether or not to save the figures:
        if svfig:
            print('saving figure as ' + figname)
            abs_svpath = os.path.join(self.path, svdir)
            print('to ' + abs_svpath)
            
            # make directory if doesn't exist already:
            if not os.path.exists(abs_svpath):
                os.makedirs(abs_svpath)
            
            abs_figname = os.path.join(abs_svpath, figname)
            fig.savefig(abs_figname)
                
                
                              
              
                
        
    
    
    def quartile_profiles(self, varname = None, colormap = 'viridis', svfig = False, svdir =None):
        # Purpose: this function will plot the verticle profiles from the four PW sorted quartiles
        fig, ax = plt.subplots(2,3, figsize=(12,10))
        
    def spectrum_evolution(self, varname = None, svfig = False, svdir =None):
        # Purpose: this function will plot the evolution of spectrum of a given quantity. (ref. Pornampai et al. 2021)
        fig, ax = plt.subplots(2,3, figsize=(12,10))
        
        
        


# In[ ]:




