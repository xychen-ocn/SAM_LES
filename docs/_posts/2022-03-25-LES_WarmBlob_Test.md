---
layout: single
title: "Research notes for an initial warm blob test case"
date: 2022-03-06 02:00:00 -0000
categories: research-notes
---

# Processing procedures:
### 1. "meridionally" and time averaged stream-wise profiles
 - Two ways:
   1. directly average the y direction regardless of the wind direction.
   2. rebase the coordinate so that the x-axis is align with the domain-, depth-,averaged stream-wise direction. Average results along the span-wise (y) direction. 
      * issue 1: slow regridding 
      * issue 2: the number of data in the span-wise direction varies across the stream-wise direction. 

### 2. spectrum and co-spectrum analysis:
| variables  | Reference |
| ---------- | --------- |
| qtot (total water) | BB2017 Fig 16 |
| ---------- | --------- |
| qtot and buoyant turbulence kinetic energy (TKE_b) |    Pornampai et al.(2021)  |
|/ buoyancy production term in momentum budget equation | Fig 6, and in Appendix. |
| ---------- | --------- |
| w, theta,       | jonker et al. (1999)  |
| top-down, bottom-up scalar  | Fig 2 & Fig 4 |
| ---------- | --------- |
| spectra and co-spectra of  | De Roode (2004) |
| td, bu scalar; w and virtual potential temperature |  Fig 8 | 


# Questions I have:
### 1. How long does it take cloud characteristics to change in response to the SST anomaly forcings?


### 2. What kind of effects do periodic boundary condition have?


### 3. How long does it take for the upstream condition to be affected by the periodic boundary condition?



# Interesting initial results I noticed:

 - Hovmoller diagram: enhanced cloud top height, surface precipitation rate, starting at perhaps 6hr around over the center of the warm anomaly. I need to extend my control run for comparison (I have already run this case for 1.5 days, so I should have enough length of data for this initial comparison.)

 - meridional,time averaged cross-section of 2D variables: cloud top height, PW (water vapor path), CWP both has maximum around the warm blob center. wind speed dropped at the warm blob center. 


# information about the control run: 
-----
| run type  | advection scheme   | domain size and resolution  | 
| control   |  SELPPM (5th order) | 128km x 128km, dx=dy = 100 |

advection scheme: 5th order scheme. 
