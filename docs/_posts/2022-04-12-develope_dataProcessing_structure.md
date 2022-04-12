---
layout: single
title: "Directory Structure for saving data"
date: 2022-04-12 10:40:00 -0000
categories: research-notes
---

```
data_output_dir: /work2 (or work1) disk space 
|
|---- caseFolder1 (naming convention: SAM version + campaign case: e.g. SAM_UW_RICO)
|     |  
|     |--- OUT_STAT 
|     |     + testIDs.stat
|     |     + testIDs.nc (post-processing by SAM)
|     |--- OUT_2D
|     |     + testIDs.2Dbin
|     |     + testIDs.2Dbin_1.nc (post-processing by SAM)
|     |--- OUT_3D
|     |     + testIDs_timestep.bin3D
|     |     + testIDs_timestep.nc (post-processing by SAM)
|     |--- RESTART
|     |     + testIDs.restart (a bunch of files for each case)
|     |--- OBJ (compiled code stored here)
|     |--- OUT_MOVIES/OUT_MOMENTS (usually empty)
|     |  
|     |--- Figs
|     |     + **caseID_figname.jpg (currently, could put them is separate subfolders)
|     |     |--- subfolder (named by caseID)
|     |     |     + figname01.jpg
|     |     |     + figname0N.jpg
|     |  
|     |--- processed_data
|     |     |--- subfolder (named by caseID)
|     |     |     + processed_type01.nc (created by user written python scripts)
|     |     |     

