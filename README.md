# seagliderHoeLegacy

This repository contains MATLAB routines to analyze two seaglider missions from summer 2015 related to the Hoe Legacy 2 cruises. Missions are g146_m11 and sg512_m06 and data were extracted and binned using the routines in the github repository https://github.com/duebi/seagliderExtraction.
The extractSg\*.m routines create the gridded variables that are used in the other routines (analyzeSg\*.m & dielSg\*.m).

Other Matlab functions that are called in the scripts are found in the common/ subfolder:
- vdist.m: computes distances on the WGS-84 ellipsoid (by Michael Kleder)
- createPatches.m: allows transparency in barplots (modified from the version by Brendan Hamm)

The data/ subfolder contains:
- sg146m11data.mat: gridded data for seaglider 146 mission 11
- sg512m106data.mat: gridded data for seaglider 512 mission 6 
- ccar2015.mat: sea surface height anomaly computed by the Colorado Center for Astrodynamics Research (CCAR) for year 2015
- aviso2015.mat: sea level anomaly computed by AVISO for year 2015
- drifter01.txt: the position of a lagrangian drifter deployed in an anticyclonic eddy North of the Hawaiian island
- hawaii.dat: coastline for the Hawaiian islands
- colorbrewer_anom: some palettes from www.ColorBrewer.org (by Cynthia A. Brewer)

Benedetto Barone - Jan 2016
