# seagliderHoeLegacy

This repository contains MATLAB routines to extract and analyze two seaglider missions from summer 2015 related to the Hoe Legacy 2 cruises. Missions are g146_m11 and sg512_m06.

Routines in the main folder require the private routines for seaglider data extraction at https://github.com/whoi-glider/glider-kit  (David Nicholson). Other Matlab functions that are called in the scripts are found in the common/ subfolder:
- betasw_ZHH2009: computes scattering by pure seawater (by Xiaodong Zhang)
- vdist: computes distances on the WGS-84 ellipsoid (by Michael Kleder)
- binning: binning routine (by Benedetto Barone)

The data/ subfolder contains:
- ccar2015.mat: sea surface height anomaly computed by the Colorado Center for Astrodynamics Research (CCAR) for year 2015
- aviso2015.mat: sea level anomaly computed by AVISO for year 2015
- drifter01.txt: the position of a lagrangian drifter deployed in an anticyclonic eddy North of the Hawaiian island
- oxy_cal.mat: Winkler oxygen measurements from several HOT cruises. These data are used for sensor calibration.

Benedetto Barone - Oct 2015
