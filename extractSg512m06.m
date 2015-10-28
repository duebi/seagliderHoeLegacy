% extractSg512m06.m
%
% script to extract data from 2015 SeaGlider mission sg512m06
% 
% Benedetto Barone - Oct 2015

mission = 'sg512_m06';
d_range = [2 581]; 
upth = userpath; 
sgpath =  [upth(1:end-1) '/Data/seaglider/' mission];
clear upth
load([sgpath '/oxy_cal'])
%load pcpn_cal

nd = length(d_range(1):d_range(2)); % total number of dives

% 1. Collect data from all dives:
[UP,DWN,~] = sg_cat(d_range,sgpath,'userVarNames',{'vmtime','press','sigmath0','salin', ...
    'tempc','oxygen','optode_oxygen','optode_dphase_oxygen','optode_calphase','wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'});

% 2. Rename WETLabs variables
DWN.Properties.VariableNames({'wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'}) = {'bb470','bb700','chl1','bb650','cdom','chl2'};
UP.Properties.VariableNames({'wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'}) = {'bb470','bb700','chl1','bb650','cdom','chl2'};

% 3. Clean backscattering data based on 3 std on a 20m grid (sets outliers = NaN)
for i = d_range(1):d_range(2)
    ind_d = DWN.divenum==i & ~isnan(DWN.bb470);
    [bin_d,bin,ind_good_d] = binning(DWN.bb470(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb470(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb470(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb470);
    [bin_u,bin,ind_good_u] = binning(UP.bb470(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb470(ind_u); temp_u(~ind_good_u) = NaN; UP.bb470(ind_u) = temp_u;
    clear ind_d ind_u ind_good_d ind_good_u temp_d temp_u bin bin_d bin_u
    ind_d = DWN.divenum==i & ~isnan(DWN.bb650);
    [bin_d,bin,ind_good_d] = binning(DWN.bb650(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb650(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb650(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb650);
    [bin_u,bin,ind_good_u] = binning(UP.bb650(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb650(ind_u); temp_u(~ind_good_u) = NaN; UP.bb650(ind_u) = temp_u;
    clear ind_d ind_u ind_good_d ind_good_u temp_d temp_u bin bin_d bin_u
    ind_d = DWN.divenum==i & ~isnan(DWN.bb700);
    [bin_d,bin,ind_good_d] = binning(DWN.bb700(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb700(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb700(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb700);
    [bin_u,bin,ind_good_u] = binning(UP.bb700(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb700(ind_u); temp_u(~ind_good_u) = NaN; UP.bb700(ind_u) = temp_u;
    %{
    plot(UP.bb470(ind_u),UP.vmdepth(ind_u),'k.',DWN.bb470(ind_d),DWN.vmdepth(ind_d),'r.')
    hold on, plot(bin_u,bin,'k',bin_d,bin,'r'), hold off
    set(gca,'ydir','rev','ylim',[0 200])
    pause
    %}
    clear ind_d ind_u ind_good_d ind_good_u temp_d temp_u bin bin_d bin_u
end

% 4. Compute particle backscattering coefficients (this could be moved into sg_divecalc.m)
chi_p = 1.1;
ldwn = height(DWN); lup = height(UP); 
betasw470d = NaN(ldwn,1); betasw650d = NaN(ldwn,1); betasw700d = NaN(ldwn,1);
betasw470u = NaN(lup,1); betasw650u = NaN(lup,1); betasw700u = NaN(lup,1);
    % compute scatter at 117 deg due to seawater
for i = 1:ldwn % downcast 
    [betasw470d(i),~,~]= betasw_ZHH2009(470,DWN.tempc(i),117,DWN.salin(i)); % unfortunately this function doesn't work with vectors in T or S
    [betasw650d(i),~,~]= betasw_ZHH2009(650,DWN.tempc(i),117,DWN.salin(i));
    [betasw700d(i),~,~]= betasw_ZHH2009(700,DWN.tempc(i),117,DWN.salin(i));
end
for i = 1:lup % upcast
    [betasw470u(i),~,~]= betasw_ZHH2009(470,UP.tempc(i),117,UP.salin(i)); % unfortunately this function doesn't work with vectors in T or S
    [betasw650u(i),~,~]= betasw_ZHH2009(650,UP.tempc(i),117,UP.salin(i));
    [betasw700u(i),~,~]= betasw_ZHH2009(700,UP.tempc(i),117,UP.salin(i));
end
    % calculation for particle backscattering coefficients
DWN.bbp470 = 2*pi*chi_p*(DWN.bb470-betasw470d);
DWN.bbp650 = 2*pi*chi_p*(DWN.bb650-betasw650d);
DWN.bbp700 = 2*pi*chi_p*(DWN.bb700-betasw700d);
UP.bbp470 = 2*pi*chi_p*(UP.bb470-betasw470u);
UP.bbp650 = 2*pi*chi_p*(UP.bb650-betasw650u);
UP.bbp700 = 2*pi*chi_p*(UP.bb700-betasw700u);
clear ldwn lup sgpath
clear chi_p betasw470d betasw470u betasw650d betasw650u betasw700d betasw700u

% 5. Save concatenated dive file
save([mission '_cat'],'UP','DWN');
 
% 6. Put data on a regular grid
datafile = [mission(1:5) mission(7)  mission(8:9) 'data'];
depth = 2:2:1000; ld = length(depth);
[UG,DG] = sg_grid([sgpath '/' mission '_cat'],depth,'gridVar','vmdepth','diveRange',d_range(1):d_range(2),'outVars',{'lon','lat','daten','salin','tempc','sigmath0','oxygen','optode_dphase_oxygen','optode_oxygen','chl1','chl2','bbp470','bbp650','bbp700','cdom'});
    % Build new table for downcast
    % Time & Position
sgd = table(depth',reshape(DG.daten,ld,nd)); sgd.Properties.VariableNames = ({'depth','date'});
sgd.date = sgd.date-10/24; %to transform in HST time
sgd.hour = sgd.date - fix(sgd.date);
sgd.lon = reshape(DG.lon,ld,nd);
sgd.lat = reshape(DG.lat,ld,nd);
    % Water column measurements
sgd.s = reshape(DG.salin,ld,nd);
sgd.t = reshape(DG.tempc,ld,nd);
sgd.sig = reshape(DG.sigmath0,ld,nd);
sgd.o = reshape(DG.oxygen,ld,nd).*((sgd.sig+1000)/1000); % From umol kg-1 to umol L-1
sgd.opt = reshape(DG.optode_dphase_oxygen,ld,nd).*((sgd.sig+1000)/1000); % From umol kg-1 to umol L-1
sgd.chl1 = reshape(DG.chl1,ld,nd);
sgd.chl2 = reshape(DG.chl2,ld,nd);
sgd.bbp470 = reshape(DG.bbp470,ld,nd);
sgd.bbp650 = reshape(DG.bbp650,ld,nd);
sgd.bbp700 = reshape(DG.bbp700,ld,nd);
sgd.cdom = reshape(DG.cdom,ld,nd);
sgd.Properties.VariableUnits = {'m','Matlab date (HST)','Decimal day (HST)','Degrees E','Degrees N'...
    'Absolute salinity', 'degrees Celsius', 'kg m-3','umol L-1','umol L-1','rel. units','rel. units' ...
    'm-1','m-1','m-1', 'rel. units'};
    % Dive information
dived.lon = nanmean(sgd.lon);
dived.lat = nanmean(sgd.lat);
dived.dive = d_range(1):d_range(2);
dived.date = nanmean(sgd.date);
dived.hour = dived.date - fix(dived.date);
dived.units = {'lon: Degrees E','lat: Degrees N','dive: #','date: Matlab date (HST)','hour: Decimal day (HST)',};

% 7. Oxygen calibration on Winkler measurements
    % build comparison array
o2_comp = oxy_cal;
o2_comp.dist = NaN(height(oxy_cal),1);
o2_comp.lag = NaN(height(oxy_cal),1);
o2_comp.zdist = NaN(height(oxy_cal),1);
o2_comp.optode = NaN(height(oxy_cal),1);
o2_comp.sbe43 = NaN(height(oxy_cal),1);
for i = 1:height(oxy_cal)
    [o2_comp.lag(i),ind_t] = min(abs(dived.date-o2_comp.date(i)));
    o2_comp.dist(i) = vdist(dived.lat(ind_t),dived.lon(ind_t),o2_comp.lat(i),-o2_comp.lon(i));
    [o2_comp.zdist(i),ind_z] = min(abs(sgd.depth-o2_comp.depth(i)));
    o2_comp.optode(i) = sgd.opt(ind_z,ind_t);
    o2_comp.sbe43(i) = sgd.o(ind_z,ind_t);
    clear ind_t ind_z
end
ind_good = o2_comp.dist<10000 & o2_comp.lag<0.2 & o2_comp.zdist<10;
ind_correction = o2_comp.dist<10000 & o2_comp.lag<0.2 & o2_comp.zdist<10 & o2_comp.depth<100;
    % Compute offset in the upper 100 m and correct oxygen from sensors
optode_offset = nanmedian(o2_comp.oxygen(ind_correction)-o2_comp.optode(ind_correction));
sbe43_offset = nanmedian(o2_comp.oxygen(ind_correction)-o2_comp.sbe43(ind_correction));
sgd.opt = sgd.opt + optode_offset;
sgd.o = sgd.o + sbe43_offset;
    % Plot comparison results
%{
subplot(2,2,1)
scatter(o2_comp.oxygen(ind_good),o2_comp.optode(ind_good),200,o2_comp.depth(ind_good),'.')
caxis([0 500]), colormap(jet), cb = colorbar, title(cb,'depth (m)','Fontsize',16)
ylim([80 240]),xlim([80 240])
hold on, plot([80 240],[80 240],'k--'), hold off
set(gca,'Fontsize',16), xlabel('Winkler oxygen (umol L-1)'), ylabel('Optode oxygen (umol L-1)')
subplot(2,2,2)
scatter(o2_comp.oxygen(ind_good),o2_comp.sbe43(ind_good),200,o2_comp.depth(ind_good),'.')
caxis([0 500]), colormap(jet), cb = colorbar, title(cb,'depth (m)','Fontsize',16)
ylim([80 240]),xlim([80 240])
hold on, plot([80 240],[80 240],'k--'), hold off
set(gca,'Fontsize',16), xlabel('Winkler oxygen (umol L-1)'), ylabel('SBE43 oxygen (umol L-1)')
subplot(2,2,3)
hist(o2_comp.oxygen(ind_correction)-o2_comp.optode(ind_correction),20)
set(gca,'Fontsize',16), xlabel('Winkler - optode (umol L-1)'), ylabel('n. obs.')
title('upper 100m')
subplot(2,2,4)
hist(o2_comp.oxygen(ind_correction)-o2_comp.sbe43(ind_correction),20)
set(gca,'Fontsize',16), xlabel('Winkler - SBE43 (umol L-1)'), ylabel('n. obs.')
title('upper 100m')
%}

% Save new variables sgd and dived
    save(datafile,'sgd','dived')