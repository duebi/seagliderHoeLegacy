% extractSg146m11.m
%
% script to extract data from 2015 SeaGlider mission sg146m11
% 
% Benedetto Barone - Oct 2015

mission = 'sg146_m11';
d_range = [2 321]; 
upth = userpath; 
sgpath =  [upth(1:end-1) '/Data/seaglider/' mission];
clear upth
%load oxy_cal
%load pcpn_cal

nd = length(d_range(1):d_range(2)); % total number of dives

% 1. Collect data from all dives:
[UP,DWN,~] = sg_cat(d_range,sgpath,'userVarNames',{'vmtime','press','sigmath0','salin', ...
    'tempc','oxygen','optode_oxygen','optode_temp','optode_dphase_oxygen','optode_dphase','wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'});

% 2. Rename WETLabs variables
DWN.Properties.VariableNames({'wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'}) = {'bb470','bb700','chl1','bb660','cdom','chl2'};
UP.Properties.VariableNames({'wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'}) = {'bb470','bb700','chl1','bb660','cdom','chl2'};

% 3. Clean backscattering data based on 3 std on a 20m grid (sets outliers = NaN)
for i = d_range(1):d_range(2)
    ind_d = DWN.divenum==i & ~isnan(DWN.bb470);
    [bin_d,bin,ind_good_d] = binning(DWN.bb470(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb470(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb470(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb470);
    [bin_u,bin,ind_good_u] = binning(UP.bb470(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb470(ind_u); temp_u(~ind_good_u) = NaN; UP.bb470(ind_u) = temp_u;
    clear ind_d ind_u ind_good_d ind_good_u temp_d temp_u bin bin_d bin_u
    ind_d = DWN.divenum==i & ~isnan(DWN.bb660);
    [bin_d,bin,ind_good_d] = binning(DWN.bb660(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb660(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb660(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb660);
    [bin_u,bin,ind_good_u] = binning(UP.bb660(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb660(ind_u); temp_u(~ind_good_u) = NaN; UP.bb660(ind_u) = temp_u;
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
betasw470d = NaN(ldwn,1); betasw660d = NaN(ldwn,1); betasw700d = NaN(ldwn,1);
betasw470u = NaN(lup,1); betasw660u = NaN(lup,1); betasw700u = NaN(lup,1);
    % compute scatter at 117 deg due to seawater
for i = 1:ldwn % downcast 
    [betasw470d(i),~,~]= betasw_ZHH2009(470,DWN.tempc(i),117,DWN.salin(i)); % unfortunately this function doesn't work with vectors in T or S
    [betasw660d(i),~,~]= betasw_ZHH2009(660,DWN.tempc(i),117,DWN.salin(i));
    [betasw700d(i),~,~]= betasw_ZHH2009(700,DWN.tempc(i),117,DWN.salin(i));
end
for i = 1:lup % upcast
    [betasw470u(i),~,~]= betasw_ZHH2009(470,UP.tempc(i),117,UP.salin(i)); % unfortunately this function doesn't work with vectors in T or S
    [betasw660u(i),~,~]= betasw_ZHH2009(660,UP.tempc(i),117,UP.salin(i));
    [betasw700u(i),~,~]= betasw_ZHH2009(700,UP.tempc(i),117,UP.salin(i));
end
    % calculation for particle backscattering coefficients
DWN.bbp470 = 2*pi*chi_p*(DWN.bb470-betasw470d);
DWN.bbp660 = 2*pi*chi_p*(DWN.bb660-betasw660d);
DWN.bbp700 = 2*pi*chi_p*(DWN.bb700-betasw700d);
UP.bbp470 = 2*pi*chi_p*(UP.bb470-betasw470u);
UP.bbp660 = 2*pi*chi_p*(UP.bb660-betasw660u);
UP.bbp700 = 2*pi*chi_p*(UP.bb700-betasw700u);
clear ldwn lup sgpath
clear chi_p betasw470d betasw470u betasw660d betasw660u betasw700d betasw700u

% 5. Correct optode oxygen for salinity and pressure (because dphase is missing in this mission)
    % Aandera salinity correction coefficients (sensor independent)
aa_sb = [-6.24097E-3 -6.93498E-3 -6.90358E-3 -4.29155E-3];
aa_sc = -3.11680E-7;
    % Salinity correction
scal_t = log((298.15 - DWN.optode_temp)./(273.15 + DWN.optode_temp));
DWN.optode_oxygen =  DWN.optode_oxygen.*exp(DWN.salin.*(aa_sb(1)+aa_sb(2)*scal_t+aa_sb(3)*scal_t.^2+aa_sb(4)*scal_t.^3)+aa_sc.*DWN.salin.^2);
scal_t = log((298.15 - UP.optode_temp)./(273.15 + UP.optode_temp));
UP.optode_oxygen =  UP.optode_oxygen.*exp(UP.salin.*(aa_sb(1)+aa_sb(2)*scal_t+aa_sb(3)*scal_t.^2+aa_sb(4)*scal_t.^3)+aa_sc.*UP.salin.^2);
    % Pressure correction
DWN.optode_oxygen = DWN.optode_oxygen.*(1+DWN.press.*0.032./1000);
UP.optode_oxygen = UP.optode_oxygen.*(1+UP.press.*0.032./1000);

% 6. Save concatenated dive file
save([mission '_cat'],'UP','DWN');

% 7. Put data on a regular grid
datafile = [mission(1:5) mission(7)  mission(8:9) 'data'];
depth = 2:2:1000; ld = length(depth);
[UG,DG] = sg_grid(['./' mission '_cat'],depth,'gridVar','vmdepth','diveRange',d_range(1):d_range(2),'outVars',{'lon','lat','daten','salin','tempc','sigmath0','oxygen','optode_dphase_oxygen','optode_oxygen','chl1','chl2','bbp470','bbp660','bbp700','cdom'});
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
sgd.o = reshape(DG.oxygen,ld,nd);
sgd.opt = reshape(DG.optode_oxygen,ld,nd); % The dphase is missing for this cruise, so I had to use the data that have been computed internally
sgd.chl1 = reshape(DG.chl1,ld,nd);
sgd.chl2 = reshape(DG.chl2,ld,nd);
sgd.bbp470 = reshape(DG.bbp470,ld,nd);
sgd.bbp660 = reshape(DG.bbp660,ld,nd);
sgd.bbp700 = reshape(DG.bbp700,ld,nd);
sgd.cdom = reshape(DG.cdom,ld,nd);
    % Dive information
dived.lon = nanmean(sgd.lon);
dived.lat = nanmean(sgd.lat);
dived.dive = d_range(1):d_range(2);
dived.date = nanmean(sgd.date);
dived.hour = dived.date - fix(dived.date);

% Save new variables sgd and dived
    save(datafile,'sgd','dived')