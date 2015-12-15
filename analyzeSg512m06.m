% analyzeSg512m06.m
%
% script to analyze data from 2015 SeaGlider mission sg512m06
% 
% Benedetto Barone - Oct 2015

mission = 'sg512_m06';
upth = userpath; 
sgpath =  [upth(1:end-1) '/Data/seaglider/' mission];
cd(sgpath);
clear upth sgpath mission
load sg512m06data
load ../../colorbrewer_anom
load ../aviso2015 % load SLA from AVISO
load ../hawaii.dat
coastline = hawaii;
coastline(:,1) = -coastline(:,1);
ind_n = find(isnan(coastline(:,1)));
clear hawaii

d1 = readtable('../drifter01.txt'); % load lagrangian drifter position
d1.date = datenum(d1.DeviceDateTime)-10/24;
% compute glider distance from drifter (ddist)
dlon = interp1(d1.date,d1.Longitude,dived.date);
dlat = interp1(d1.date,d1.Latitude,dived.date);
ddist = vdist(dived.lat,dived.lon,dlat,dlon);

% indeces of different transects
merid1 = dived.dive >= 113 & dived.dive <= 178 & dived.dive ~= 143; % first meridional transect
shortz1 = dived.dive >= 178 & dived.dive <= 206; % short zonal transect
zonal1 = dived.dive >= 206 & dived.dive <= 242; % first zonal transect
lagr1 = dived.dive >= 245 & dived.dive <= 405; % lagrangian period (following drifter)

% Find mixed layer depth and density at the base of the mixed layer (0.03 kg m-3 difference from sigma at 4 m depth)
mld = sgd.sig - repmat(sgd.sig(2,:),height(sgd),1) - 0.03;
mld(mld<0) = NaN; mld(1,:) = NaN;
[sig003,ind003] = nanmin(mld);
mld003 = sgd.depth(ind003);
mld003sig = sig003 + sgd.sig(2,:) + 0.03;

% Compute Brunt-Vaisala frequency
rho_z = (sgd.sig(3:end,:)-sgd.sig(1:end-2,:))./repmat(sgd.depth(3:end)-sgd.depth(1:end-2),1,nd);
rho_z = [NaN(1,nd); rho_z; NaN(1,nd)];
rho_z(rho_z<0) = NaN;
bvf = sqrt((9.8*rho_z)./(1000+sgd.sig)); % s-1


%% Lagrangian period
ind_part = lagr1;
subplot(7,1,2:3)
contourf(dived.date(ind_part),sgd.depth,sgd.opt(:,ind_part),140:1.25:230,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
datetick('x','mm/dd'), xlim([min(dived.date(ind_part)) max(dived.date(ind_part))])
ylim([0 200])
hold on, contour(dived.date(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Date mm/dd 2015'),ylabel('Depth (m)')
caxis([192.5 227.5])
cb = colorbar, title(cb,'Oxygen (umol L-1)'), set(cb,'Fontsize',16);
subplot(7,1,4:5)
contourf(dived.date(ind_part),sgd.depth,sgd.s(:,ind_part),30:0.025:36,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
datetick('x','mm/dd'), xlim([min(dived.date(ind_part)) max(dived.date(ind_part))])
ylim([0 200])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.date(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Date mm/dd 2015'),ylabel('Depth (m)')
caxis([34.8 35.4])
cb = colorbar, title(cb,'Salinity'), set(cb,'Fontsize',16);
subplot(7,1,6:7)
contourf(dived.date(ind_part),sgd.depth,sgd.chl2(:,ind_part),0:25:1000,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
datetick('x','mm/dd'), xlim([min(dived.date(ind_part)) max(dived.date(ind_part))])
ylim([0 200])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.date(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Date mm/dd 2015'),ylabel('Depth (m)')
caxis([50 500])
cb = colorbar, title(cb,'Chlorophyll a (ng L-1)'), set(cb,'Fontsize',16);
pp = get(gca,'Position');
%{
subplot(7,1,6:7)
contourf(dived.date(ind_part),sgd.depth,sgd.bbp470(:,ind_part),0:5e-4:1e-2,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
datetick('x','mm/dd'), xlim([min(dived.date(ind_part)) max(dived.date(ind_part))])
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.date(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Date mm/dd 2015'),ylabel('Depth (m)')
caxis([5e-4 4e-3])
cb = colorbar, title(cb,'bbp 470 nm (m-1)'), set(cb,'Fontsize',16);
pp = get(gca,'Position');
%}
subplot(7,1,1)
plot(dived.date(ind_part),dived.sla(ind_part),'k-'), xlabel('date mm/dd 2015'),ylabel('SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
datetick('x','mm/dd'), xlim([min(dived.date(ind_part)) max(dived.date(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)]);
colormap(jet(200))

%% Long meridional transect
ind_part = merid1;
subplot(7,1,2:3)
contourf(dived.lat(ind_part),sgd.depth,sgd.opt(:,ind_part),140:1.25:230,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([192.5 227.5])
cb = colorbar, title(cb,'Oxygen (umol L-1)'), set(cb,'Fontsize',16);
subplot(7,1,4:5)
contourf(dived.lat(ind_part),sgd.depth,sgd.s(:,ind_part),30:0.05:36,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([34.1 35.4])
cb = colorbar, title(cb,'Salinity'), set(cb,'Fontsize',16);
subplot(7,1,6:7)
contourf(dived.lat(ind_part),sgd.depth,sgd.chl2(:,ind_part),0:25:1000,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([50 450])
cb = colorbar, title(cb,'chlorophyll a (ng L-1)'), set(cb,'Fontsize',16);
pp = get(gca,'Position');
%{
subplot(7,1,6:7)
contourf(dived.lat(ind_part),sgd.depth,sgd.bbp470(:,ind_part),0:5e-4:1e-2,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([5e-4 4e-3])
cb = colorbar, title(cb,'bbp 470 nm (m-1)'), set(cb,'Fontsize',16);
pp = get(gca,'Position');
%}
subplot(7,1,1)
plot(dived.lat(ind_part),dived.sla(ind_part),'k-'), xlabel('Longitude E'),ylabel('SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
xlim([min(dived.lat(ind_part)) max(dived.lat(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)]);
colormap(jet(200))


%% Long zonal transect
ind_part = zonal1;
subplot(7,1,2:3)
contourf(dived.lon(ind_part),sgd.depth,sgd.o(:,ind_part),140:2.5:230,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
hold on, contour(dived.lon(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Longitude E'),ylabel('Depth (m)')
caxis([192.5 225])
cb = colorbar, title(cb,'Oxygen (umol L-1)'), set(cb,'Fontsize',16);
subplot(7,1,4:5)
contourf(dived.lon(ind_part),sgd.depth,sgd.s(:,ind_part),30:0.05:36,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lon(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Longitude E'),ylabel('Depth (m)')
caxis([34.1 35.45])
cb = colorbar, title(cb,'Salinity'), set(cb,'Fontsize',16);
subplot(7,1,6:7)
contourf(dived.lon(ind_part),sgd.depth,sgd.bbp470(:,ind_part),0:5e-4:1e-2,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lon(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Longitude E'),ylabel('Depth (m)')
caxis([5e-4 4e-3])
cb = colorbar, title(cb,'bbp 470 nm (m-1)'), set(cb,'Fontsize',16);
pp = get(gca,'Position');
subplot(7,1,1)
plot(dived.lon(ind_part),dived.sla(ind_part),'k-'), xlabel('Longitude E'),ylabel('SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
xlim([min(dived.lon(ind_part)) max(dived.lon(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)]);
colormap(jet(200))


%% Map of lagrangian period

% indeces for the central dates of the transects
ind_aug29 = sla.date == datenum(2015,8,29);
ind_sep18 = sla.date == datenum(2015,9,18);

ind_drift_lag = d1.date >= datenum(2015,8,29) & d1.date <= datenum(2015,9,18);

subplot(1,2,1)
contourf(sla.lon,sla.lat,sla.sla(:,:,ind_aug29)*100,-26:4:26,'edgecolor','none')
title('Start of Lagrangian period - Aug 29')
colormap(anom_map2)
cb = colorbar; set(cb,'Fontsize',16), title(cb,'SLA (cm)')
caxis([-22 22])
set(cb,'Ticks',[-20:4:20])
xlabel('Longitude (deg E)'), ylabel('Latitude (deg N)')
% Station ALOHA & coastline
hold on
    viscircles([-158 22.75],6/59.79,'edgecolor','k','linewidth',2);
    text(-159,22.5,'  Sta. ALOHA','color','k','Fontsize',16);
    for i = 1:(length(ind_n)-1)
        patch(-coastline((ind_n(i)+1):(ind_n(i+1)-1),1),coastline((ind_n(i)+1):(ind_n(i+1)-1),2),'k','edgecolor','none')
    end
    plot(dived.lon,dived.lat,'Color',[0.7 0.7 0.7])
    plot(d1.Longitude(ind_drift_lag),d1.Latitude(ind_drift_lag),'w-','Linewidth',2)
    plot(dived.lon(lagr1),dived.lat(lagr1),'k--','LineWidth',3)
hold off
set(gca,'Fontsize',16)
axis equal
xlim([-160 -156]),ylim([21.5 26])

subplot(1,2,2)
contourf(sla.lon,sla.lat,sla.sla(:,:,ind_sep18)*100,-26:4:26,'edgecolor','none')
title('End of Lagrangian period - Sep 18')
colormap(anom_map2)
cb = colorbar; set(cb,'Fontsize',16), title(cb,'SLA (cm)')
caxis([-22 22])
set(cb,'Ticks',[-20:4:20])
xlabel('Longitude (deg E)'), ylabel('Latitude (deg N)')
% Station ALOHA & coastline
hold on
    viscircles([-158 22.75],6/59.79,'edgecolor','k','linewidth',2);
    text(-159,22.5,'  Sta. ALOHA','color','k','Fontsize',16);
    for i = 1:(length(ind_n)-1)
        patch(-coastline((ind_n(i)+1):(ind_n(i+1)-1),1),coastline((ind_n(i)+1):(ind_n(i+1)-1),2),'k','edgecolor','none')
    end
    plot(dived.lon,dived.lat,'Color',[0.7 0.7 0.7])
    plot(d1.Longitude(ind_drift_lag),d1.Latitude(ind_drift_lag),'w-','Linewidth',2)
    plot(dived.lon(lagr1),dived.lat(lagr1),'k--','LineWidth',3)
hold off
set(gca,'Fontsize',16)
axis equal
xlim([-160 -156]),ylim([21.5 26])