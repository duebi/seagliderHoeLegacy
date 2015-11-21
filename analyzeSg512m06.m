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
load ../aviso2015 % load SLA from AVISO

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

%% Long meridional transect
ind_part = lagr1;
subplot(7,1,2:3)
contourf(dived.lat(ind_part),sgd.depth,sgd.o(:,ind_part),140:2.5:230,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([192.5 225])
cb = colorbar, title(cb,'Oxygen (umol L-1)'), set(cb,'Fontsize',16);
subplot(7,1,4:5)
contourf(dived.lat(ind_part),sgd.depth,sgd.s(:,ind_part),30:0.05:36,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([34.1 35.45])
cb = colorbar, title(cb,'Salinity'), set(cb,'Fontsize',16);
subplot(7,1,6:7)
contourf(dived.lat(ind_part),sgd.depth,sgd.chl1(:,ind_part),-0.1:0.025:0.9,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([0 0.6])
cb = colorbar, title(cb,'bbp 470 nm (m-1)'), set(cb,'Fontsize',16);
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

%% Lagrangian period
ind_part = lagr1;
subplot(7,1,2:3)
contourf(dived.date(ind_part),sgd.depth,sgd.o(:,ind_part),140:2.5:230,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
datetick('x','mm/dd'), xlim([min(dived.date(ind_part)) max(dived.date(ind_part))])
ylim([0 350])
hold on, contour(dived.date(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Date mm/dd 2015'),ylabel('Depth (m)')
caxis([192.5 225])
cb = colorbar, title(cb,'Oxygen (umol L-1)'), set(cb,'Fontsize',16);
subplot(7,1,4:5)
contourf(dived.date(ind_part),sgd.depth,sgd.s(:,ind_part),30:0.05:36,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
datetick('x','mm/dd'), xlim([min(dived.date(ind_part)) max(dived.date(ind_part))])
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.date(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Date mm/dd 2015'),ylabel('Depth (m)')
caxis([34.1 35.45])
cb = colorbar, title(cb,'Salinity'), set(cb,'Fontsize',16);
subplot(7,1,6:7)
contourf(dived.date(ind_part),sgd.depth,sgd.chl1(:,ind_part),-0.1:0.025:0.9,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
datetick('x','mm/dd'), xlim([min(dived.date(ind_part)) max(dived.date(ind_part))])
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.date(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Date mm/dd 2015'),ylabel('Depth (m)')
caxis([0 0.6])
cb = colorbar, title(cb,'Chlorophyll fluorescence'), set(cb,'Fontsize',16);
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