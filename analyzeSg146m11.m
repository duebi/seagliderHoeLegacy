% analyzeSg146m11.m
%
% script to analyze data from 2015 SeaGlider mission sg146m11
% 
% Benedetto Barone - Oct 2015

mission = 'sg146_m11';
upth = userpath; 
sgpath =  [upth(1:end-1) '/Data/seaglider/' mission];
cd(sgpath);
clear upth sgpath mission

load sg146m11data
load ../ccar2015
load ../aviso2015
load ../../colorbrewer_anom
load ../hawaii.dat
coastline = hawaii;
coastline(:,1) = -coastline(:,1);
ind_n = find(isnan(coastline(:,1)));
clear hawaii

% Horizontal length of the arrays
nd = length(dived.dive);

% interpolate ssh on each glider dive
sshsg = interp3(ssh.lon_g,ssh.lat_g,ssh.date_g,ssh.ssh,dived.lon,dived.lat,dived.date);
slasg = interp3(sla.lon_g,sla.lat_g,sla.date_g,sla.sla,dived.lon,dived.lat,dived.date);

% indeces of different transects
merid1 = dived.dive >= 50 & dived.dive <= 137 & dived.dive ~= 106; % first meridional transect
shortz1 = dived.dive >= 137 & dived.dive <= 151; % short zonal transect
shortm1 = dived.dive >= 151 & dived.dive <= 155; % short meridional transect
zonal1 = dived.dive >= 155 & dived.dive <= 226 &  dived.dive ~= 189; % first zonal transect
zonal2 = dived.dive >= 226 & dived.dive <= 265; % second zonal transect
merid2 = dived.dive >= 265 & dived.dive <= 318; % second meridional transect

% Find mixed layer depth and density at the base of the mixed layer (0.03 kg m-3 difference from sigma at 4 m depth)
mld = sgd.sig - repmat(sgd.sig(2,:),height(sgd),1) - 0.03;
mld(mld<0) = NaN; mld(1,:) = NaN;
[sig003,ind003] = nanmin(mld);
mld003 = sgd.depth(ind003);
mld003sig = sig003 + sgd.sig(2,:) + 0.03;

%% Long zonal transect
ind_part = zonal1; %shortz1;
subplot(7,1,2:3)
contourf(dived.lon(ind_part),sgd.depth,sgd.o(:,ind_part),140:2.5:230,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
hold on, contour(dived.lon(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Longitude E'),ylabel('Depth (m)')
caxis([190 224])
cb = colorbar, title(cb,'Oxygen (umol L-1)'), set(cb,'Fontsize',16)
subplot(7,1,4:5)
contourf(dived.lon(ind_part),sgd.depth,sgd.s(:,ind_part),30:0.05:36,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lon(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Longitude E'),ylabel('Depth (m)')
caxis([34.1 35.45])
cb = colorbar, title(cb,'Salinity (g kg-1)'), set(cb,'Fontsize',16)
subplot(7,1,6:7)
contourf(dived.lon(ind_part),sgd.depth,sgd.chl2(:,ind_part),0:0.25:10,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lon(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Longitude E'),ylabel('Depth (m)')
caxis([0 7.5])
cb = colorbar, title(cb,'Chl-a fluorescence'), set(cb,'Fontsize',16)
pp = get(gca,'Position')
%{
subplot(7,1,6:7)
contourf(dived.lon(ind_part),sgd.depth,sgd.bbp470(:,ind_part),0:5e-4:1e-2,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lon(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Longitude E'),ylabel('Depth (m)')
caxis([5e-4 4e-3])
cb = colorbar, title(cb,'bbp 470 nm (m-1)'), set(cb,'Fontsize',16)
pp = get(gca,'Position')
%}
subplot(7,1,1)
plot(dived.lon(ind_part),sshsg(ind_part),'k--',dived.lon(ind_part),slasg(ind_part)*100,'k-'), xlabel('Longitude E'),ylabel('SSHA/SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
xlim([min(dived.lon(ind_part)) max(dived.lon(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)])
lg = legend('ccar','aviso'), set(lg,'Fontsize',16), legend('boxoff')

%% Long meridional transect
ind_part = merid1;
subplot(7,1,2:3)
contourf(dived.lat(ind_part),sgd.depth,sgd.o(:,ind_part),140:2.5:230,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([192 227])
cb = colorbar, title(cb,'Oxygen (umol L-1)'), set(cb,'Fontsize',16)
subplot(7,1,4:5)
contourf(dived.lat(ind_part),sgd.depth,sgd.s(:,ind_part),30:0.05:36,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([34.1 35.45])
cb = colorbar, title(cb,'Salinity (g kg-1)'), set(cb,'Fontsize',16)
subplot(7,1,6:7)
contourf(dived.lat(ind_part),sgd.depth,sgd.chl2(:,ind_part),0:0.25:10,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
%hold on, plot(mdate(1:end_1),mld003(1:end_1),'k'), hold off
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
%datetick('x','mm/dd')
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([0 7.5])
cb = colorbar, title(cb,'chl (m-1)'), set(cb,'Fontsize',16)
pp = get(gca,'Position')
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
cb = colorbar, title(cb,'bbp 470 nm (m-1)'), set(cb,'Fontsize',16)
pp = get(gca,'Position')
%}
subplot(7,1,1)
plot(dived.lat(ind_part),sshsg(ind_part),'k--',dived.lat(ind_part),slasg(ind_part)*100,'k-'), xlabel('Latitude E'),ylabel('SSHA/SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
xlim([min(dived.lat(ind_part)) max(dived.lat(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)])
lg = legend('ccar','aviso'), set(lg,'Fontsize',16), legend('boxoff')

%% Isopycnal zonal transect of anomalies

ind_part = zonal1; %shortz1;
subplot(7,1,2:3)
contourf(dived.lon(ind_part),isod.depth,isod.o(:,ind_part)-repmat(nanmedian(isod.o,2),1,sum(ind_part)),-40:1:40,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
xlabel('Longitude E'),ylabel('Average isopycnal depth (m)')
caxis([-18 18])
cb = colorbar, title(cb,'Oxygen anomaly (umol L-1)'), set(cb,'Fontsize',16)
subplot(7,1,4:5)
contourf(dived.lon(ind_part),isod.depth,isod.s(:,ind_part)-repmat(nanmedian(isod.s,2),1,sum(ind_part)),-0.6:0.01:0.6,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
xlabel('Longitude E'),ylabel('Average isopycnal depth (m)')
caxis([-0.15 0.15])
cb = colorbar, title(cb,'Salinity anomaly (g kg-1)'), set(cb,'Fontsize',16)
%{
subplot(7,1,6:7)
contourf(dived.lon(ind_part),isod.depth,isod.chl2(:,ind_part)-repmat(nanmedian(isod.chl2,2),1,sum(ind_part)),-4:0.1:4,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
xlabel('Longitude E'),ylabel('Average isopycnal depth (m)')
caxis([-3 3])
cb = colorbar, title(cb,'Chl-a anomaly'), set(cb,'Fontsize',16)
pp = get(gca,'Position')
%}
subplot(7,1,6:7)
contourf(dived.lon(ind_part),isod.depth,isod.bbp470(:,ind_part)-repmat(nanmedian(isod.bbp470,2),1,sum(ind_part)),-0.005:5e-5:0.005,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
xlabel('Longitude E'),ylabel('Average isopycnal depth (m)')
caxis([-0.0012 0.0012])
cb = colorbar, title(cb,'bbp 470 anomaly (m-1)'), set(cb,'Fontsize',16)
pp = get(gca,'Position')
subplot(7,1,1)
plot(dived.lon(ind_part),sshsg(ind_part),'k--',dived.lon(ind_part),slasg(ind_part)*100,'k-'), xlabel('Longitude E'),ylabel('SSHA/SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
xlim([min(dived.lon(ind_part)) max(dived.lon(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)])
lg = legend('ccar','aviso'), set(lg,'Fontsize',16), legend('boxoff')
colormap(anom_map2)

%% Isopycnal meridional transect of anomalies

clf
ind_part = merid1; %shortz1;
subplot(7,1,2:3)
contourf(dived.lat(ind_part),isod.depth,isod.o(:,ind_part)-repmat(nanmedian(isod.o,2),1,sum(ind_part)),-40:1:40,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
xlabel('Latitude N'),ylabel('Average isopycnal depth (m)')
caxis([-18 18])
cb = colorbar, title(cb,'Oxygen anomaly (umol L-1)'), set(cb,'Fontsize',16)
subplot(7,1,4:5)
contourf(dived.lat(ind_part),isod.depth,isod.s(:,ind_part)-repmat(nanmedian(isod.s,2),1,sum(ind_part)),-0.6:0.01:0.6,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
xlabel('Latitude N'),ylabel('Average isopycnal depth (m)')
caxis([-0.15 0.15])
cb = colorbar, title(cb,'Salinity anomaly (g kg-1)'), set(cb,'Fontsize',16)
%{
subplot(7,1,6:7)
contourf(dived.lat(ind_part),isod.depth,isod.chl2(:,ind_part)-repmat(nanmedian(isod.chl2,2),1,sum(ind_part)),-4:0.1:4,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
xlabel('Latitude N'),ylabel('Average isopycnal depth (m)')
caxis([-3 3])
cb = colorbar, title(cb,'Chl-a anomaly'), set(cb,'Fontsize',16)
pp = get(gca,'Position')
%}
subplot(7,1,6:7)
contourf(dived.lat(ind_part),isod.depth,isod.bbp470(:,ind_part)-repmat(nanmedian(isod.bbp470,2),1,sum(ind_part)),-0.005:5e-5:0.005,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
xlabel('Latitude N'),ylabel('Average isopycnal depth (m)')
caxis([-0.0012 0.0012])
cb = colorbar, title(cb,'bbp 470 anomaly (m-1)'), set(cb,'Fontsize',16)
pp = get(gca,'Position')
subplot(7,1,1)
plot(dived.lat(ind_part),sshsg(ind_part),'k--',dived.lat(ind_part),slasg(ind_part)*100,'k-'), xlabel('Latitude N'),ylabel('SSHA/SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
xlim([min(dived.lat(ind_part)) max(dived.lat(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)])
lg = legend('ccar','aviso'), set(lg,'Fontsize',16), legend('boxoff')
colormap(anom_map2)

%% Map of transects

% indeces for the central dates of the transects
ind_aug24 = sla.date == datenum(2015,8,24);
ind_aug2 = sla.date == datenum(2015,8,2);

subplot(1,2,1)
contourf(sla.lon,sla.lat,sla.sla(:,:,ind_aug2)*100,-24:4:24,'edgecolor','none')
title('Meridional transect - July 24 to Aug 12')
colormap(anom_map2)
cb = colorbar; set(cb,'Fontsize',16), title(cb,'SLA for Aug 2 (cm)')
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
    plot(dived.lon,dived.lat,'w-')
    plot(dived.lon(merid1),dived.lat(merid1),'k--','LineWidth',3)
hold off
set(gca,'Fontsize',16)
axis equal
xlim([-161 -154]),ylim([19 26])

subplot(1,2,2)
contourf(sla.lon,sla.lat,sla.sla(:,:,ind_aug24)*100,-24:4:24,'edgecolor','none')
title('Meridional transect - Aug 16 to Sep 1')
colormap(anom_map2)
cb = colorbar; set(cb,'Fontsize',16), title(cb,'SLA for Aug 24 (cm)')
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
    plot(dived.lon,dived.lat,'w-')
    plot(dived.lon(zonal1),dived.lat(zonal1),'k--','LineWidth',3)
hold off
set(gca,'Fontsize',16)
axis equal
xlim([-161 -154]),ylim([19 26])