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
load ccar2015

% interpolate ssh on each glider dive
sshsg = interp3(ssh.lon_g,ssh.lat_g,ssh.date_g,ssh.ssh,dived.lon,dived.lat,dived.date);

% indeces of different transects
merid1 = dived.dive >= 50 & dived.dive <= 137 & dived.dive ~= 106; % first meridional transect
shortz1 = dived.dive >= 137 & dived.dive <= 151; % short zonal transect
shortm1 = dived.dive >= 151 & dived.dive <= 155; % short meridional transect
zonal1 = dived.dive >= 155 & dived.dive <= 226 &  dived.dive ~= 189; % first zonal transect
zonal2 = dived.dive >= 226 & dived.dive <= 265; % second zonal transect
merid2 = dived.dive >= 265 & dived.dive <= 318; % second meridional transect


%% Long zonal transect
ind_part = zonal1;
subplot(7,1,2:3)
contourf(dived.lon(ind_part),sgd.depth,sgd.o(:,ind_part),140:2.5:220,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
hold on, contour(dived.lon(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Longitude E'),ylabel('Depth (m)')
caxis([177.5 215])
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
cb = colorbar, title(cb,'Salinity'), set(cb,'Fontsize',16)
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
subplot(7,1,1)
plot(dived.lon(ind_part),sshsg(ind_part),'k'), xlabel('Longitude E'),ylabel('SSHA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
xlim([min(dived.lon(ind_part)) max(dived.lon(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)])

%% Long meridional transect
ind_part = merid1;
subplot(7,1,2:3)
contourf(dived.lat(ind_part),sgd.depth,sgd.o(:,ind_part),140:2.5:220,'edgecolor','none')
set(gca,'ydir','rev','Fontsize',16)
ylim([0 350])
hold on, contour(dived.lat(ind_part),sgd.depth,sgd.sig(:,ind_part),[23.5 24.3 25.3],'k--'), hold off
xlabel('Latitude N'),ylabel('Depth (m)')
caxis([177.5 215])
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
cb = colorbar, title(cb,'Salinity'), set(cb,'Fontsize',16)
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
subplot(7,1,1)
plot(dived.lat(ind_part),sshsg(ind_part),'k'), xlabel('Longitude E'),ylabel('SSHA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
xlim([min(dived.lat(ind_part)) max(dived.lat(ind_part))])
pp1 = get(gca,'Position'), set(gca,'Position',[pp(1) pp1(2) pp(3) pp1(4)])
