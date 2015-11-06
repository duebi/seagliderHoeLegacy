% dielSg146m11.m
%
% script to analyze diel cycles from 2015 SeaGlider mission sg146m11
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

%% Single Day fits of oxygen and bbp to estimate GOP/GPP

% Optimisation options
opts = optimset('Algorithm','trust-region-reflective','TolFun',1e-9,'TolX',1e-9,'MaxIter',40000,'MaxFunEval',20000);
% Days to perform the fit
days = fix(dived.date(1)):1:fix(dived.date(end));
srs = NaN(length(days),2);
mldsig_max = NaN*days;
Po = NaN*days; Ro = NaN*days;
Pbb = NaN*days; Rbb = NaN*days;
% Initialize variables for time serie in ml
dateo_ml = []; datebb_ml = []; opt_ml = []; bb_ml = []; srs_ml = [];

for i = 1:length(days)
        % isolate data from a single day
    ind_dd = fix(dived.date) == days(i);
    o_dd = sgd.opt(:,ind_dd); bb_dd = sgd.bbp660(:,ind_dd);
    sig_dd = sgd.sig(:,ind_dd);
    mld_dd = mld003(ind_dd);
    mldsig_dd = mld003sig(ind_dd);
    depth_dd = repmat(sgd.depth,1,sum(ind_dd));
    % Fs_dd = Fs(ind_dd);
    date_dd = dived.date(ind_dd);
    srs_temp = suncycle(mean(dived.lat(ind_dd)),mean(dived.lon(ind_dd)),days(i)); srs_temp = (srs_temp-10)/24; srs_temp(srs_temp<0) = srs_temp(srs_temp<0)+1;
    srs(i,:) = srs_temp; clear srs_temp
    x_fit = date_dd-days(i); 
    % Integration till the max daily density at the base of mixed layer
    mldsig_max(i) = max(mldsig_dd); % maximum daily density at the base of the mixed layer
    indsig_max = sig_dd<mldsig_max(i); % indeces of the values with density lower than mldsig_max(i)
    depth_dd(~indsig_max) = NaN;
    depthsig_max = nanmax(depth_dd);
    meandepthsig_max = mean(depthsig_max);
        % oxygen data for fit
    o_fit = o_dd; o_fit(~indsig_max) = NaN; % o_fitd(1:4,:) = NaN;
    ind_ofit = ~isnan(nanmean(o_fit));
    yo_fit = nanmean(o_fit(:,ind_ofit));%-cumtrapz(x_fit(ind_ofit),1000*86400*Fs_dd(ind_ofit)/meandepthsig_max);
    xo_fit = x_fit(ind_ofit);
        % backscattering data for fit
    bb_fit = bb_dd; bb_fit(~indsig_max) = NaN; %bb_fitd(1:5,:) = NaN;
    ybb_fit = nanmean(bb_fit);
    ind_bbfit = ~isnan(ybb_fit);
    ybb_fit = ybb_fit(ind_bbfit);
    xbb_fit = x_fit(ind_bbfit);
        % Save average mixed layer properties 
    opt_ml= [opt_ml yo_fit];
    bb_ml = [bb_ml ybb_fit];
    dateo_ml = [dateo_ml date_dd(ind_ofit)];
    datebb_ml = [datebb_ml date_dd(ind_bbfit)];
    % ------ Fitting theoretical model with lsqnonlin and constrainded parameters
    [tt,Pout,Rout] = diel_PR2(mean(dived.lat(ind_dd)),mean(dived.lon(ind_dd)),days(i),-10,1,0);
    tt3 = [tt(1:end-1)-1; tt; tt(2:end)+1]; Pout3 = [Pout(1:end-1); Pout; Pout(2:end)]; Rout3 = [Rout(1:end-1); Rout; Rout(2:end)];
    % ----- Oxygen (optode) -----------
    if length(yo_fit)>=4
        clf
        %costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),x_fit+param(4))-y_fitd;
        costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),xo_fit)-yo_fit;
        amp_fit_o = lsqnonlin(costfun,[1 1 1],[-Inf 0 0],[Inf Inf Inf],opts);
        Po(i) = amp_fit_o(2);
        Ro(i) = amp_fit_o(3);
        % PLOT
            subplot(2,2,1)
            contourf(xo_fit,sgd.depth,o_dd(:,ind_ofit),190:1:230,'edgecolor','none')
            hold on, l1 = plot(x_fit,depthsig_max,'k-'); hold off
            set(gca,'ydir','rev')
            xlim([0,1]),ylim([0 60])
            caxis([205 230])
            xlabel('Decimal day'), ylabel('Depth (m)')
            cb = colorbar('Location','NorthOutside','Fontsize',16);
            title(cb,'Oxygen (umol L-1)')
            subplot(2,2,3)
            patch([0 srs(i,1) srs(i,1) 0],[250 250 180 180],[0.9 0.9 0.9],'edgecolor','none')
            hold on, patch([srs(i,2) srs(i,1)+1 srs(i,1)+1 srs(i,2)],[250 250 180 180],[0.9 0.9 0.9],'edgecolor','none'), hold off
            hold on, p1 = plot(xo_fit,yo_fit,'o'); hold off
            hold on, l1 = plot(tt3,[amp_fit_o(1)+cumtrapz(tt3,amp_fit_o(2)*Pout3+amp_fit_o(3)*Rout3)],'k-'); hold off
            xlim([0 1]), ylim([min(yo_fit)-0.4 max(yo_fit)+0.4])
            legend([p1 l1],{'Datapoint','Fit'},'Fontsize',16)
            xlabel('Decimal day'), ylabel('Adjusted oxygen (umol L-1)')
            text(0.05,213,{['GOP = ' num2str(Po(i),2)];['R = ' num2str(Ro(i),2)]})
            title(datestr(days(i)))
    end
    % ----- Particle backscattering coefficient (bbp660) -----------
    if length(ybb_fit)>=4
        %costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),x_fit+param(4))-y_fitd;
        costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),xbb_fit)-ybb_fit;
        amp_fit_bb = lsqnonlin(costfun,[1 1 1],[-Inf 0 0],[Inf Inf Inf],opts);
        Pbb(i) = amp_fit_bb(2);
        Rbb(i) = amp_fit_bb(3);
        % PLOT
        subplot(2,2,2)
        contourf(xbb_fit,sgd.depth,bb_dd(:,ind_bbfit),0:1.25e-4:8e-3,'edgecolor','none')
        hold on, l1=plot(x_fit,depthsig_max,'k-'); hold off
        set(gca,'ydir','rev')
        xlim([0,1]),ylim([0 60])
        caxis([0 1e-3]) %caxis([1e-3 5e-3])
        xlabel('Decimal day'), ylabel('Depth (m)')
        cb = colorbar('Location','NorthOutside','Fontsize',16);
        title(cb,'bbp470 (m-1)')
        subplot(2,2,4)
        patch([0 srs(i,1) srs(i,1) 0],[0.08 0.08 0 0],[0.9 0.9 0.9],'edgecolor','none')
        hold on, patch([srs(i,2) srs(i,1)+1 srs(i,1)+1 srs(i,2)],[0.08 0.08 0 0],[0.9 0.9 0.9],'edgecolor','none'), hold off
        hold on, p1 = plot(xbb_fit,ybb_fit,'o'); hold off
        hold on, l1 = plot(tt3,[amp_fit_bb(1)+cumtrapz(tt3,amp_fit_bb(2)*Pout3+amp_fit_bb(3)*Rout3)],'k-'); hold off
        xlim([0 1]), ylim([min(ybb_fit)-0.0005 max(ybb_fit)+0.0005])
        legend([p1 l1],{'Datapoint','Fit'},'Fontsize',16)
        xlabel('Decimal day'), ylabel('Particle backscattering coeff. (m-1)')
        text(0.1,0.002,{['GCP = ' num2str(Pbb(i),2)];['R = ' num2str(Rbb(i),2)]})
        title(datestr(days(i)))        
    end
    %pause(0.1)
    clear ind_dd o_dd bb_dd sig_dd mld_dd mldsig_dd depth_dd date_dd
    clear o_fit bb_fit x_fit xo_fit yo_fit xbb_fit ybb_fit indsig_max depthsig_max
    clear ind_ofit ind_bbfit
    clear tt Pout Rout tt3 costfun 
end

%% Fit on variable time window
mldsig_max2 = NaN*days;
P2o = NaN*days; R2o = NaN*days;
P2bb = NaN*days; R2bb = NaN*days;
% Time before and after midnight to consider for the fit (to increase n)
deltaday = 0.5;

for i = 1:length(days)
        % isolate data from a single day
    ind_dd = dived.date < (days(i)+1+deltaday) & dived.date > (days(i)-deltaday);
    o_dd = sgd.opt(:,ind_dd);
    bb_dd = sgd.bbp660(:,ind_dd);
    sig_dd = sgd.sig(:,ind_dd);
    mld_dd = mld003(ind_dd);
    mldsig_dd = mld003sig(ind_dd);
    depth_dd = repmat(sgd.depth,1,sum(ind_dd));
    % Fs_dd = Fs(ind_dd);
    date_dd = dived.date(ind_dd);
    x_fit = date_dd-days(i); 
    % Integration till the max daily density at the base of mixed layer
    mldsig_max2(i) = max(mldsig_dd); % maximum daily density at the base of the mixed layer
    indsig_max = sig_dd<mldsig_max2(i); % indeces of the values with density lower than mldsig_max2(i)
    depth_dd(~indsig_max) = NaN;
    depthsig_max = nanmax(depth_dd);
    meandepthsig_max = mean(depthsig_max);
        % oxygen data for fit
    o_fit = o_dd; o_fit(~indsig_max) = NaN; % o_fitd(1:4,:) = NaN;
    ind_ofit = ~isnan(nanmean(o_fit));
    yo_fit = nanmean(o_fit(:,ind_ofit));%-cumtrapz(x_fit(ind_ofit),1000*86400*Fs_dd(ind_ofit)/meandepthsig_max);
    xo_fit = x_fit(ind_ofit);
            % backscattering data for fit
    bb_fit = bb_dd; bb_fit(~indsig_max) = NaN; %bb_fitd(1:5,:) = NaN;
    ybb_fit = nanmean(bb_fit);
    ind_bbfit = ~isnan(ybb_fit);
    ybb_fit = ybb_fit(ind_bbfit);
    xbb_fit = x_fit(ind_bbfit);
    % ------ Fitting theoretical model with lsqnonlin and constrainded parameters
    [tt,Pout,Rout] = diel_PR2(mean(dived.lat(ind_dd)),mean(dived.lon(ind_dd)),days(i),-10,1,0);
    tt3 = [tt(1:end-1)-1; tt; tt(2:end)+1]; Pout3 = [Pout(1:end-1); Pout; Pout(2:end)]; Rout3 = [Rout(1:end-1); Rout; Rout(2:end)];
    % ----- Oxygen (optode) -----------
    if length(yo_fit)>=4
        clf
        %costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),x_fit+param(4))-y_fitd;
        costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),xo_fit)-yo_fit;
        amp_fit_o = lsqnonlin(costfun,[1 1 1],[-Inf 0 0],[Inf Inf Inf],opts);
        P2o(i) = amp_fit_o(2);
        R2o(i) = amp_fit_o(3);
        % PLOT
            subplot(2,2,1)
            contourf(xo_fit,sgd.depth,o_dd(:,ind_ofit),190:1:230,'edgecolor','none')
            hold on, l1 = plot(x_fit,depthsig_max,'k-'); hold off
            set(gca,'ydir','rev')
            xlim([0-deltaday,1+deltaday]),ylim([0 60])
            caxis([205 230])
            xlabel('Decimal day'), ylabel('Depth (m)')
            cb = colorbar('Location','NorthOutside','Fontsize',16);
            title(cb,'Oxygen (umol L-1)')
            subplot(2,2,3)
            patch([srs(i,2)-1 srs(i,1) srs(i,1) srs(i,2)-1],[250 250 180 180],[0.9 0.9 0.9],'edgecolor','none')
            hold on, patch([srs(i,2) srs(i,1)+1 srs(i,1)+1 srs(i,2)],[250 250 180 180],[0.9 0.9 0.9],'edgecolor','none'), hold off
            hold on, p1 = plot(xo_fit,yo_fit,'o'); hold off
            hold on, l1 = plot(tt3,[amp_fit_o(1)+cumtrapz(tt3,amp_fit_o(2)*Pout3+amp_fit_o(3)*Rout3)],'k-'); hold off
            xlim([0-deltaday,1+deltaday]), ylim([min(yo_fit)-0.4 max(yo_fit)+0.4])
            legend([p1 l1],{'Datapoint','Fit'},'Fontsize',16)
            xlabel('Decimal day'), ylabel('Adjusted oxygen (umol L-1)')
            text(0.05,213,{['GOP = ' num2str(Po(i),2)];['R = ' num2str(Ro(i),2)]})
            title(datestr(days(i)))
    end
    % ----- Particle backscattering coefficient (bbp660) -----------
    if length(ybb_fit)>=4
        %costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),x_fit+param(4))-y_fitd;
        costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),xbb_fit)-ybb_fit;
        amp_fit_bb = lsqnonlin(costfun,[1 1 1],[-Inf 0 0],[Inf Inf Inf],opts);
        P2bb(i) = amp_fit_bb(2);
        R2bb(i) = amp_fit_bb(3);
        % PLOT
        subplot(2,2,2)
        contourf(xbb_fit,sgd.depth,bb_dd(:,ind_bbfit),0:1.25e-4:8e-3,'edgecolor','none')
        hold on, l1=plot(x_fit,depthsig_max,'k-'); hold off
        set(gca,'ydir','rev')
        xlim([0-deltaday,1+deltaday]),ylim([0 60])
        caxis([0 1e-3]) % caxis([1e-3 5e-3]) 
        xlabel('Decimal day'), ylabel('Depth (m)')
        cb = colorbar('Location','NorthOutside','Fontsize',16);
        title(cb,'bbp660 (m-1)')
        subplot(2,2,4)
        patch([srs(i,2)-1 srs(i,1) srs(i,1) srs(i,2)-1],[0.08 0.08 0 0],[0.9 0.9 0.9],'edgecolor','none')
        hold on, patch([srs(i,2) srs(i,1)+1 srs(i,1)+1 srs(i,2)],[0.08 0.08 0 0],[0.9 0.9 0.9],'edgecolor','none'), hold off
        hold on, p1 = plot(xbb_fit,ybb_fit,'o'); hold off
        hold on, l1 = plot(tt3,[amp_fit_bb(1)+cumtrapz(tt3,amp_fit_bb(2)*Pout3+amp_fit_bb(3)*Rout3)],'k-'); hold off
        xlim([0-deltaday,1+deltaday]), ylim([min(ybb_fit)-0.0001 max(ybb_fit)+0.0001])
        legend([p1 l1],{'Datapoint','Fit'},'Fontsize',16)
        xlabel('Decimal day'), ylabel('Particle backscattering coeff. (m-1)')
        text(0.1,0.002,{['GCP = ' num2str(Pbb(i),2)];['R = ' num2str(Rbb(i),2)]})
        title(datestr(days(i)))        
    end
    pause(0.1)
    clear ind_dd o_dd bb_dd sig_dd mld_dd mldsig_dd depth_dd date_dd
    clear o_fit bb_fit x_fit xo_fit yo_fit xbb_fit ybb_fit indsig_max depthsig_max
    clear ind_ofit ind_bbfit
    clear tt Pout Rout tt3 costfun 
end

%% Plot mixed layer properties and rate estimates from oxygen

clf
subplot(4,1,1)
plot(dived.date,sshsg,'k--',dived.date,slasg*100,'k-'),
ylabel('SSHA/SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)])
lg = legend('ccar','aviso'), set(lg,'Fontsize',16), legend('boxoff')
subplot(4,1,2)
hold on
for i = 2:length(days)
    patch(days(i)+[srs(i-1,2)-1 srs(i-1,2)-1 srs(i,1) srs(i,1)],[0 300 300 0],[0.9 0.9 0.9],'edgecolor','none')
end
plot(dateo_ml,opt_ml,'color',[0 0 0])
hold off
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)])
ylim([nanmin(opt_ml)-0.25 nanmax(opt_ml)+0.25])
set(gca,'Fontsize',16) 
ylabel('Average ML O2 (mmol m-3)')
subplot(4,1,3)
bar(days+0.5,Po,'k')
hold on, bar(days+0.5,-Ro,'r'), hold off
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)])
ylim([-4.5 4.5])
legend('GOP','R')
legend('boxoff')
ylabel('O2 rates (mmol O2 m-3 d-1)')
set(gca,'Fontsize',16)
subplot(4,1,4)
bar(days+0.5,P2o,'k')
hold on, bar(days+0.5,-R2o,'r'), hold off
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)])
ylim([-4.5 4.5])
legend('GOP','R')
legend('boxoff')
ylabel('O2 rates (mmol O2 m-3 d-1)')
xlabel('mm/dd 2015')
set(gca,'Fontsize',16)

%% Plot mixed layer properties and rate estimates from bbp

clf
subplot(4,1,1)
plot(dived.date,sshsg,'k--',dived.date,slasg*100,'k-'),
ylabel('SSHA/SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)])
lg = legend('ccar','aviso'), set(lg,'Fontsize',16), legend('boxoff')
subplot(4,1,2)
hold on
for i = 2:length(days)
    patch(days(i)+[srs(i-1,2)-1 srs(i-1,2)-1 srs(i,1) srs(i,1)],[0 300 300 0],[0.9 0.9 0.9],'edgecolor','none')
end
plot(datebb_ml,bb_ml,'color',[0 0 0])
hold off
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([datebb_ml(1) datebb_ml(end)])
ylim([nanmin(bb_ml)-0.0001 nanmax(bb_ml)+0.0001])
ylim([1e-4 8e-4])
set(gca,'Fontsize',16) 
ylabel('Average ML bbp660 (m-1)')
subplot(4,1,3)
bar(days+0.5,Pbb,'b')
hold on, bar(days+0.5,-Rbb,'g'), hold off
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([datebb_ml(1) datebb_ml(end)])
ylim([-5e-4 5e-4])
legend('GPP','R')
legend('boxoff')
ylabel('bbp rates (m-1 d-1)')
set(gca,'Fontsize',16)
subplot(4,1,4)
bar(days+0.5,P2bb,'b')
hold on, bar(days+0.5,-R2bb,'g'), hold off
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([datebb_ml(1) datebb_ml(end)])
ylim([-5e-4 5e-4])
legend('GPP','R')
legend('boxoff')
ylabel('bbp rates (m-1 d-1)')
xlabel('mm/dd 2015')
set(gca,'Fontsize',16)

%% Whole cruise diel averages in the mixed layer
clear bins diel_opt diel_bbp660 md_opt md_bbp660 anom_opt anom_bbp660

binwidth = 1/12; %1/24;
binrange = [binwidth/2 1-binwidth/2];

[diel_opt,bins] = dielbinning(opt_ml,dateo_ml-fix(dateo_ml),binwidth,binrange,'mean');
[diel_bbp660,bins] = dielbinning(bb_ml,datebb_ml-fix(datebb_ml),binwidth,binrange,'mean');

subplot(1,2,1)
plot([bins(end)-1 bins bins(1)+1],[diel_opt(end) diel_opt diel_opt(1)],'ko--')
xlim([0 1])
subplot(1,2,2)
plot([bins(end)-1 bins bins(1)+1],[diel_opt(end) diel_opt diel_opt(1)],'ko--')