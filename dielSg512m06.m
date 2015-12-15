% dielSg512m06.m
%
% script to analyze diel cycles from 2015 SeaGlider mission sg512m06
% 
% Benedetto Barone - Nov 2015

mission = 'sg512_m06';
upth = userpath; 
sgpath =  [upth(1:end-1) '/Data/seaglider/' mission];
cd(sgpath);
clear upth sgpath mission

load sg512m06data

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

% Optimisation options
opts = optimset('Algorithm','trust-region-reflective','TolFun',1e-9,'TolX',1e-9,'MaxIter',40000,'MaxFunEval',20000);
% Days to perform the fit + sunrise and sunset time
days = fix(dived.date(1)):1:fix(dived.date(end));
srs = NaN(length(days),2);
for i = 1:length(days)
        % isolate data from a single day
    ind_dd = fix(dived.date) == days(i);
    srs_temp = suncycle(mean(dived.lat(ind_dd)),mean(dived.lon(ind_dd)),days(i)); srs_temp = (srs_temp-10)/24; srs_temp(srs_temp<0) = srs_temp(srs_temp<0)+1;
    srs(i,:) = srs_temp;
    clear ind_dd srs_temp
end

% Interpolate missing values from air sea flux

dived.Fs(isnan(dived.Fs)) = interp1(dived.date(~isnan(dived.Fs)),dived.Fs(~isnan(dived.Fs)),dived.date(isnan(dived.Fs)));

%% Single Day fits of oxygen and bbp to estimate GOP/GPP & R

mldsig_max = NaN*days;
Po = NaN*days; Ro = NaN*days;
rsq_o = NaN*days; pval_o = NaN*days;
Pbb = NaN*days; Rbb = NaN*days;
rsq_bb = NaN*days; pval_bb = NaN*days;
% Initialize variables for time serie in ml
dateo_ml = []; datebb_ml = []; opt_ml = []; bb_ml = []; srs_ml = [];  chl_ml = NaN*days; 
% Average daily mixed layer characteristics
daychl_ml = NaN*days; dayopt_ml = NaN*days; daybb_ml = NaN*days;

for i = 1:length(days)
        % isolate data from a single day
    ind_dd = fix(dived.date) == days(i);
    daylat(i) = nanmean(dived.lat(ind_dd));
    daylon(i) = nanmean(dived.lon(ind_dd));
    o_dd = sgd.opt(:,ind_dd); bb_dd = sgd.bbp650(:,ind_dd);
    chl_dd = sgd.chl1(:,ind_dd);
    sig_dd = sgd.sig(:,ind_dd);
    mld_dd = mld003(ind_dd);
    mldsig_dd = mld003sig(ind_dd);
    depth_dd = repmat(sgd.depth,1,sum(ind_dd));
    Fs_dd = dived.Fs(ind_dd);
    date_dd = dived.date(ind_dd);
    x_fit = date_dd-days(i); 
    % Integration till the max daily density at the base of mixed layer
    mldsig_max(i) = max(mldsig_dd); % maximum daily density at the base of the mixed layer
    indsig_max = sig_dd<mldsig_max(i); % indeces of the values with density lower than mldsig_max(i)
    depth_dd(~indsig_max) = NaN;
    depthsig_max = nanmax(depth_dd);
    meandepthsig_max = nanmean(depthsig_max);
    % daily averages of mixed layer chl, oxyg. and bbp
    daychl_ml(i) = nanmedian(nanmedian(chl_dd(indsig_max))); dayopt_ml(i) = nanmedian(nanmedian(o_dd(indsig_max))); daybb_ml(i) = nanmedian(nanmedian(bb_dd(indsig_max)));
        % oxygen data for fit
    o_fit = o_dd; o_fit(~indsig_max) = NaN; % o_fitd(1:4,:) = NaN;
    ind_ofit = ~isnan(nanmedian(o_fit));    
    xo_fit = x_fit(ind_ofit);
        % backscattering data for fit
    bb_fit = bb_dd; bb_fit(~indsig_max) = NaN; %bb_fitd(1:5,:) = NaN;
    ybb_fit = nanmedian(bb_fit);
    ind_bbfit = ~isnan(ybb_fit);
    ybb_fit = ybb_fit(ind_bbfit);
    xbb_fit = x_fit(ind_bbfit);
        % Save average mixed layer properties 
    opt_ml= [opt_ml nanmedian(o_fit(:,ind_ofit))];
    bb_ml = [bb_ml ybb_fit];
    dateo_ml = [dateo_ml date_dd(ind_ofit)];
    datebb_ml = [datebb_ml date_dd(ind_bbfit)];
    % ------ Fitting theoretical model with lsqnonlin and constrainded parameters
    [tt,Pout,Rout] = diel_PR2(mean(dived.lat(ind_dd)),mean(dived.lon(ind_dd)),days(i),-10,1,0);
    tt3 = [tt(1:end-1)-1; tt; tt(2:end)+1]; Pout3 = [Pout(1:end-1); Pout; Pout(2:end)]; Rout3 = [Rout(1:end-1); Rout; Rout(2:end)];
    % ----- Oxygen (optode) -----------
    if sum(ind_ofit)>=4
        yo_fit = nanmedian(o_fit(:,ind_ofit))-cumtrapz(x_fit(ind_ofit),1000*Fs_dd(ind_ofit)/meandepthsig_max);
        clf
        %costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),x_fit+param(4))-y_fitd;
        costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),xo_fit)-yo_fit;
        amp_fit_o = lsqnonlin(costfun,[1 1 1],[-Inf 0 0],[Inf Inf Inf],opts);
        yo_model = interp1(tt3,amp_fit_o(1)+cumtrapz(tt3,amp_fit_o(2)*Pout3+amp_fit_o(3)*Rout3),xo_fit);
        [r_temp, p_temp] = corrcoef(yo_fit,yo_model);
        rsq_o(i) = r_temp(2)^2; pval_o(i) = p_temp(2);
        clear r_temp p_temp
        Po(i) = amp_fit_o(2);
        Ro(i) = amp_fit_o(3);
        % PLOT
            subplot(2,2,1)
            contourf(xo_fit,sgd.depth,o_dd(:,ind_ofit),190:1:230,'edgecolor','none')
            hold on, l1 = plot(x_fit,depthsig_max,'k-'); hold off
            set(gca,'ydir','rev')
            xlim([0,1]),ylim([0 80])
            caxis([205 230])
            xlabel('Decimal day'), ylabel('Depth (m)')
            cb = colorbar('Location','NorthOutside','Fontsize',16);
            title(cb,'Oxygen (umol L-1)')
            subplot(2,2,3)
            patch([0 srs(i,1) srs(i,1) 0],[250 250 180 180],[0.9 0.9 0.9],'edgecolor','none')
            hold on, patch([srs(i,2) srs(i,1)+1 srs(i,1)+1 srs(i,2)],[250 250 180 180],[0.9 0.9 0.9],'edgecolor','none'), hold off
            hold on, p1 = plot(xo_fit,yo_fit,'o'); hold off
            hold on, p2 = plot(xo_fit,nanmean(o_fit(:,ind_ofit)),'ro'); hold off
            hold on, l1 = plot(tt3,[amp_fit_o(1)+cumtrapz(tt3,amp_fit_o(2)*Pout3+amp_fit_o(3)*Rout3)],'k-'); hold off
            xlim([0 1]), ylim([min(yo_fit)-0.4 max(yo_fit)+0.4])
            legend([p1 l1],{'Datapoint','Fit'},'Fontsize',16)
            xlabel('Decimal day'), ylabel('Adjusted oxygen (umol L-1)')
            text(0.05,213,{['GOP = ' num2str(Po(i),2)];['R = ' num2str(Ro(i),2)]})
            title(datestr(days(i)))
    end
    % ----- Particle backscattering coefficient (bbp650) -----------
    if length(ybb_fit)>=4
        %costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),x_fit+param(4))-y_fitd;
        costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),xbb_fit)-ybb_fit;
        amp_fit_bb = lsqnonlin(costfun,[1 1 1],[-Inf 0 0],[Inf Inf Inf],opts);
        ybb_model = interp1(tt3,amp_fit_bb(1)+cumtrapz(tt3,amp_fit_bb(2)*Pout3+amp_fit_bb(3)*Rout3),xbb_fit);
        [r_temp, p_temp] = corrcoef(ybb_fit,ybb_model);
        rsq_bb(i) = r_temp(2)^2; pval_bb(i) = p_temp(2);
        Pbb(i) = amp_fit_bb(2);
        Rbb(i) = amp_fit_bb(3);
        % PLOT
        subplot(2,2,2)
        contourf(xbb_fit,sgd.depth,bb_dd(:,ind_bbfit),0:1.25e-4:8e-3,'edgecolor','none')
        hold on, l1=plot(x_fit,depthsig_max,'k-'); hold off
        set(gca,'ydir','rev')
        xlim([0,1]),ylim([0 80])
        caxis([0 1e-3]) %caxis([1e-3 5e-3])
        xlabel('Decimal day'), ylabel('Depth (m)')
        cb = colorbar('Location','NorthOutside','Fontsize',16);
        title(cb,'bbp470 (m-1)')
        subplot(2,2,4)
        patch([0 srs(i,1) srs(i,1) 0],[0.08 0.08 0 0],[0.9 0.9 0.9],'edgecolor','none')
        hold on, patch([srs(i,2) srs(i,1)+1 srs(i,1)+1 srs(i,2)],[0.08 0.08 0 0],[0.9 0.9 0.9],'edgecolor','none'), hold off
        hold on, p1 = plot(xbb_fit,ybb_fit,'o'); hold off
        hold on, l1 = plot(tt3,[amp_fit_bb(1)+cumtrapz(tt3,amp_fit_bb(2)*Pout3+amp_fit_bb(3)*Rout3)],'k-'); hold off
        xlim([0 1]), ylim([min(ybb_fit)-0.00001 max(ybb_fit)+0.00001])
        legend([p1 l1],{'Datapoint','Fit'},'Fontsize',16)
        xlabel('Decimal day'), ylabel('Particle backscattering coeff. (m-1)')
        text(0.1,0.002,{['GCP = ' num2str(Pbb(i),2)];['R = ' num2str(Rbb(i),2)]})
        title(datestr(days(i)))        
    end
    %pause
    clear ind_dd o_dd bb_dd sig_dd mld_dd mldsig_dd depth_dd date_dd Fs_dd
    clear o_fit bb_fit x_fit xo_fit yo_fit xbb_fit ybb_fit indsig_max depthsig_max
    clear ind_ofit ind_bbfit
    clear tt Pout Rout tt3 costfun 
end

%% Fit on variable time window
mldsig_max2 = NaN*days;
P2o = NaN*days; R2o = NaN*days;
rsq2_o = NaN*days; pval2_o = NaN*days;
P2bb = NaN*days; R2bb = NaN*days;
rsq2_bb = NaN*days; pval2_bb = NaN*days;

% Time before and after midnight to consider for the fit (to increase n)
deltaday = 0.5;

for i = 1:length(days)
        % isolate data from a single day
    ind_dd = dived.date < (days(i)+1+deltaday) & dived.date > (days(i)-deltaday);
    o_dd = sgd.opt(:,ind_dd);
    bb_dd = sgd.bbp650(:,ind_dd);
    sig_dd = sgd.sig(:,ind_dd);
    mld_dd = mld003(ind_dd);
    mldsig_dd = mld003sig(ind_dd);
    depth_dd = repmat(sgd.depth,1,sum(ind_dd));
    Fs_dd = dived.Fs(ind_dd);
    date_dd = dived.date(ind_dd);
    x_fit = date_dd-days(i); 
    % Integration till the max daily density at the base of mixed layer
    mldsig_max2(i) = max(mldsig_dd); % maximum daily density at the base of the mixed layer
    indsig_max = sig_dd<mldsig_max2(i); % indeces of the values with density lower than mldsig_max2(i)
    depth_dd(~indsig_max) = NaN;
    depthsig_max = nanmax(depth_dd);
    meandepthsig_max = nanmean(depthsig_max);
        % oxygen data for fit
    o_fit = o_dd; o_fit(~indsig_max) = NaN; % o_fitd(1:4,:) = NaN;
    ind_ofit = ~isnan(nanmedian(o_fit));
    xo_fit = x_fit(ind_ofit);
            % backscattering data for fit
    bb_fit = bb_dd; bb_fit(~indsig_max) = NaN; %bb_fitd(1:5,:) = NaN;
    ybb_fit = nanmedian(bb_fit);
    ind_bbfit = ~isnan(ybb_fit);
    ybb_fit = ybb_fit(ind_bbfit);
    xbb_fit = x_fit(ind_bbfit);
    % ------ Fitting theoretical model with lsqnonlin and constrainded parameters
    [tt,Pout,Rout] = diel_PR2(mean(dived.lat(ind_dd)),mean(dived.lon(ind_dd)),days(i),-10,1,0);
    tt3 = [tt(1:end-1)-1; tt; tt(2:end)+1]; Pout3 = [Pout(1:end-1); Pout; Pout(2:end)]; Rout3 = [Rout(1:end-1); Rout; Rout(2:end)];
    % ----- Oxygen (optode) -----------
    if sum(ind_ofit)>=4
        yo_fit = nanmedian(o_fit(:,ind_ofit))-cumtrapz(x_fit(ind_ofit),1000*Fs_dd(ind_ofit)/meandepthsig_max);
        clf
        %costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),x_fit+param(4))-y_fitd;
        costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),xo_fit)-yo_fit;
        amp_fit_o = lsqnonlin(costfun,[1 1 1],[-Inf 0 0],[Inf Inf Inf],opts);
        yo_model = interp1(tt3,amp_fit_o(1)+cumtrapz(tt3,amp_fit_o(2)*Pout3+amp_fit_o(3)*Rout3),xo_fit);
        [r_temp, p_temp] = corrcoef(yo_fit,yo_model);
        rsq2_o(i) = r_temp(2)^2; pval2_o(i) = p_temp(2);
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
            hold on, p2 = plot(xo_fit,nanmean(o_fit(:,ind_ofit)),'ro'); hold off
            hold on, l1 = plot(tt3,[amp_fit_o(1)+cumtrapz(tt3,amp_fit_o(2)*Pout3+amp_fit_o(3)*Rout3)],'k-'); hold off
            xlim([0-deltaday,1+deltaday]), ylim([min(yo_fit)-0.4 max(yo_fit)+0.4])
            legend([p1 l1],{'Datapoint','Fit'},'Fontsize',16)
            xlabel('Decimal day'), ylabel('Adjusted oxygen (umol L-1)')
            text(0.05,213,{['GOP = ' num2str(P2o(i),2)];['R = ' num2str(R2o(i),2)]})
            title(datestr(days(i)))
    end
    % ----- Particle backscattering coefficient (bbp650) -----------
    if length(ybb_fit)>=4
        %costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),x_fit+param(4))-y_fitd;
        costfun = @(param) interp1(tt3,param(1)+cumtrapz(tt3,param(2)*Pout3+param(3)*Rout3),xbb_fit)-ybb_fit;
        amp_fit_bb = lsqnonlin(costfun,[1 1 1],[-Inf 0 0],[Inf Inf Inf],opts);
        ybb_model = interp1(tt3,amp_fit_bb(1)+cumtrapz(tt3,amp_fit_bb(2)*Pout3+amp_fit_bb(3)*Rout3),xbb_fit);
        [r_temp, p_temp] = corrcoef(ybb_fit,ybb_model);
        rsq2_bb(i) = r_temp(2)^2; pval2_bb(i) = p_temp(2);
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
        title(cb,'bbp650 (m-1)')
        subplot(2,2,4)
        patch([srs(i,2)-1 srs(i,1) srs(i,1) srs(i,2)-1],[0.08 0.08 0 0],[0.9 0.9 0.9],'edgecolor','none')
        hold on, patch([srs(i,2) srs(i,1)+1 srs(i,1)+1 srs(i,2)],[0.08 0.08 0 0],[0.9 0.9 0.9],'edgecolor','none'), hold off
        hold on, p1 = plot(xbb_fit,ybb_fit,'o'); hold off
        hold on, l1 = plot(tt3,[amp_fit_bb(1)+cumtrapz(tt3,amp_fit_bb(2)*Pout3+amp_fit_bb(3)*Rout3)],'k-'); hold off
        xlim([0-deltaday,1+deltaday]), ylim([min(ybb_fit)-0.0001 max(ybb_fit)+0.0001])
        legend([p1 l1],{'Datapoint','Fit'},'Fontsize',16)
        xlabel('Decimal day'), ylabel('Particle backscattering coeff. (m-1)')
        text(0.1,0.002,{['GCP = ' num2str(P2bb(i),2)];['R = ' num2str(R2bb(i),2)]})
        title(datestr(days(i)))        
    end
    %pause(0.1)
    clear ind_dd o_dd bb_dd sig_dd mld_dd mldsig_dd depth_dd date_dd Fs_dd
    clear o_fit bb_fit x_fit xo_fit yo_fit xbb_fit ybb_fit indsig_max depthsig_max
    clear ind_ofit ind_bbfit
    clear tt Pout Rout tt3 costfun 
    clear r_temp p_temp 
end

%% Plot mixed layer properties and rate estimates from oxygen and bbp

clf
subplot(5,1,1)
plot(dived.date,dived.sla,'k-'),
ylabel('SLA (cm)')
set(gca,'Fontsize',16,'box','off','Color','none')
set(gca,'xtick',days(1:4:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)+0.5])
subplot(5,1,2)
hold on
for i = 2:length(days)
    patch(days(i)+[srs(i-1,2)-1 srs(i-1,2)-1 srs(i,1) srs(i,1)],[0 300 300 0],[0.9 0.9 0.9],'edgecolor','none')
end
plot(dateo_ml,opt_ml,'color',[0 0 0])
hold off
set(gca,'xtick',days(1:4:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)+0.5])
ylim([204 215])%ylim([nanmin(opt_ml)-0.25 nanmax(opt_ml)+0.25])
set(gca,'Fontsize',16) 
ylabel('Median ML O2 (mmol m-3)')
subplot(5,1,3)
%{
bar(days(pval_o<0.05)+0.5,Po(pval_o<0.05),'k')
hold on, bar(days(pval_o<0.05)+0.5,-Ro(pval_o<0.05),'r'), hold off
hold on, createPatches(days(pval_o>=0.05)+0.5,Po(pval_o>=0.05),0.4,'k',0.15), hold off
hold on, createPatches(days(pval_o>=0.05)+0.5,-Ro(pval_o>=0.05),0.4,'r',0.15), hold off
ylim([-5.5 5.5])
title('1-day fit')
legend('GOP (p<0.05)','R (p<0.05)')
legend('boxoff')
ylabel('O2 rates (mmol O2 m-3 d-1)')
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)])
set(gca,'Fontsize',16)
%}
bar(days(pval2_o<0.05)+0.5,P2o(pval2_o<0.05),'k') % Significant GOP
hold on, bar(days(pval2_o<0.05)+0.5,-R2o(pval2_o<0.05),'r'), hold off % Significant R
hold on, createPatches(days(pval2_o>=0.05)+0.5,P2o(pval2_o>=0.05),0.4,'k',0.15), hold off
hold on, createPatches(days(pval2_o>=0.05)+0.5,-R2o(pval2_o>=0.05),0.4,'r',0.15), hold off
set(gca,'xtick',days(1:4:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)+0.5])
ylim([-5 5])
title('2-days fit')
legend('GOP (p<0.05)','R (p<0.05)')
legend('boxoff')
ylabel('O2 rates (mmol O2 m-3 d-1)')
xlabel('mm/dd 2015')
set(gca,'Fontsize',16)
subplot(5,1,4)
hold on
for i = 2:length(days)
    patch(days(i)+[srs(i-1,2)-1 srs(i-1,2)-1 srs(i,1) srs(i,1)],[0 1 1 0],[0.9 0.9 0.9],'edgecolor','none')
end
plot(datebb_ml,bb_ml,'color',[0 0 0])
hold off
set(gca,'xtick',days(1:4:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)+0.5])
ylim([nanmin(bb_ml)-0.0001 nanmax(bb_ml)+0.0001])
ylim([1e-4 6e-4])
set(gca,'Fontsize',16) 
ylabel('Median ML bbp650 (m-1)')
subplot(5,1,5)
%{
bar(days(pval_bb<0.05)+0.5,Pbb(pval_bb<0.05),'b') % Significant GOP
hold on, bar(days(pval_bb<0.05)+0.5,-Rbb(pval_bb<0.05),'g'), hold off % Significant R
hold on, createPatches(days(pval_bb>=0.05)+0.5,Pbb(pval_bb>=0.05),0.4,'b',0.15), hold off
hold on, createPatches(days(pval_bb>=0.05)+0.5,-Rbb(pval_bb>=0.05),0.4,'g',0.15), hold off
set(gca,'xtick',days(1:3:end))
datetick('x','mm/dd','keepticks')
xlim([datebb_ml(1) datebb_ml(end)+0.5])
ylim([-5e-4 5e-4])
title('1-day fit')
legend('GPP','R')
legend('boxoff')
ylabel('bbp rates (m-1 d-1)')
set(gca,'Fontsize',16)
%}
bar(days(pval2_bb<0.05)+0.5,P2bb(pval2_bb<0.05),'b') % Significant GOP
hold on, bar(days(pval2_bb<0.05)+0.5,-R2bb(pval2_bb<0.05),'g'), hold off % Significant R
hold on, createPatches(days(pval2_bb>=0.05)+0.5,P2bb(pval2_bb>=0.05),0.4,'b',0.15), hold off
hold on, createPatches(days(pval2_bb>=0.05)+0.5,-R2bb(pval2_bb>=0.05),0.4,'g',0.15), hold off
set(gca,'xtick',days(1:4:end))
datetick('x','mm/dd','keepticks')
xlim([dateo_ml(1) dateo_ml(end)+0.5])
ylim([-3.5e-4 3.5e-4])
title('2-days fit')
legend('GPP (p<0.05)','R (p<0.05)')
legend('boxoff')
ylabel('bbp rates (m-1 d-1)')
xlabel('mm/dd 2015')
set(gca,'Fontsize',16)

%% Whole cruise diel averages in the mixed layer
clear bins diel_opt diel_bbp650 md_opt md_bbp650 anom_opt anom_bbp650
%{
binwidth = 1/12; %1/24;
binrange = [binwidth/2 1-binwidth/2];

[diel_opt,bins] = dielbinning(opt_ml,dateo_ml-fix(dateo_ml),binwidth,binrange,'median');
[diel_bbp650,bins] = dielbinning(bb_ml,datebb_ml-fix(datebb_ml),binwidth,binrange,'median');

subplot(1,2,1)
plot([bins(end)-1 bins bins(1)+1],[diel_opt(end) diel_opt diel_opt(1)],'ko--')
xlim([0 1])
subplot(1,2,2)
plot([bins(end)-1 bins bins(1)+1],[diel_opt(end) diel_opt diel_opt(1)],'ko--')
%}

%% Plot spatial variability of GOP and R from optode data (2-days fit)

load ../hawaii.dat
coastline = hawaii;
coastline(:,1) = -coastline(:,1);
ind_n = find(isnan(coastline(:,1)));
clear hawaii
clf

subplot(2,1,1)
% Station ALOHA & coastline
    raloha = 6/59.79; 
    plot3(-158 + raloha*cos(0:0.01:2*pi),22.75 + raloha*sin(0:0.01:2*pi),0*(0:0.01:2*pi),'k','linewidth',1);
    text(-159,22.75,0,'Sta. ALOHA','color','k','Fontsize',16);
hold on
    for i = 1:(length(ind_n)-1)
        patch(-coastline((ind_n(i)+1):(ind_n(i+1)-1),1),coastline((ind_n(i)+1):(ind_n(i+1)-1),2),'k','edgecolor','none')
    end
    plot(dived.lon,dived.lat,'k--')
set(gca,'Fontsize',16)
axis equal
scatterbar3(daylon(pval2_o<0.05),daylat(pval2_o<0.05),P2o(pval2_o<0.05),0.05,'k',0.3)
hold off
xlim([-159 -156]),ylim([21 25.5])
daspect([min(daspect)*[1 1] 8])
set(gca,'Color','none','Box','off  ','CameraPosition',[-171.1661 5.1333 136.4140])
xlabel('Longitude (deg E)'), ylabel('Latitude (deg N)'),zlabel('GOP (mmol O2 m-3 d-1)')
subplot(2,1,2)
% Station ALOHA & coastline
    raloha = 6/59.79; 
    plot3(-158 + raloha*cos(0:0.01:2*pi),22.75 + raloha*sin(0:0.01:2*pi),0*(0:0.01:2*pi),'k','linewidth',1);
    text(-159,22.75,0,'Sta. ALOHA','color','k','Fontsize',16);
hold on
    for i = 1:(length(ind_n)-1)
        patch(-coastline((ind_n(i)+1):(ind_n(i+1)-1),1),coastline((ind_n(i)+1):(ind_n(i+1)-1),2),'k','edgecolor','none')
    end
    plot(dived.lon,dived.lat,'k--')
set(gca,'Fontsize',16)
axis equal
scatterbar3(daylon(pval2_o<0.05),daylat(pval2_o<0.05),R2o(pval2_o<0.05),0.05,'r',0.3)
hold off
xlim([-159 -156]),ylim([21 25.5])
daspect([min(daspect)*[1 1] 8])
set(gca,'Color','none','Box','off','CameraPosition',[-171.1661 5.1333 136.4140])
xlabel('Longitude (deg E)'), ylabel('Latitude (deg N)'),zlabel('R (mmol O2 m-3 d-1)')

%% Plot spatial variability of GPP and R from backscatter data (2-days fit)

load ../hawaii.dat
coastline = hawaii;
coastline(:,1) = -coastline(:,1);
ind_n = find(isnan(coastline(:,1)));
clear hawaii
clf

subplot(1,2,1)
% Station ALOHA & coastline
    raloha = 6/59.79; 
    plot3(-158 + raloha*cos(0:0.01:2*pi),22.75 + raloha*sin(0:0.01:2*pi),0*(0:0.01:2*pi),'k','linewidth',1);
    text(-159,22.75,0,'Sta. ALOHA','color','k','Fontsize',16);
hold on
    for i = 1:(length(ind_n)-1)
        patch(-coastline((ind_n(i)+1):(ind_n(i+1)-1),1),coastline((ind_n(i)+1):(ind_n(i+1)-1),2),'k','edgecolor','none')
    end
    plot(dived.lon,dived.lat,'k--')
set(gca,'Fontsize',16)
axis equal
scatterbar3(daylon(pval2_bb<0.05),daylat(pval2_bb<0.05),P2bb(pval2_bb<0.05),0.05,'b',0.3)
hold off
xlim([-159 -156]),ylim([21 25.5])
daspect([min(daspect)*[1 1] 4.5e-4])
set(gca,'Color','none','Box','off  ','CameraPosition',[-149.5743 4.8230 0.0072])
xlabel('Longitude (deg E)'), ylabel('Latitude (deg N)'),zlabel('GPP (m-1 d-1)')
subplot(1,2,2)
% Station ALOHA & coastline
    raloha = 6/59.79; 
    plot3(-158 + raloha*cos(0:0.01:2*pi),22.75 + raloha*sin(0:0.01:2*pi),0*(0:0.01:2*pi),'k','linewidth',1);
    text(-159,22.75,0,'Sta. ALOHA','color','k','Fontsize',16);
hold on
    for i = 1:(length(ind_n)-1)
        patch(-coastline((ind_n(i)+1):(ind_n(i+1)-1),1),coastline((ind_n(i)+1):(ind_n(i+1)-1),2),'k','edgecolor','none')
    end
    plot(dived.lon,dived.lat,'k--')
set(gca,'Fontsize',16)
axis equal
scatterbar3(daylon(pval2_bb<0.05),daylat(pval2_bb<0.05),R2bb(pval2_bb<0.05),0.05,'g',0.3)
hold off
xlim([-159 -156]),ylim([21 25.5])
daspect([min(daspect)*[1 1] 4.5e-4])
set(gca,'Color','none','Box','off','CameraPosition',[-149.5743 4.8230 0.0072])
xlabel('Longitude (deg E)'), ylabel('Latitude (deg N)'),zlabel('R (m-1 d-1)')