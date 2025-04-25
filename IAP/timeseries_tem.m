addpath E:\1_matlab\help\EEMD_code\;
amI01 = IndianA(Tempa,latData,1); % Indian Ocean surface 0000m
rslt=eemd(amI01,0.2,2000);
IndianAplot(amI01,rslt,'0000m',0.6,0.2,0.3,0.1,'Indian Ocean','IndianOcean')
%%
amI02 = IndianA(Tempa,latData,7); % Indian Ocean surface 0050m
rslt=eemd(amI02,0.2,2000);
IndianAplot(amI02,rslt,'0050m',0.6,0.2,0.3,0.1,'Indian Ocean','IndianOcean')
amI03 = IndianA(Tempa,latData,17); % Indian Ocean surface 0200m
rslt=eemd(amI03,0.2,2000);
IndianAplot(amI03,rslt,'0200m',0.4,0.1,0.2,0.1,'Indian Ocean','IndianOcean')
amI04 = IndianA(Tempa,latData,23); % Indian Ocean surface 0500m
rslt=eemd(amI04,0.2,2000);
IndianAplot(amI04,rslt,'0500m',0.2,0.1,0.1,0.05,'Indian Ocean','IndianOcean')
amI05 = IndianA(Tempa,latData,32); % Indian Ocean surface 1000m
rslt=eemd(amI05,0.2,2000);
IndianAplot(amI05,rslt,'1000m',0.2,0.1,0.1,0.05,'Indian Ocean','IndianOcean')
amI06 = IndianA(Tempa,latData,41); % Indian Ocean surface 2000m
rslt=eemd(amI06,0.2,2000);
IndianAplot(amI06,rslt,'2000m',0.03,0.01,0.03,0.01,'Indian Ocean','IndianOcean')
%% SEIO (100-116E,35-15S) Li et al. (2019)
Tempaz = permute(Tempa,[1 2 4 3]); % transfer
[amI11z amI11] = areamean(Tempaz(:,:,:,1),[100:116],[55:75],latData); 
rslt=eemd(amI11,0.2,2000);
IndianAplot(amI11,rslt,'0000m',1,0.2,0.8,0.2,'Southeast IO','IndianPart1.')
[amI12z amI12] = areamean(Tempaz(:,:,:,7),[100:116],[55:75],latData); 
rslt=eemd(amI12,0.2,2000);
IndianAplot(amI12,rslt,'0050m',1,0.2,0.8,0.2,'Southeast IO','IndianPart1.')
[amI13z amI13] = areamean(Tempaz(:,:,:,17),[100:116],[55:75],latData); 
rslt=eemd(amI13,0.2,2000);
IndianAplot(amI13,rslt,'0200m',0.8,0.2,0.6,0.2,'Southeast IO','IndianPart1.')
[amI14z amI14] = areamean(Tempaz(:,:,:,23),[100:116],[55:75],latData); 
rslt=eemd(amI14,0.2,2000);
IndianAplot(amI14,rslt,'0500m',0.2,0.05,0.2,0.05,'Southeast IO','IndianPart1.')
[amI15z amI15] = areamean(Tempaz(:,:,:,32),[100:116],[55:75],latData); 
rslt=eemd(amI15,0.2,2000);
IndianAplot(amI15,rslt,'1000m',0.4,0.1,0.4,0.1,'Southeast IO','IndianPart1.')
[amI16z amI16] = areamean(Tempaz(:,:,:,41),[100:116],[55:75],latData); 
rslt=eemd(amI16,0.2,2000);
IndianAplot(amI16,rslt,'2000m',0.03,0.01,0.03,0.01,'Southeast IO','IndianPart1.')
%% SWIO (30-80E,65-35S) 
[amI21z amI21] = areamean(Tempaz(:,:,:,1),[30:80],[25:55],latData); 
rslt=eemd(amI21,0.2,2000);
IndianAplot(amI21,rslt,'0000m',0.6,0.2,0.6,0.2,'Southwest IO','IndianPart2.')
[amI22z amI22] = areamean(Tempaz(:,:,:,7),[30:80],[25:55],latData); 
rslt=eemd(amI22,0.2,2000);
IndianAplot(amI22,rslt,'0050m',0.8,0.2,0.6,0.2,'Southwest IO','IndianPart2.')
[amI23z amI23] = areamean(Tempaz(:,:,:,17),[30:80],[25:55],latData); 
rslt=eemd(amI23,0.2,2000);
IndianAplot(amI23,rslt,'0200m',0.6,0.2,0.4,0.2,'Southwest IO','IndianPart2.')
[amI24z amI24] = areamean(Tempaz(:,:,:,23),[30:80],[25:55],latData); 
rslt=eemd(amI24,0.2,2000);
IndianAplot(amI24,rslt,'0500m',0.4,0.2,0.2,0.1,'Southwest IO','IndianPart2.')
[amI25z amI25] = areamean(Tempaz(:,:,:,32),[30:80],[25:55],latData); 
rslt=eemd(amI25,0.2,2000);
IndianAplot(amI25,rslt,'1000m',0.2,0.1,0.2,0.1,'Southwest IO','IndianPart2.')
[amI26z amI26] = areamean(Tempaz(:,:,:,41),[30:80],[25:55],latData); 
rslt=eemd(amI26,0.2,2000);
IndianAplot(amI26,rslt,'2000m',0.06,0.02,0.06,0.02,'Southwest IO','IndianPart2.')
%% NIO (40-100E,0-25N) 
[amI31z amI31] = areamean(Tempaz(:,:,:,1),[40:100],[91:116],latData); 
rslt=eemd(amI31,0.2,2000);
IndianAplot(amI31,rslt,'0000m',1,0.2,0.8,0.2,'North IO','IndianPart3.')
[amI32z amI32] = areamean(Tempaz(:,:,:,7),[40:100],[91:116],latData); 
rslt=eemd(amI32,0.2,2000);
IndianAplot(amI32,rslt,'0050m',1,0.2,0.8,0.2,'North IO','IndianPart3.')
[amI33z amI33] = areamean(Tempaz(:,:,:,17),[40:100],[91:116],latData); 
rslt=eemd(amI33,0.2,2000);
IndianAplot(amI33,rslt,'0200m',0.6,0.2,0.4,0.2,'North IO','IndianPart3.')
[amI34z amI34] = areamean(Tempaz(:,:,:,23),[40:100],[91:116],latData); 
rslt=eemd(amI34,0.2,2000);
IndianAplot(amI34,rslt,'0500m',0.2,0.05,0.2,0.05,'North IO','IndianPart3.')
[amI35z amI35] = areamean(Tempaz(:,:,:,32),[40:100],[91:116],latData); 
rslt=eemd(amI35,0.2,2000);
IndianAplot(amI35,rslt,'1000m',0.2,0.1,0.2,0.1,'North IO','IndianPart3.')
[amI36z amI36] = areamean(Tempaz(:,:,:,41),[40:100],[91:116],latData); 
rslt=eemd(amI36,0.2,2000);
IndianAplot(amI36,rslt,'2000m',0.06,0.02,0.06,0.02,'North IO','IndianPart3.')
%% Pacific
amP01 = PacificA(Tempa,latData,1); % Pacific surface 0000m
rslt=eemd(amP01,0.2,2000);
IndianAplot(amP01,rslt,'0000m',0.6,0.2,0.4,0.1,'Pacific','Pacific')
amP02 = PacificA(Tempa,latData,7); % Pacific surface 0050m
rslt=eemd(amP02,0.2,2000);
IndianAplot(amP02,rslt,'0050m',0.6,0.2,0.4,0.1,'Pacific','Pacific')
amP03 = PacificA(Tempa,latData,17); % Pacific surface 0200m
rslt=eemd(amP03,0.2,2000);
IndianAplot(amP03,rslt,'0200m',0.4,0.1,0.3,0.1,'Pacific','Pacific')
amP04 = PacificA(Tempa,latData,23); % Pacific surface 0500m
rslt=eemd(amP04,0.2,2000);
IndianAplot(amP04,rslt,'0500m',0.1,0.05,0.1,0.05,'Pacific','Pacific')
amP05 = PacificA(Tempa,latData,32); % Pacific surface 1000m
rslt=eemd(amP05,0.2,2000);
IndianAplot(amP05,rslt,'1000m',0.08,0.02,0.05,0.025,'Pacific','Pacific')
amP06 = PacificA(Tempa,latData,41); % Pacific surface 2000m
rslt=eemd(amP06,0.2,2000);
IndianAplot(amP06,rslt,'2000m',0.02,0.01,0.02,0.01,'Pacific','Pacific')
%% NP (150E-150W,25-50N) 
[amP11z amP11] = areamean(Tempaz(:,:,:,1),[150:210],[116:141],latData); 
rslt=eemd(amP11,0.2,2000);
IndianAplot(amP11,rslt,'0000m',1.2,0.2,1,0.2,'North Pacific','PacificPart1.')
[amP12z amP12] = areamean(Tempaz(:,:,:,7),[150:210],[116:141],latData); 
rslt=eemd(amP12,0.2,2000);
IndianAplot(amP12,rslt,'0050m',1,0.2,0.8,0.2,'North Pacific','PacificPart1.')
[amP13z amP13] = areamean(Tempaz(:,:,:,17),[150:210],[116:141],latData); 
rslt=eemd(amP13,0.2,2000);
IndianAplot(amP13,rslt,'0200m',0.6,0.2,0.6,0.2,'North Pacific','PacificPart1.')
[amP14z amP14] = areamean(Tempaz(:,:,:,23),[150:210],[116:141],latData); 
rslt=eemd(amP14,0.2,2000);
IndianAplot(amP14,rslt,'0500m',0.3,0.1,0.3,0.1,'North Pacific','PacificPart1.')
[amP15z amP15] = areamean(Tempaz(:,:,:,32),[150:210],[116:141],latData); 
rslt=eemd(amP15,0.2,2000);
IndianAplot(amP15,rslt,'1000m',0.1,0.05,0.1,0.05,'North Pacific','PacificPart1.')
[amP16z amP16] = areamean(Tempaz(:,:,:,41),[150:210],[116:141],latData); 
rslt=eemd(amP16,0.2,2000);
IndianAplot(amP16,rslt,'2000m',0.02,0.01,0.02,0.01,'North Pacific','PacificPart1.')
%% TP (160W-90W,10S-10N) 
[amP21z amP21] = areamean(Tempaz(:,:,:,1),[200:270],[80:101],latData); 
rslt=eemd(amP21,0.2,2000);
IndianAplot(amP21,rslt,'0000m',1.5,0.5,1.2,0.2,'Tropical Pacific','PacificPart2.')
[amP22z amP22] = areamean(Tempaz(:,:,:,7),[200:270],[80:101],latData); 
rslt=eemd(amP22,0.2,2000);
IndianAplot(amP22,rslt,'0050m',1.5,0.5,1.2,0.2,'Tropical Pacific','PacificPart2.')
[amP23z amP23] = areamean(Tempaz(:,:,:,17),[200:270],[80:101],latData); 
rslt=eemd(amP23,0.2,2000);
IndianAplot(amP23,rslt,'0200m',0.6,0.2,0.4,0.2,'Tropical Pacific','PacificPart2.')
[amP24z amP24] = areamean(Tempaz(:,:,:,23),[200:270],[80:101],latData); 
rslt=eemd(amP24,0.2,2000);
IndianAplot(amP24,rslt,'0500m',0.2,0.05,0.2,0.05,'Tropical Pacific','PacificPart2.')
[amP25z amP25] = areamean(Tempaz(:,:,:,32),[200:270],[80:101],latData); 
rslt=eemd(amP25,0.2,2000);
IndianAplot(amP25,rslt,'1000m',0.1,0.05,0.1,0.05,'Tropical Pacific','PacificPart2.')
[amP26z amP26] = areamean(Tempaz(:,:,:,41),[200:270],[80:101],latData); 
rslt=eemd(amP26,0.2,2000);
IndianAplot(amP26,rslt,'2000m',0.06,0.02,0.06,0.02,'Tropical Pacific','PacificPart2.')
%% SP (170E-120W,45S-30S) 
[amP31z amP31] = areamean(Tempaz(:,:,:,1),[170:240],[45:60],latData); 
rslt=eemd(amP31,0.2,2000);
IndianAplot(amP31,rslt,'0000m',1,0.2,0.8,0.2,'South Pacific','PacificPart3.')
[amP32z amP32] = areamean(Tempaz(:,:,:,7),[170:240],[45:60],latData); 
rslt=eemd(amP32,0.2,2000);
IndianAplot(amP32,rslt,'0050m',1,0.2,0.8,0.2,'South Pacific','PacificPart3.')
[amP33z amP33] = areamean(Tempaz(:,:,:,17),[170:240],[45:60],latData); 
rslt=eemd(amP33,0.2,2000);
IndianAplot(amP33,rslt,'0200m',0.6,0.2,0.6,0.2,'South Pacific','PacificPart3.')
[amP34z amP34] = areamean(Tempaz(:,:,:,23),[170:240],[45:60],latData); 
rslt=eemd(amP34,0.2,2000);
IndianAplot(amP34,rslt,'0500m',0.3,0.1,0.3,0.1,'South Pacific','PacificPart3.')
[amP35z amP35] = areamean(Tempaz(:,:,:,32),[170:240],[45:60],latData); 
rslt=eemd(amP35,0.2,2000);
IndianAplot(amP35,rslt,'1000m',0.2,0.1,0.2,0.1,'South Pacific','PacificPart3.')
[amP36z amP36] = areamean(Tempaz(:,:,:,41),[170:240],[45:60],latData); 
rslt=eemd(amP36,0.2,2000);
IndianAplot(amP36,rslt,'2000m',0.04,0.02,0.04,0.02,'South Pacific','PacificPart3.')
%% Atlantic
TemAtl = cat(1,Tempa(255:360,:,:,:),Tempa(1:20,:,:,:));
lonAtl = [lonData(255:360)-360;lonData(1:20)];
%% Atlantic
amA01 = AtlanticA(TemAtl,latData,1); % Atlantic surface 0000m
rslt=eemd(amA01,0.2,2000);
%%
IndianAplot(amA01,rslt,'0000m',0.6,0.2,0.4,0.1,'Atlantic','Atlantic')
amA02 = AtlanticA(TemAtl,latData,7); % Atlantic surface 0050m
rslt=eemd(amA02,0.2,2000);
IndianAplot(amA02,rslt,'0050m',0.6,0.2,0.4,0.1,'Atlantic','Atlantic')
amA03 = AtlanticA(TemAtl,latData,17); % Atlantic surface 0200m
rslt=eemd(amA03,0.2,2000);
IndianAplot(amA03,rslt,'0200m',0.4,0.1,0.2,0.1,'Atlantic','Atlantic')
amA04 = AtlanticA(TemAtl,latData,23); % Atlantic surface 0500m
rslt=eemd(amA04,0.2,2000);
IndianAplot(amA04,rslt,'0500m',0.4,0.1,0.1,0.05,'Atlantic','Atlantic')
amA05 = AtlanticA(TemAtl,latData,32); % Atlantic surface 1000m
rslt=eemd(amA05,0.2,2000);
IndianAplot(amA05,rslt,'1000m',0.16,0.04,0.05,0.025,'Atlantic','Atlantic')
amA06 = AtlanticA(TemAtl,latData,41); % Atlantic surface 2000m
rslt=eemd(amA06,0.2,2000);
IndianAplot(amA06,rslt,'2000m',0.1,0.05,0.06,0.02,'Atlantic','Atlantic')
%% NA (60W-0,55-65N) 
[amA11z amA11] = areamean(Tempaz(:,:,:,1),[300:360],[146:156],latData); 
rslt=eemd(amA11,0.2,2000);
IndianAplot(amA11,rslt,'0000m',1.2,0.2,1,0.2,'North Atlantic','AtlanticPart1.')
[amA12z amA12] = areamean(Tempaz(:,:,:,7),[300:360],[146:156],latData); 
rslt=eemd(amA12,0.2,2000);
IndianAplot(amA12,rslt,'0050m',1,0.2,1,0.2,'North Atlantic','AtlanticPart1.')
[amA13z amA13] = areamean(Tempaz(:,:,:,17),[300:360],[146:156],latData); 
rslt=eemd(amA13,0.2,2000);
IndianAplot(amA13,rslt,'0200m',1,0.2,1,0.2,'North Atlantic','AtlanticPart1.')
[amA14z amA14] = areamean(Tempaz(:,:,:,23),[300:360],[146:156],latData); 
rslt=eemd(amA14,0.2,2000);
IndianAplot(amA14,rslt,'0500m',0.6,0.2,0.6,0.2,'North Atlantic','AtlanticPart1.')
[amA15z amA15] = areamean(Tempaz(:,:,:,32),[300:360],[146:156],latData); 
rslt=eemd(amA15,0.2,2000);
IndianAplot(amA15,rslt,'1000m',0.6,0.2,0.6,0.2,'North Atlantic','AtlanticPart1.')
[amA16z amA16] = areamean(Tempaz(:,:,:,41),[300:360],[146:156],latData); 
rslt=eemd(amA16,0.2,2000);
IndianAplot(amA16,rslt,'2000m',0.4,0.2,0.4,0.2,'North Atlantic','AtlanticPart1.')
%% TA (67W-0,10S-20N) 
[amA21z amA21] = areamean(Tempaz(:,:,:,1),[293:360],[80:101],latData); 
rslt=eemd(amA21,0.2,2000);
IndianAplot(amA21,rslt,'0000m',1.2,0.2,1,0.2,'Tropical Atlantic','AtlanticPart2.')
[amA22z amA22] = areamean(Tempaz(:,:,:,7),[293:360],[80:101],latData); 
rslt=eemd(amA22,0.2,2000);
IndianAplot(amA22,rslt,'0050m',1,0.2,0.8,0.2,'Tropical Atlantic','AtlanticPart2.')
[amA23z amA23] = areamean(Tempaz(:,:,:,17),[293:360],[80:101],latData); 
rslt=eemd(amA23,0.2,2000);
IndianAplot(amA23,rslt,'0200m',0.4,0.2,0.4,0.2,'Tropical Atlantic','AtlanticPart2.')
[amA24z amA24] = areamean(Tempaz(:,:,:,23),[293:360],[80:101],latData); 
rslt=eemd(amA24,0.2,2000);
IndianAplot(amA24,rslt,'0500m',0.3,0.1,0.2,0.1,'Tropical Atlantic','AtlanticPart2.')
[amA25z amA25] = areamean(Tempaz(:,:,:,32),[293:360],[80:101],latData); 
rslt=eemd(amA25,0.2,2000);
IndianAplot(amA25,rslt,'1000m',0.2,0.1,0.1,0.05,'Tropical Atlantic','AtlanticPart2.')
[amA26z amA26] = areamean(Tempaz(:,:,:,41),[293:360],[80:101],latData); 
rslt=eemd(amA26,0.2,2000);
IndianAplot(amA26,rslt,'2000m',0.1,0.05,0.08,0.02,'Tropical Atlantic','AtlanticPart2.')
%% SA (67W-0,30S-50S) 
[amA31z amA31] = areamean(Tempaz(:,:,:,1),[293:360],[40:60],latData); 
rslt=eemd(amA31,0.2,2000);
IndianAplot(amA31,rslt,'0000m',1,0.2,1,0.2,'South Atlantic','AtlanticPart3.')
[amA32z amA32] = areamean(Tempaz(:,:,:,7),[293:360],[40:60],latData); 
rslt=eemd(amA32,0.2,2000);
IndianAplot(amA32,rslt,'0050m',0.8,0.2,0.8,0.2,'South Atlantic','AtlanticPart3.')
[amA33z amA33] = areamean(Tempaz(:,:,:,17),[293:360],[40:60],latData); 
rslt=eemd(amA33,0.2,2000);
IndianAplot(amA33,rslt,'0200m',0.6,0.2,0.6,0.2,'South Atlantic','AtlanticPart3.')
[amA34z amA34] = areamean(Tempaz(:,:,:,23),[293:360],[40:60],latData); 
rslt=eemd(amA34,0.2,2000);
IndianAplot(amA34,rslt,'0500m',0.5,0.1,0.3,0.1,'South Atlantic','AtlanticPart3.')
[amA35z amA35] = areamean(Tempaz(:,:,:,32),[293:360],[40:60],latData); 
rslt=eemd(amA35,0.2,2000);
IndianAplot(amA35,rslt,'1000m',0.3,0.1,0.2,0.1,'South Atlantic','AtlanticPart3.')
[amA36z amA36] = areamean(Tempaz(:,:,:,41),[293:360],[40:60],latData); 
rslt=eemd(amA36,0.2,2000);
IndianAplot(amA36,rslt,'2000m',0.12,0.04,0.12,0.04,'South Atlantic','AtlanticPart3.')







function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end
function [amI01] = IndianA(Tempa,latData,nlevel)
    weit1 = Tempa(40:100,83:116,nlevel,1); 
    weit1(find(isnan(weit1) == 0)) = 1; % weight
    weit2 = Tempa(20:146,20:82,nlevel,1);
    weit2(find(isnan(weit2) == 0)) = 1; % weight
    amI01 = permute((nansum((cos(latData(83:116)'/180*pi)).*Tempa(40:100,83:116,nlevel,:),[1 2])+nansum((cos(latData(20:82)'/180*pi)).*Tempa(20:146,20:82,nlevel,:),[1 2]))...
    ./(nansum(cos(latData(83:116)'/180*pi).*weit1,'all')+nansum(cos(latData(20:82)'/180*pi).*weit2,'all')),[4 1 2 3]);
end
function IndianAplot(amI01,rslt,lvlstr,ylim1,yint1,ylim2,yint2,name1,namep)
    close all;
    ftsz = 14;
    Fig = figure('Position',[100 100 800 800])
    subplot(2,1,1)
    plot(amI01,'k','LineWidth',1.5)
    hold on
    plot(rslt(:,7),'r','LineWidth',1.5)
    plot(amI01-detrend(amI01),'b','LineWidth',1.5)
    plot(zeros(1,81))
    set(gca,'XLim',[1,82]);
    set(gca,'XTick',[1:10:81]);
    set(gca,'XTickLabel',[1940:10:2020],'FontSize',ftsz);
    set(gca,'YLim',[-ylim1,ylim1]);
    set(gca,'YTick',[-ylim1:yint1:ylim1]);
    xlabel('YEAR','FontSize',ftsz),ylabel('Temperature','FontSize',ftsz);
    legend([name1,' ', lvlstr],'Non-linear trend','Linear trend','Location','north','Orientation','horizontal')
    legend('boxoff')
    subplot(2,1,2)
    plot(amI01-rslt(:,7),'r','LineWidth',1.5)
    hold on
    plot(detrend(amI01),'b','LineWidth',1.5)
    plot(zeros(1,81))
    set(gca,'XLim',[1,82]);
    set(gca,'XTick',[1:10:81]);
    set(gca,'XTickLabel',[1940:10:2020],'FontSize',ftsz);
    set(gca,'YLim',[-ylim2,ylim2]);
    set(gca,'YTick',[-ylim2:yint2:ylim2]);
    xlabel('YEAR','FontSize',ftsz),ylabel('Temperature','FontSize',ftsz);
    legend('Non-linear','Linear','Location','north','Orientation','horizontal')
    legend('boxoff')
    title('Detrended'); 
    print(Fig,['E:\figures\IAP\Yearly\20220713TimeSeries\',namep,lvlstr,'.png'],'-dpng','-r600')

end
function [amP01] = PacificA(Tempa,latData,nlevel)
    weit1 = Tempa(146:293,20:82,nlevel,1); weit1(find(isnan(weit1) == 0)) = 1; % weight
    weit2 = Tempa(100:293,83:101,nlevel,1); weit2(find(isnan(weit2) == 0)) = 1; % weight
    weit3 = Tempa(100:270,102:109,nlevel,1); weit3(find(isnan(weit3) == 0)) = 1; % weight
    weit4 = Tempa(100:255,110:156,nlevel,1); weit4(find(isnan(weit4) == 0)) = 1; % weight
    amP01 = permute((nansum((cos(latData(20:82)'/180*pi)).*Tempa(146:293,20:82,nlevel,:),[1 2])...
        +nansum((cos(latData(83:101)'/180*pi)).*Tempa(100:293,83:101,nlevel,:),[1 2])...
        +nansum((cos(latData(102:109)'/180*pi)).*Tempa(100:270,102:109,nlevel,:),[1 2])...
        +nansum((cos(latData(110:156)'/180*pi)).*Tempa(100:255,110:156,nlevel,:),[1 2]))...
        ./(nansum(cos(latData(20:82)'/180*pi).*weit1,'all')...
        +nansum(cos(latData(83:101)'/180*pi).*weit2,'all')...
        +nansum(cos(latData(102:109)'/180*pi).*weit3,'all')...
        +nansum(cos(latData(110:156)'/180*pi).*weit4,'all')),[4 1 2 3]);
end
function [amA01] = AtlanticA(Tempa,latData,nlevel)
    weit1 = Tempa(39:end,20:101,nlevel,1); weit1(find(isnan(weit1) == 0)) = 1; % weight
    weit2 = Tempa(16:end,102:109,nlevel,1); weit2(find(isnan(weit2) == 0)) = 1; % weight
    weit3 = Tempa(:,110:161,nlevel,1); weit3(find(isnan(weit3) == 0)) = 1; % weight
    amA01 = permute((nansum((cos(latData(20:101)'/180*pi)).*Tempa(39:end,20:101,nlevel,:),[1 2])...
        +nansum((cos(latData(102:109)'/180*pi)).*Tempa(16:end,102:109,nlevel,:),[1 2])...
        +nansum((cos(latData(110:161)'/180*pi)).*Tempa(:,110:161,nlevel,:),[1 2]))...
        ./(nansum(cos(latData(20:101)'/180*pi).*weit1,'all')...
        +nansum(cos(latData(102:109)'/180*pi).*weit2,'all')...
        +nansum(cos(latData(110:161)'/180*pi).*weit3,'all')),[4 1 2 3]);
end

