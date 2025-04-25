load('MatFile/Tempadlf');
Tempadnlf(:,:,:,1:4) = [];
Tempadnlf(:,:,:,end-3:end) = []; % 1944-2016
% Phase 1:  1944-1976; Phase 2: 1977-1999; Phase 3: 2000-2016
Temcom1 = mean(Tempadnlf(:,:,:,1:33),4); % Phase 1:  1944-1976
Temcom2 = mean(Tempadnlf(:,:,:,34:56),4); % Phase 2: 1977-1999
Temcom3 = mean(Tempadnlf(:,:,:,57:73),4); % Phase 3: 2000-2016
Temcom = cat(4,Temcom1,Temcom2,Temcom3);
%% longitude mean (zonal mean)
% Indian Ocean
phases = {'1944-1976','1977-1999','2000-2016'};
lonmI1 = permute(nanmean(Temcom(40:100,83:116,:,:),1),[2 3 4 1]); % 40E-100E, 8S-25N
lonmI2 = permute(nanmean(Temcom(20:146,20:82,:,:),1),[2 3 4 1]);
lonmI = cat(1,lonmI2,lonmI1);
lonmI1m = permute(nanmean(Temp(40:100,83:116,:,:),1),[2 3 4 1]); % climotology
lonmI2m = permute(nanmean(Temp(20:146,20:82,:,:),1),[2 3 4 1]);
lonmIm = cat(1,lonmI2m,lonmI1m);
%
close all;
for k = 1:3
    Fig = figure('position',[100 100 800 400]);
    contourf(latData(20:116,1),-depthData,lonmI(:,:,k)',[-2:0.02:2],'linestyle','none')
    hold on
    caxis([-0.1,0.1])
%     [ac ah] = contour(latData(20:116,1),-depthData,(nanmean(lonmIm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    ch = colorbar;
    set(ch,'Ticks',[-1:0.02:0.1]);
    set(gca,'XLim',[-70,20]);
    set(gca,'XTick',[-60:20:20]);
    set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN'},'FontSize',14);
    set(gca,'YLim',[-2000,0]);
    set(gca,'YTick',[-2000:200:0],'FontSize',14);
    set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
    ylabel('Depth (m)')
    title([phases{k},' Com.zonal mean temperature in Indian Ocean'])
    print(Fig,['E:\figures\IAP\Yearly\8yrs_nonlinearDetrend\IPO.phase',num2str(k),'.Indian.zonalmean.tem.png'],'-dpng','-r600')
end
%% Pacific 
lonmP1 = permute(nanmean(Temcom(146:293,20:82,:,:),1),[2 3 4 1]); % 146E-67W, 70S-8S
lonmP2 = permute(nanmean(Temcom(100:293,83:101,:,:),1),[2 3 4 1]); % 100E-67W, 8S-10N
lonmP3 = permute(nanmean(Temcom(100:270,102:109,:,:),1),[2 3 4 1]); % 100E-90W, 10N-18N
lonmP4 = permute(nanmean(Temcom(100:255,110:156,:,:),1),[2 3 4 1]); % 100E-105W, 18N-65N
lonmP = cat(1,lonmP1,lonmP2,lonmP3,lonmP4);
max(lonmP,[],'all')
min(lonmP,[],'all')
lonmP1m = permute(nanmean(Temp(146:293,20:82,:,:),1),[2 3 4 1]); % climotology
lonmP2m = permute(nanmean(Temp(100:293,83:101,:,:),1),[2 3 4 1]); 
lonmP3m = permute(nanmean(Temp(100:270,102:109,:,:),1),[2 3 4 1]);
lonmP4m = permute(nanmean(Temp(100:255,110:156,:,:),1),[2 3 4 1]); 
lonmPm = cat(1,lonmP1m,lonmP2m,lonmP3m,lonmP4m);
%
close all;
for k = 1:3
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:156,1),-depthData,lonmP(:,:,k)',[-0.5:0.01:0.5],'linestyle','none')
caxis([-0.1,0.1])
hold on
[ac ah] = contour(latData(20:156,1),-depthData,(nanmean(lonmPm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.1:0.02:0.1]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:200:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
ylabel('Depth (m)')
title([phases{k},' Com.zonal mean temperature in Pacific'])
print(Fig,['E:\figures\IAP\Yearly\8yrs_nonlinearDetrend\IPO.phase',num2str(k),'.Pacific.zonalmean.tem.png'],'-dpng','-r600')
end
%% Atlantic
TemAtl = cat(1,Temcom(255:360,:,:,:),Temcom(1:20,:,:,:));
lonAtl = [lonData(255:360)-360;lonData(1:20)];
lonmA1 = permute(nanmean(TemAtl(39:end,20:101,:,:),1),[2 3 4 1]); % 67W-20E, -70S-10N
lonmA2 = permute(nanmean(TemAtl(16:end,102:109,:,:),1),[2 3 4 1]); % 90W-20E, 11N-18N 
lonmA3 = permute(nanmean(TemAtl(:,110:161,:,:),1),[2 3 4 1]); % 105W-20E, 19N-70N
lonmA = cat(1,lonmA1,lonmA2,lonmA3);
max(lonmA,[],'all')
min(lonmA,[],'all')
TemAtlm = cat(1,Temp(255:360,:,:,:),Temp(1:20,:,:,:)); % climotology
lonmA1m = permute(nanmean(TemAtlm(39:end,20:101,:,:),1),[2 3 4 1]); 
lonmA2m = permute(nanmean(TemAtlm(16:end,102:109,:,:),1),[2 3 4 1]); 
lonmA3m = permute(nanmean(TemAtlm(:,110:161,:,:),1),[2 3 4 1]); 
lonmAm = cat(1,lonmA1m,lonmA2m,lonmA3m);
%
close all;
for k = 1:3
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:161,1),-depthData,lonmA(:,:,k)',[-1:0.01:1],'linestyle','none')
caxis([-0.1,0.1])
hold on
[ac ah] = contour(latData(20:161,1),-depthData,(nanmean(lonmAm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.1:0.02:0.1]);
set(gca,'XLim',[-70,70]);
set(gca,'XTick',[-60:20:70]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:200:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
ylabel('Depth (m)')
title([phases{k},' Com.zonal mean temperature in Atlantic'])
print(Fig,['E:\figures\IAP\Yearly\8yrs_nonlinearDetrend\IPO.phase',num2str(k),'.Atlantic.zonalmean.tem.png'],'-dpng','-r600')
end
%% meridional mean 
for k = 1:3
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[50:111],'40S-20N',k,phases,TPI_filtered,lonData,depthData);
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[20:40],'70S-50S',k,phases,TPI_filtered,lonData,depthData);
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[40:60],'50S-30S',k,phases,TPI_filtered,lonData,depthData);
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[60:116],'30S-25N',k,phases,TPI_filtered,lonData,depthData);
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[116:141],'25N-50N',k,phases,TPI_filtered,lonData,depthData);
end
%% different depth
for k = 1:3
    diffdepth(Temcom(:,:,:,k),1,k,phases,TPI_filtered,lonData,latData,'0000',1.2,0.4)
    diffdepth(Temcom(:,:,:,k),7,k,phases,TPI_filtered,lonData,latData,'0050',1.2,0.4)
    diffdepth(Temcom(:,:,:,k),12,k,phases,TPI_filtered,lonData,latData,'0100',1.2,0.4)
    diffdepth(Temcom(:,:,:,k),17,k,phases,TPI_filtered,lonData,latData,'0200',1.2,0.3)
    diffdepth(Temcom(:,:,:,k),23,k,phases,TPI_filtered,lonData,latData,'0500',1.2,0.2)
    diffdepth(Temcom(:,:,:,k),29,k,phases,TPI_filtered,lonData,latData,'0800',1.2,0.1)
    diffdepth(Temcom(:,:,:,k),32,k,phases,TPI_filtered,lonData,latData,'1000',1.2,0.1)
    diffdepth(Temcom(:,:,:,k),37,k,phases,TPI_filtered,lonData,latData,'1500',1.2,0.05)
    diffdepth(Temcom(:,:,:,k),41,k,phases,TPI_filtered,lonData,latData,'2000',1.2,0.05)
end



function [ts] = lonmean(var,lats,latData)
% average along latitude (all longitudes)
% var is lon*lat*depth*time
% ts is lat*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,1)./nansum(cos(latData(lats)'/180*pi).*var2,1),[2 3 4 1]);
end
function [ts] = latmean(var,lats,latData)
% average along longitude (all latitudes)
% var is lon*lat*depth*time
% ts is lon*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,2)./nansum(cos(latData(lats)'/180*pi).*var2,2),[1 3 4 2]);
end
function contourVARra(var,val,ctr,lon1,lon2,laTempa,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[laTempa lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
end
function meridionalmean(Tempadlf,Temp,latData,lats,names,k,phases,index,lonData,depthData)
    latm = latmean(Tempadlf,lats,latData);
    size(latm)
    latmc = latmean(Temp,lats,latData); % climotology
    max(latm,[],'all')
    min(latm,[],'all')
    %
    close all;
    Fig = figure('position',[100 100 800 400]);
    contourf(lonData,-depthData,latm',[-0.5:0.01:0.5],'linestyle','none')
    caxis([-0.1,0.1])
    hold on
    [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    ch = colorbar;
    set(ch,'Ticks',[-1:0.02:1]);
    set(gca,'XLim',[0,360]);
    set(gca,'XTick',[0:60:360]);
    set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
    set(gca,'YLim',[-2000,0]);
    set(gca,'YTick',[-2000:200:0],'FontSize',14);
    set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
    ylabel('Depth (m)')
    title([phases{k},' Com.meridional mean temperature (',names,')'])
    print(Fig,['E:\figures\IAP\Yearly\8yrs_nonlinearDetrend\IPO.phase',num2str(k),'.meridionalmean',names,'.tem.png'],'-dpng','-r600')
end
function diffdepth(Temcom,nlevel,k,phases,index,lonData,latData,names,val,ticks)
    close all;
    Fig = figure('position',[100 100 600 300]);
    contourVARra(Temcom(:,:,nlevel),[-val:ticks/10:val],1,0,360,-90,90,lonData,latData);
    caxis([-ticks,ticks]);
    ch = colorbar;
    set(ch,'Ticks',[-ticks:ticks/5:ticks]);
    % testdots(h05,'k',lonData,latData)
    title([phases{k},' Com.T',' ',names,'m']);
    print(Fig,['E:\figures\IAP\Yearly\8yrs_nonlinearDetrend\IPO.phase',num2str(k),'.diffdepth.',names,'m.tem.png'],'-dpng','-r600')
end





    