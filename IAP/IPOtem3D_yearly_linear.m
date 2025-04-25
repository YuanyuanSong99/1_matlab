addpath E:\1_matlab\help;
load("MatFile\Tempadlf.mat"); % detrended and 8-yr filtered temperature
Tempadlf(:,:,:,1:4) = []; Tempadlf(:,:,:,end-3:end) = []; 
load("MatFile\TPI_filtered.mat");
TPI_filtered(1:4) = []; TPI_filtered(end-3:end) = [];
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
%%
close all;
Fig = figure('Position',[100 100 1000 500]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[20 380]);
hold on
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
% Indian Ocean
m_line([20 146 146 100 100 40 40 20 20],[-70 -70 -8 -8 25 25 -8 -8 -70],'linewidth',1.5,'color','b')
% Pacific 
m_line([146 293 293 270 270 255 255 100 100],[-70 -70 10 10 18 18 65 65 25],'linewidth',1.5,'color','r')
% Atlantic
m_line([293 380 380 255 255],[-70 -70 70 70 60],'linewidth',1.5,'color','g')
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\Indian_Pacific_Atlantic.png'],'-dpng','-r300')
%% longitude mean (zonal mean)
% Indian Ocean
lonmI1 = permute(nanmean(Tempadlf(40:100,83:116,:,:),1),[2 3 4 1]); % 40E-100E, 8S-25N
lonmI2 = permute(nanmean(Tempadlf(20:146,20:82,:,:),1),[2 3 4 1]);
lonmI = cat(1,lonmI2,lonmI1);
% corr with TPI
var = permute(lonmI,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(TPI_filtered,var(:,i,j),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
lonmI1m = permute(nanmean(Temp(40:100,83:116,:,:),1),[2 3 4 1]); % climotology
lonmI2m = permute(nanmean(Temp(20:146,20:82,:,:),1),[2 3 4 1]);
lonmIm = cat(1,lonmI2m,lonmI1m);
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:116,1),-depthData,par11',[-1:0.1:1],'linestyle','none')
hold on
caxis([-1,1])
[ac ah] = contour(latData(20:116,1),-depthData,(nanmean(lonmIm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-70,20]);
set(gca,'XTick',[-60:20:20]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:200:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.IPO.zonal mean temperature in Indian Ocean')
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.corr.Indian.zonalmean.tem.png'],'-dpng','-r300')
%% Pacific 
lonmP1 = permute(nanmean(Tempadlf(146:293,20:82,:,:),1),[2 3 4 1]); % 146E-67W, 70S-8S
lonmP2 = permute(nanmean(Tempadlf(100:293,83:101,:,:),1),[2 3 4 1]); % 100E-67W, 8S-10N
lonmP3 = permute(nanmean(Tempadlf(100:270,102:109,:,:),1),[2 3 4 1]); % 100E-90W, 10N-18N
lonmP4 = permute(nanmean(Tempadlf(100:255,110:156,:,:),1),[2 3 4 1]); % 100E-105W, 18N-65N
lonmP = cat(1,lonmP1,lonmP2,lonmP3,lonmP4);
% corr with TPI
var = permute(lonmP,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(TPI_filtered,var(:,i,j),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
lonmP1m = permute(nanmean(Temp(146:293,20:82,:,:),1),[2 3 4 1]); % climotology
lonmP2m = permute(nanmean(Temp(100:293,83:101,:,:),1),[2 3 4 1]); 
lonmP3m = permute(nanmean(Temp(100:270,102:109,:,:),1),[2 3 4 1]);
lonmP4m = permute(nanmean(Temp(100:255,110:156,:,:),1),[2 3 4 1]); 
lonmPm = cat(1,lonmP1m,lonmP2m,lonmP3m,lonmP4m);
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:156,1),-depthData,par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
[ac ah] = contour(latData(20:156,1),-depthData,(nanmean(lonmPm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:200:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.IPO.zonal mean temperature in Pacific')
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.corr.Pacific.zonalmean.tem.png'],'-dpng','-r300')
%% Atlantic
TemAtl = cat(1,Tempadlf(255:360,:,:,:),Tempadlf(1:20,:,:,:));
lonAtl = [lonData(255:360)-360;lonData(1:20)];
lonmA1 = permute(nanmean(TemAtl(39:end,20:101,:,:),1),[2 3 4 1]); % 67W-20E, -70S-10N
lonmA2 = permute(nanmean(TemAtl(16:end,102:109,:,:),1),[2 3 4 1]); % 90W-20E, 11N-18N 
lonmA3 = permute(nanmean(TemAtl(:,110:161,:,:),1),[2 3 4 1]); % 105W-20E, 19N-70N
lonmA = cat(1,lonmA1,lonmA2,lonmA3);
% corr with TPI
var = permute(lonmA,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(TPI_filtered,var(:,i,j),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
TemAtlm = cat(1,Temp(255:360,:,:,:),Temp(1:20,:,:,:)); % climotology
lonmA1m = permute(nanmean(TemAtlm(39:end,20:101,:,:),1),[2 3 4 1]); 
lonmA2m = permute(nanmean(TemAtlm(16:end,102:109,:,:),1),[2 3 4 1]); 
lonmA3m = permute(nanmean(TemAtlm(:,110:161,:,:),1),[2 3 4 1]); 
lonmAm = cat(1,lonmA1m,lonmA2m,lonmA3m);
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:161,1),-depthData,par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
[ac ah] = contour(latData(20:161,1),-depthData,(nanmean(lonmAm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-70,70]);
set(gca,'XTick',[-60:20:70]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:200:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.IPO.zonal mean temperature in Atlantic')
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.corr.Atlantic.zonalmean.tem.png'],'-dpng','-r300')
%% meridional mean 
meridionalmean(Tempadlf,Temp,latData,[50:111],'40S-20N',TPI_filtered,lonData,depthData);
%%
meridionalmean(Tempadlf,Temp,latData,[20:40],'70S-50S',TPI_filtered,lonData,depthData);
%%
meridionalmean(Tempadlf,Temp,latData,[40:60],'50S-30S',TPI_filtered,lonData,depthData);
%%
meridionalmean(Tempadlf,Temp,latData,[60:116],'30S-25N',TPI_filtered,lonData,depthData);
%%
meridionalmean(Tempadlf,Temp,latData,[116:141],'25N-50N',TPI_filtered,lonData,depthData);
%% different depth
%
diffdepth(Tempadlf,1,TPI_filtered,lonData,latData,'0000')
diffdepth(Tempadlf,7,TPI_filtered,lonData,latData,'0050')
diffdepth(Tempadlf,12,TPI_filtered,lonData,latData,'0100')
diffdepth(Tempadlf,17,TPI_filtered,lonData,latData,'0200')
diffdepth(Tempadlf,23,TPI_filtered,lonData,latData,'0500')
diffdepth(Tempadlf,29,TPI_filtered,lonData,latData,'0800')
diffdepth(Tempadlf,32,TPI_filtered,lonData,latData,'1000')
diffdepth(Tempadlf,37,TPI_filtered,lonData,latData,'1500')
diffdepth(Tempadlf,41,TPI_filtered,lonData,latData,'2000')
%% 100-400m
Tm1 = nanmean(Tempadlf(:,:,12:21,:),3);
diffdepth(Tm1,1,TPI_filtered,lonData,latData,'100-400')
%% zonal mean 140W-20W
lonm1 = permute(nanmean(Tempadlf(220:340,:,1:23,:),1),[2 3 4 1]);
lonm1m = permute(nanmean(Temp(220:340,:,1:23,:),1),[2 3 4 1]); % climotology
% corr with TPI
var = permute(lonm1,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(TPI_filtered,var(:,i,j),0.05); % filtered
    end
end
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData,-depthData(1:23),par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-70,70]);
set(gca,'XTick',[-60:20:70]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.IPO.zonal mean temperature 140-20W')
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.corr.zonalmean.140W-20W.tem.png'],'-dpng','-r300')
%% meridional mean 50S-30N
lats = [40:121]; names = '50S-30N'; index = TPI_filtered;
latm = latmean(Tempadlf,lats,latData);
latmc = latmean(Temp,lats,latData); % climotology
% corr with TPI
var = permute(latm,[3 1 2]);
clear par11 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(index,var(:,i,j),0.05);
    end
end
max(par11,[],'all')
min(par11,[],'all')
%%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(lonData,-depthData,par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
[ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[0,360]);
set(gca,'XTick',[0:60:360]);
set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title(['Corr.IPO.meridional mean temperature (',names,')'])
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.corr.meridionalmean',names,'.tem.png'],'-dpng','-r300')






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
function meridionalmean(Tempadlf,Temp,latData,lats,names,index,lonData,depthData)
    latm = latmean(Tempadlf,lats,latData);
    latmc = latmean(Temp,lats,latData); % climotology
    % corr with TPI
    var = permute(latm,[3 1 2]);
    clear par11 h05
    for i = 1:size(var,2);
        for j = 1:size(var,3);
            [par11(i,j),h05(i,j),t] = corr_eff(index,var(:,i,j),0.05); 
        end
    end
    max(par11,[],'all')
    min(par11,[],'all')
    %
    close all;
    Fig = figure('position',[100 100 800 400]);
    contourf(lonData,-depthData,par11',[-1:0.1:1],'linestyle','none')
    caxis([-1,1])
    hold on
    [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    ch = colorbar;
    set(ch,'Ticks',[-1:0.2:1]);
    set(gca,'XLim',[0,360]);
    set(gca,'XTick',[0:60:360]);
    set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
    set(gca,'YLim',[-2000,0]);
    set(gca,'YTick',[-2000:200:0],'FontSize',14);
    set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
    ylabel('Depth (m)')
    title(['Corr.IPO.meridional mean temperature (',names,')'])
    print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.corr.meridionalmean',names,'.tem.png'],'-dpng','-r300')
end
function diffdepth(Tempadlf,nlevel,index,lonData,latData,names)
    Tsdf = reshape(Tempadlf(:,:,nlevel,:),[360*180,73]);  % sst
    clear par11 h05
    parfor i = 1:360*180;
            [par11(i),h05(i),t] = corr_eff(index,Tsdf(i,:),0.05); 
    end
    par11 = reshape(par11,360,180);
    max(par11,[],'all')
    min(par11,[],'all')
    %
    close all;
    Fig = figure('position',[100 100 600 300]);
    contourVARra(par11,[-1:0.1:1],1,0,360,-90,90,lonData,latData);
    ch = colorbar;
    set(ch,'Ticks',[-1:0.2:1]);
    % testdots(h05,'k',lonData,latData)
    title(['Corr.on.T',' ',names,'m']);
    print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\TPIcorrT',names,'m.png'],'-dpng','-r300')
end

