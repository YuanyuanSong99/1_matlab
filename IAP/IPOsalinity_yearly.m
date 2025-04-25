clc,clear,close all;
addpath E:\1_matlab\help;
datadir='E:\data\IAP\salinity\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
fileSalt=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(fileSalt,'/lon');
latData = h5read(fileSalt,'/lat');
depthData = h5read(fileSalt,'/level');
% anomaly 
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Salt(:,:,:,s) = ncread(filename1,'absolute salinity'); %读入变量
end
%%
Salta = Salt - nanmean(Salt(:,:,:,42:71),4); % remove climatology from 1981-2010 
% Ts1 = reshape(Salta(1,:,:,:,:),[360 180 12*82]);
% Ts_globalmean = areamean(Ts1,1:360,1:180,latData); 
% Tsd = Ts1-permute(Ts_globalmean,[2 3 1]); % remove global mean SST each month (NOAA NCL scripts)
% detrend and filter 13 years each grid
Salta1 = double(permute(Salta,[4 1 2 3]));
dT = 1; % interval
cf = 1/13;
for i = 1:size(Salta1,2);
    for j = 1:size(Salta1,3);
        for k = 1:size(Salta1,4);
            Saltad(:,i,j,k) = detrend(Salta1(:,i,j,k));
            Saltadf(:,i,j,k) = lanczosfilter(Saltad(:,i,j,k),dT,cf,[],'low'); % 13 year filtered 
        end
    end
end
Saltad = permute(Saltad,[2 3 4 1]);  % detrended      
save('MatFile/Saltad.mat','Saltad');
Saltadf = permute(Saltadf,[2 3 4 1]);   % detrended and 13-year filtered  
save('MatFile/Saltadf.mat','Saltadf');
sssdf = permute(Saltadf(:,:,1,:),[1 2 4 3]);  % sss
%%
addpath E:\1_matlab\help;
load("MatFile\Saltadf.mat"); % detrended and 13-yr filtered temperature
load("MatFile\TPI_filtered.mat");
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
%% longitude mean (zonal mean)
% Indian Ocean
lonmI1 = permute(nanmean(Saltadf(40:100,83:116,:,:),1),[2 3 4 1]);
lonmI2 = permute(nanmean(Saltadf(20:146,20:82,:,:),1),[2 3 4 1]);
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
lonmI1m = permute(nanmean(Salt(40:100,83:116,:,:),1),[2 3 4 1]); % climotology
lonmI2m = permute(nanmean(Salt(20:146,20:82,:,:),1),[2 3 4 1]);
lonmIm = cat(1,lonmI2m,lonmI1m);
%%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:116,1),-depthData,par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
[ac ah] = contour(latData(20:116,1),-depthData,(nanmean(lonmIm,3))',14,'ShowText','on','LevelList',[32:0.2:38],'LineColor','k','LabelSpacing',240); 
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
title('Corr.IPO.zonal mean salinity in Indian Ocean')
print(Fig,['E:\figures\IAP\Yearly\IPO.corr.Indian.zonalmean.salt.png'],'-dpng','-r600')
%% Pacific 
lonmP1 = permute(nanmean(Saltadf(146:293,20:82,:,:),1),[2 3 4 1]); % 146E-67W, 70S-8S
lonmP2 = permute(nanmean(Saltadf(100:293,83:101,:,:),1),[2 3 4 1]); % 100E-67W, 8S-10N
lonmP3 = permute(nanmean(Saltadf(100:270,102:109,:,:),1),[2 3 4 1]); % 100E-90W, 10N-18N
lonmP4 = permute(nanmean(Saltadf(100:255,110:156,:,:),1),[2 3 4 1]); % 100E-105W, 18N-65N
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
lonmP1m = permute(nanmean(Salt(146:293,20:82,:,:),1),[2 3 4 1]); % climotology
lonmP2m = permute(nanmean(Salt(100:293,83:101,:,:),1),[2 3 4 1]); % 
lonmP3m = permute(nanmean(Salt(100:270,102:109,:,:),1),[2 3 4 1]); % 
lonmP4m = permute(nanmean(Salt(100:255,110:156,:,:),1),[2 3 4 1]); % 
lonmPm = cat(1,lonmP1m,lonmP2m,lonmP3m,lonmP4m);
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:156,1),-depthData,par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
[ac ah] = contour(latData(20:156,1),-depthData,(nanmean(lonmPm,3))',14,'ShowText','on','LevelList',[32:0.2:38],'LineColor','k','LabelSpacing',240); 
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
title('Corr.IPO.zonal mean salinity in Pacific')
print(Fig,['E:\figures\IAP\Yearly\IPO.corr.Pacific.zonalmean.salt.png'],'-dpng','-r600')
%% Atlantic
SaltAtl = cat(1,Saltadf(255:360,:,:,:),Saltadf(1:20,:,:,:));
lonAtl = [lonData(255:360)-360;lonData(1:20)];
lonmA1 = permute(nanmean(SaltAtl(39:end,20:101,:,:),1),[2 3 4 1]); % 67W-20E, -70S-10N
lonmA2 = permute(nanmean(SaltAtl(16:end,102:109,:,:),1),[2 3 4 1]); % 90W-20E, 11N-18N 
lonmA3 = permute(nanmean(SaltAtl(:,110:161,:,:),1),[2 3 4 1]); % 105W-20E, 19N-70N
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
SaltAtlm = cat(1,Salt(255:360,:,:,:),Salt(1:20,:,:,:));
lonmA1m = permute(nanmean(SaltAtlm(39:end,20:101,:,:),1),[2 3 4 1]); % 67W-20E, -70S-10N
lonmA2m = permute(nanmean(SaltAtlm(16:end,102:109,:,:),1),[2 3 4 1]); % 90W-20E, 11N-18N 
lonmA3m = permute(nanmean(SaltAtlm(:,110:161,:,:),1),[2 3 4 1]); % 105W-20E, 19N-70N
lonmAm = cat(1,lonmA1m,lonmA2m,lonmA3m);
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:161,1),-depthData,par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
[ac ah] = contour(latData(20:161,1),-depthData,(nanmean(lonmAm,3))',14,'ShowText','on','LevelList',[32:0.2:38],'LineColor','k','LabelSpacing',240); 
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
title('Corr.IPO.zonal mean salinity in Atlantic')
print(Fig,['E:\figures\IAP\Yearly\IPO.corr.Atlantic.zonalmean.salt.png'],'-dpng','-r600')
%% meridional mean 
%
meridionalmean(Saltadf,Salt,latData,[50:111],'40S-20N',TPI_filtered,lonData,depthData);
%
meridionalmean(Saltadf,Salt,latData,[20:40],'70S-50S',TPI_filtered,lonData,depthData);
%
meridionalmean(Saltadf,Salt,latData,[40:60],'50S-30S',TPI_filtered,lonData,depthData);
%
meridionalmean(Saltadf,Salt,latData,[60:116],'30S-25N',TPI_filtered,lonData,depthData);
%
meridionalmean(Saltadf,Salt,latData,[116:141],'25N-50N',TPI_filtered,lonData,depthData);
%% different depth
%
diffdepth(Saltadf,1,TPI_filtered,lonData,latData,'0000')
%
diffdepth(Saltadf,7,TPI_filtered,lonData,latData,'0050')
diffdepth(Saltadf,12,TPI_filtered,lonData,latData,'0100')
diffdepth(Saltadf,17,TPI_filtered,lonData,latData,'0200')
diffdepth(Saltadf,23,TPI_filtered,lonData,latData,'0500')
diffdepth(Saltadf,29,TPI_filtered,lonData,latData,'0800')
diffdepth(Saltadf,32,TPI_filtered,lonData,latData,'1000')
diffdepth(Saltadf,37,TPI_filtered,lonData,latData,'1500')
diffdepth(Saltadf,41,TPI_filtered,lonData,latData,'2000')




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
function meridionalmean(Tempadf,Temp,latData,lats,names,index,lonData,depthData)
    latm = latmean(Tempadf,lats,latData);
    latmc = latmean(Temp,lats,latData); % climotology
    % corr with TPI
    var = permute(latm,[3 1 2]);
    clear par11 h05
    for i = 1:size(var,2);
        for j = 1:size(var,3);
            [par11(i,j),h05(i,j),t] = corr_eff(index,var(:,i,j),0.05); % filtered
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
    [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[30:0.2:40],'LineColor','k','LabelSpacing',240);
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
    title(['Corr.IPO.meridional mean salinity (',names,')'])
    print(Fig,['E:\figures\IAP\Yearly\IPO.corr.meridionalmean',names,'.salt.png'],'-dpng','-r600')
end
function diffdepth(Tempadf,nlevel,index,lonData,latData,names)
    Tsdf = reshape(Tempadf(:,:,nlevel,:),[360*180,81]);  % sss
    clear par11 h05
    parfor i = 1:360*180;
            [par11(i),h05(i),t] = corr_eff(index,Tsdf(i,:),0.05); % filtered
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
    title(['Corr.on.S',' ',names,'m']);
    print(Fig,['E:\figures\IAP\Yearly\TPIcorrS',names,'m.png'],'-dpng','-r600')
end



