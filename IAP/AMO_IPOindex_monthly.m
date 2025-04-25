clc,clear,close all;
addpath G:\1_matlab\help;
datadir='G:\data\IAP\temperature\monthly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'CZ16*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE');
lonData = ncread(filetemp,'lon');
latData = ncread(filetemp,'lat');
depthData = ncread(filetemp,'depth_std');
% save('MatFile/lonData.mat','lonData');
% save('MatFile/latData.mat','latData');
% save('MatFile/depthData.mat','depthData');
%% anomaly 
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp(:,:,:,s) = ncread(filename1,'temp'); %读入变量
    netcdf.close(ncid);  %关闭文件    
end
%%
Temp1 = reshape(Temp,[41 360 180 12 984/12]);
Tempa = Temp1 - nanmean(Temp1(:,:,:,:,42:71),5); % remove climatology from 1981-2010 
Ts1 = reshape(Tempa(1,:,:,:,:),[360 180 12*82]);
Ts_globalmean = areamean(Ts1,1:360,1:180,latData); 
Tsd = Ts1-permute(Ts_globalmean,[2 3 1]); % remove global mean SST each month (NOAA NCL scripts)
% IPO index
[ts1_zs ts1] = areamean(Tsd,140:215,116:136,latData); % 25N-45N,140E-145W
[ts2_zs ts2] = areamean(Tsd,170:270,80:101,latData); % 10S-10N,170E-90W
[ts3_zs ts3] = areamean(Tsd,150:200,40:75,latData); % 50S-15S,150E-160W
TPI1 = ts2-(ts1+ts3)/2; % unfiltered IPO index
dT = 1; % interval
% cf = 1/13;
cf = 1/(13*12);
TPI2 = lanczosfilter(TPI1,dT,cf,[],'low'); % 13 year filtered IPO index

figure(1)
plot(TPI1)
hold on
plot(TPI2)
plot(zeros(1,12*82))
set(gca,'XLim',[1,12*82]);
set(gca,'XTick',[1:10*12:12*82]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',12);
set(gca,'YLim',[-4,4]);
set(gca,'YTick',[-4:1:4]);
xlabel('YEAR','FontSize',12),ylabel('TPI','FontSize',12);
title('TPI');

figure(2)
plot(TPI2)
hold on
plot(zeros(1,12*82))
% plot([11:82],TPI2ersst,'k')
set(gca,'XLim',[1,12*82]);
set(gca,'XTick',[1:10*12:12*82]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',12);
set(gca,'YLim',[-1,1]);
set(gca,'YTick',[-1:0.2:1]);
%% regress on sst
for i = 1:size(Tsd,1);
    for j = 1:size(Tsd,2);
        [par11(i,j),h05(i,j),t] = reg1_ttest(TPI1,Ts1(i,j,:),0.05); % sat
    end
end
max(par11,[],'all')
min(par11,[],'all')
%%
close all;
Fig = figure('position',[100 100 600 300]);
contourVARra(par11,[-1.6:0.2:1.6],1.2,0,360,-90,90,lonData,latData);
ch = colorbar;
set(ch,'Ticks',[-1.2:0.4:1.2]);
title('Reg.on.SST');
print(Fig,['G:\博后入站\TPIregSST.png'],'-dpng','-r600')
%% TPI from NOAA
%TPI ERSST V5 from NOAA
load('G:\data\NCEP\indices\TPI(from NOAA ERSST V5)unfiltered 1981-2010climo.txt');
TPIersst_yr = mean(TPI_from_NOAA_ERSST_V5_unfiltered_1981_2010climo(:,2:13),2); % year mean
dT = 1; % interval
cf = 1/13;
TPIersst_yrfilter = lanczosfilter(TPIersst_yr,dT,cf,[],'low'); % 13 year filtered IPO indexfigure(4)
TPIersst_yrint = interp(TPIersst_yrfilter,12);
% TPI HADISST1.1 from NOAA
load('G:\data\NCEP\indices\TPI(from HADISST1.1)unfiltered Standard PSL Format.txt');
TPIhadisst_yr = mean(TPI_from_HADISST1_1_unfiltered_Standard_PSL_Format(:,2:13),2);
dT = 1; % interval
cf = 1/13;
TPIhadisst_yrfilter = lanczosfilter(TPIhadisst_yr,dT,cf,[],'low'); % 13 year filtered IPO indexfigure(4)
TPIhadisst_yrint = interp(TPIhadisst_yrfilter,12);
% TPI COBE from NOAA 
load('G:\data\NCEP\indices\TPI(from COBE)filtered Standard PSL Format.txt');
TPIcobe_yr = mean(TPI_from_COBE_filtered_Standard_PSL_Format(:,2:13),2);
dT = 1; % interval
cf = 1/13;
TPIcobe_yrfilter = lanczosfilter(TPIcobe_yr,dT,cf,[],'low'); % 13 year filtered IPO indexfigure(4)
TPIcobe_yrint = interp(TPIcobe_yrfilter,12);
%% HADISST calculate
ncid=netcdf.open('G:\data\hadley\HadISST_sst.nc','NOWRITE');
ncdisp('G:\data\hadley\HadISST_sst.nc');
latDataH = ncread('G:\data\hadley\HadISST_sst.nc','latitude');
lonDataH = ncread('G:\data\hadley\HadISST_sst.nc','longitude');
hadisst = ncread('G:\data\hadley\HadISST_sst.nc','sst');
hadisst(find(hadisst == -1000)) = nan;
Temphadi1 = reshape(hadisst,[360 180 12 1824/12]);
Temphadi2 = Temphadi1 - nanmean(Temphadi1(:,:,:,112:141),4); % remove climatology from 1981-2010 
% detrend 
Tshadi = reshape(Temphadi2,[360 180 12*152]);
Tshadi_globalmean = areamean(Tshadi,1:360,1:180,latDataH); 
Tsdhadi = Tshadi-permute(Tshadi_globalmean,[2 3 1]); % remove global mean SST each month
% IPO index
Tsrvs1 = cat(1,Tsdhadi(181:360,:,:),Tsdhadi(1:180,:,:)); % note different longitude and latitude; reverse same as IAPocean
Tsrvs = fliplr(Tsrvs1);
[ts1_zsh ts1h] = areamean(Tsrvs,140:215,116:136,latDataH); % 25N-45N,140E-145W
[ts2_zsh ts2h] = areamean(Tsrvs,170:270,80:101,latDataH); % 10S-10N,170E-90W
[ts3_zsh ts3h] = areamean(Tsrvs,150:200,40:75,latDataH); % 50S-15S,150E-160W
TPI1h = ts2h-(ts1h+ts3h)/2; % unfiltered IPO index
dT = 1; % interval
% cf = 1/13;
cf = 1/(13*12);
TPI2h = lanczosfilter(TPI1h,dT,cf,[],'low'); % 13 year filtered IPO index
%%
close all;
Fig = figure('position',[100 100 1200 400])
plot(TPIersst_yrint,'color',[67 189 82]/255,'linewidth',1.5)
hold on
plot([12*16+1:2016],TPIhadisst_yrint,'b','linewidth',1.5);
plot([12*37+1:2016],TPIcobe_yrint,'color',[118 240 240]/255,'linewidth',1.5);
plot([12*16+1:2016],TPI2h,'color',[251 93 5]/255,'linewidth',2.5);
plot([12*86+1:2016],TPI2,'k','linewidth',2.5);
plot(zeros(1,2016),'k')
set(gca,'XLim',[1,2016]);
set(gca,'XTick',[6*12:10*12:2016]);
set(gca,'XTickLabel',[1860:10:2020],'FontSize',14);
set(gca,'YLim',[-1,1]);
set(gca,'YTick',[-1:0.2:1],'FontSize',14);
legend('ERSST(NOAA)','HADISST1.1(NOAA)','COBE(NOAA)','HADISST','IAPoceanT','Location','south','Orientation','horizontal')
legend('boxoff');
title('TPI','fontsize',14);
ylabel('Temperature Anomaly (^oC)','fontsize',14);
print(Fig,['G:\博后入站\TPIs.png'],'-dpng','-r600')
%% AMO index
Ts1d = reshape(Ts1,[360 180 12 82]);
for i = 1:360
    for j = 1:180
        for m = 1:12
            Tsdetrend(i,j,m,:) = permute(detrend(permute(Ts1d(i,j,m,:),[4 3 2 1])),[4 3 2 1]); % detrend
        end
    end
end
%%
Tsdetrendr = reshape(Tsdetrend,[360 180 12*82]);
[AMO_zs AMO] = areamean(Tsdetrend,290:360,91:161,latData); % 70W-0,0-70N
dT = 1; % interval
cf = 1/(8*12);
AMOf = lanczosfilter(AMO,dT,cf,[],'low'); % 8 year filtered AMO index
close all;
figure(1)
plot(AMO)
hold on
plot(AMOf)
plot(zeros(1,12*82))
set(gca,'XLim',[1,12*82]);
set(gca,'XTick',[1:10*12:12*82]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',12);
set(gca,'YLim',[-4,4]);
set(gca,'YTick',[-4:1:4]);
xlabel('YEAR','FontSize',12),ylabel('TPI','FontSize',12);
title('AMO');

figure(2)
plot(AMOf)
hold on
plot(zeros(1,12*82))
% plot([11:82],TPI2ersst,'k')
set(gca,'XLim',[1,12*82]);
set(gca,'XTick',[1:10*12:12*82]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',12);
set(gca,'YLim',[-1,1]);
set(gca,'YTick',[-1:0.2:1]);
%% regress on sst
for i = 1:size(Tsdetrend,1);
    for j = 1:size(Tsdetrend,2);
        [par11(i,j),h05(i,j),t] = reg1_ttest(AMO,Tsdetrend(i,j,:),0.05); % sat
    end
end
max(par11,[],'all')
min(par11,[],'all')
%%
close all;
Fig = figure('position',[100 100 600 300]);
contourVARra(par11,[-1.6:0.2:1.6],1.2,0,360,-90,90,lonData,latData);
ch = colorbar;
set(ch,'Ticks',[-1.2:0.4:1.2]);
title('Reg.on.SST');
print(Fig,['G:\博后入站\AMOregSST.png'],'-dpng','-r600')
%% HadIsst AMO index
for i = 1:360
    for j = 1:180
        for m = 1:12
            Tsdetrend_hadi(i,j,m,:) = permute(detrend(permute(Temphadi2(i,j,m,:),[4 3 2 1])),[4 3 2 1]); % detrend each month
        end
    end
end
%%
Tsdrvs1 = cat(1,Tsdetrend_hadi(181:360,:,:),Tsdetrend_hadi(1:180,:,:)); % note different longitude and latitude; reverse same as IAPocean
Tsdrvs = fliplr(Tsdrvs1);
[AMOh_zs AMOh] = areamean(Tsdrvs,290:360,91:161,flipud(latData)); % 70W-0,0-70N
dT = 1; % interval
cf = 1/(8*12);
AMOhf = lanczosfilter(AMOh,dT,cf,[],'low'); % 13 year filtered AMO index
%%
%AMO(short) Kaplan SST from NOAA
load('G:\data\NCEP\indices\AMO(unsmoothed).txt');
AMOn_yr = mean(AMO_unsmoothed_(:,2:13),2); % year mean
dT = 1; % interval
cf = 1/8;
AMOn_yrfilter = lanczosfilter(AMOn_yr,dT,cf,[],'low'); % 13 year filtered IPO indexfigure(4)
AMOn_yrint = interp(AMOn_yrfilter,12);
close all;
Fig = figure('position',[100 100 1200 400]);
plot([12*78+1:1824],AMOn_yrint,'color',[251 93 5]/255,'linewidth',2.5);
hold on
plot([1:1824],AMOhf,'color',[118 240 240]/255,'linewidth',2.5);
plot([12*70+1:1824],AMOf,'k','linewidth',2.5);
plot(zeros(1,1824),'k')
set(gca,'XLim',[1,1824]);
set(gca,'XTick',[1:10*12:1824]);
set(gca,'XTickLabel',[1870:10:2020],'FontSize',14);
set(gca,'YLim',[-0.6,0.6]);
set(gca,'YTick',[-0.6:0.2:0.6],'FontSize',14);
legend('KaplanSST(NOAA)','HADISST','IAPoceanT','Location','south','Orientation','horizontal')
legend('boxoff');
title('AMO index','fontsize',14);
ylabel('Temperature Anomaly (^oC)','fontsize',14);
print(Fig,['G:\博后入站\AMOs.png'],'-dpng','-r600')


function contourVARra(var,val,ctr,lon1,lon2,laTempa,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[laTempa lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re3_16_02.mat');
    bl_re3_16_02([1,8,9,16],:) = [];
    colormap(bl_re3_16_02);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
end
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end
function rectbox(lonW,lonE,latS,latN,colorset)
    m_line([lonW lonE lonE lonW lonW],[latS latS latN latN latS],'linewidth',1.5,'color',colorset)
end
function testdots(h,clor,lonData,latData)
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',4,'color',clor); % F-test dots
end
function [ts] = lonmean(var,lats,latData)
% var is lon*lat*depth*time
% ts is lon*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,2)./nansum(cos(latData(lats)'/180*pi).*var2,2),[1 3 4 2]);
end