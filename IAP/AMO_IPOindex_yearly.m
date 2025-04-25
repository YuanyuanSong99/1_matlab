clc,clear,close all;
addpath D:\1_matlab\help;
datadir='D:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% save('MatFile/lonData.mat','lonData');
% save('MatFile/latData.mat','latData');
% save('MatFile/depthData.mat','depthData');
% anomaly 
k=length(filelist);

for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
%%
Tempa = double(Temp - nanmean(Temp(:,:,:,42:71),4)); % remove climatology from 1981-2010 
%% Ts1 = reshape(Tempa(1,:,:,:,:),[360 180 12*82]);
% Ts_globalmean = areamean(Ts1,1:360,1:180,latData); 
% Tsd = Ts1-permute(Ts_globalmean,[2 3 1]); % remove global mean SST each month (NOAA NCL scripts)
% detrend and filter 13 years each grid
Tempa1 = double(permute(Tempa,[4 1 2 3]));
dT = 1; % interval
cf = 1/8;
Tempa1 = reshape(Tempa1,[81,360*180*41]);
parfor i = 1:size(Tempa1,2);
    Tempadl(:,i) = detrend(Tempa1(:,i)); % linear detrend
    Tempadlf(:,i) = lanczosfilter(Tempadl(:,i),dT,cf,[],'low'); % 8 year filtered
end
Tempadl = permute(reshape(Tempadl,[81,360,180,41]),[2 3 4 1]);  % linear detrended      
save('MatFile/Tempadl.mat','Tempadl');
Tempadlf = permute(reshape(Tempadlf,[81,360,180,41]),[2 3 4 1]);   % linear detrended and 8-year filtered  
save('MatFile/Tempadlf.mat','Tempadlf');
Tsdf = permute(Tempadlf(:,:,1,:),[1 2 4 3]);  % sst
%% nonlinear detrend and filter 8 years each grid
addpath D:\1_matlab\help\EEMD_code\;
Tempa1 = double(permute(Tempa,[4 1 2 3]));
Tempa1 = reshape(Tempa1,[81,360*180*41]);
dT = 1; cf = 1/8; % 8-yr filter
parfor i = 1:2656800
    i
    if sum(ismissing(Tempa1(:,i))) == 0
    Tempeemd = eemd(Tempa1(:,i),0.2,1000); % EEMD
    Tempadnl(:,i) = Tempeemd(:,7); % nonlinear detrend
    Tempadnlf(:,i) = lanczosfilter(Tempadnl(:,i),dT,cf,[],'low'); % 8 year filtered
    end
end
Tempadnl = permute(reshape(Tempadnl,[81,360,180,41]),[2 3 4 1]);  % nonlinear detrended      
save('MatFile/Tempadnl.mat','Tempadnl');
Tempadnlf = permute(reshape(Tempadnlf,[81,360,180,41]),[2 3 4 1]);   % nonlinear detrended and 8-year filtered  
save('MatFile/Tempadnlf.mat','Tempadnlf');

%% IPO index
Tsd = permute(Tempadl(:,:,1,:),[1 2 4 3]);  % sst
[ts1_zs ts1] = areamean(Tsd,140:215,116:136,latData); % 25N-45N,140E-145W
[ts2_zs ts2] = areamean(Tsd,170:270,80:101,latData); % 10S-10N,170E-90W
[ts3_zs ts3] = areamean(Tsd,150:200,40:75,latData); % 50S-15S,150E-160W
TPI1 = ts2-(ts1+ts3)/2; % unfiltered IPO index///
dT = 1; % interval
cf = 1/8;
TPI2 = lanczosfilter(TPI1,dT,cf,[],'low'); % 8 year filtered IPO index
TPI_13yr = lanczosfilter(TPI1,dT,1/13,[],'low'); % 13 year filtered IPO index
%% 

figure(1)
plot(TPI1)
hold on
plot(TPI2)
plot(zeros(1,81))
set(gca,'XLim',[1,82]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',12);
set(gca,'YLim',[-4,4]);
set(gca,'YTick',[-4:1:4]);
xlabel('YEAR','FontSize',12),ylabel('TPI','FontSize',12);
title('TPI');

% close all;
Fig = figure('position',[100 100 800 400])
plot(TPI2,'k','linewidth',1.5);
hold on
plot(TPI_13yr,'r','linewidth',1.5)
plot(zeros(1,81),'k','linewidth',1.5)
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[-1,1]);
set(gca,'YTick',[-1:0.2:1],'FontSize',14);
legend('TPI-8yr','TPI-13yr','Location','south','Orientation','horizontal')
legend('boxoff');
title('TPI','fontsize',14);
ylabel(' (^oC)','fontsize',14);
% print(Fig,['D:\figures\IAP\Yearly\TPIs_8yr_13yr.png'],'-dpng','-r300')

%% correlation with sst
for i = 1:size(Tsd,1);
    for j = 1:size(Tsd,2);
        [par11(i,j),h05(i,j),t] = corr_eff(TPI1,Tsd(i,j,:),0.05); % unfiltered
        [par12(i,j),h052(i,j),t] = corr_eff(TPI2,Tsdf(i,j,:),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
%%
close all;
Fig = figure('position',[100 100 600 300]);
contourVARra(par11,[-1:0.1:1],1,0,360,-90,90,lonData,latData);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
% testdots(h05,'k',lonData,latData)
title('Corr.on.SST');
% print(Fig,['D:\figures\IAP\Yearly\TPIcorrSSTunfiltered1.png'],'-dpng','-r300')
%%
close all;
Fig = figure('position',[100 100 600 300]);
contourVARra(par12,[-1:0.1:1],1,0,360,-90,90,lonData,latData);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
% testdots(h05,'k',lonData,latData)
title('Corr.on.SST');
print(Fig,['D:\figures\IAP\Yearly\TPIcorrSSTfiltered1.png'],'-dpng','-r300')
%%
close all;
Fig = figure('position',[100 100 600 300]);
map1 = par11;
map1(find(h05 == 0)) = nan;
contourVARra(map1,[-1:0.1:1],1,0,360,-90,90,lonData,latData);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
title('Corr.on.SST');
print(Fig,['D:\figures\IAP\Yearly\TPIcorrSSTunfiltered2.png'],'-dpng','-r300')
%%
close all;
Fig = figure('position',[100 100 600 300]);
map1 = par12;
map1(find(h052 == 0)) = nan;
contourVARra(map1,[-1:0.1:1],1,0,360,-90,90,lonData,latData);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
title('Corr.on.SST');
print(Fig,['D:\figures\IAP\Yearly\TPIcorrSSTfiltered2.png'],'-dpng','-r300')
%% TPI from NOAA
%TPI ERSST V5 from NOAA
load('D:\data\NCEP\indices\TPI(from NOAA ERSST V5)unfiltered 1981-2010climo.txt');
TPIersst_yr = mean(TPI_from_NOAA_ERSST_V5_unfiltered_1981_2010climo(:,2:13),2); % year mean
dT = 1; % interval
cf = 1/13;
TPIersst_yrfilter = lanczosfilter(TPIersst_yr,dT,cf,[],'low'); % 13 year filtered IPO indexfigure(4)
% TPI HADISST1.1 from NOAA
load('D:\data\NCEP\indices\TPI(from HADISST1.1)unfiltered Standard PSL Format.txt');
TPIhadisst_yr = mean(TPI_from_HADISST1_1_unfiltered_Standard_PSL_Format(:,2:13),2);
dT = 1; % interval
cf = 1/13;
TPIhadisst_yrfilter = lanczosfilter(TPIhadisst_yr,dT,cf,[],'low'); % 13 year filtered IPO indexfigure(4)
% TPI COBE from NOAA 
load('D:\data\NCEP\indices\TPI(from COBE)filtered Standard PSL Format.txt');
TPIcobe_yr = mean(TPI_from_COBE_filtered_Standard_PSL_Format(:,2:13),2);
dT = 1; % interval
cf = 1/13;
TPIcobe_yrfilter = lanczosfilter(TPIcobe_yr,dT,cf,[],'low'); % 13 year filtered IPO indexfigure(4)
%% HADISST calculate 1870-2021 
ncid=netcdf.open('D:\data\hadley\HadISST_sst.nc','NOWRITE');
ncdisp('D:\data\hadley\HadISST_sst.nc');
latDataH = ncread('D:\data\hadley\HadISST_sst.nc','latitude');
lonDataH = ncread('D:\data\hadley\HadISST_sst.nc','longitude');
hadisst = ncread('D:\data\hadley\HadISST_sst.nc','sst');
hadisst(find(hadisst == -1000)) = nan;
Temphadi1 = reshape(hadisst,[360 180 12 1836/12]);
Temphadi2 = permute(nanmean(Temphadi1,3),[1 2 4 3]);
Temphadi3 = Temphadi2 - nanmean(Temphadi2(:,:,112:141),3); % remove climatology from 1981-2010 
% detrend 
Temphadi4 = permute(Temphadi3,[3 1 2]);
for i = 1:size(Temphadi4,2);
    for j = 1:size(Temphadi4,3);
        Tsdhadi1(:,i,j) = detrend(Temphadi4(:,i,j)); 
    end
end
Tsdhadi = permute(Tsdhadi1,[2 3 1]);
% IPO index
Tsrvs1 = cat(1,Temphadi3(181:360,:,:),Temphadi3(1:180,:,:)); % note different longitude and latitude; reverse same as IAPocean
Tsrvs = fliplr(Tsrvs1);
[ts1_zsh ts1h] = areamean(Tsrvs,140:215,116:136,latData); % 25N-45N,140E-145W
[ts2_zsh ts2h] = areamean(Tsrvs,170:270,80:101,latData); % 10S-10N,170E-90W
[ts3_zsh ts3h] = areamean(Tsrvs,150:200,40:75,latData); % 50S-15S,150E-160W
TPI1h = ts2h-(ts1h+ts3h)/2; % unfiltered IPO index
dT = 1; % interval
cf = 1/13;
TPI_filtered13_hadley = lanczosfilter(TPI1h,dT,cf,[],'low'); % 13 year filtered IPO index
cf = 1/8;
TPI_filtered8_hadley = lanczosfilter(TPI1h,dT,cf,[],'low'); % 8 year filtered IPO index
% save('MatFile/TPI_filtered8_hadley.mat','TPI_filtered8_hadley');
save('MatFile/TPI_filtered13_hadley.mat','TPI_filtered13_hadley');
%%
close all;
Fig = figure('position',[100 100 1200 400])
plot(TPIersst_yrfilter(8:end-7),'color',[67 189 82]/255,'linewidth',1.5)
hold on
plot([16+1:154],TPIhadisst_yrfilter(8:end-7),'b','linewidth',1.5);
plot([37+1:154],TPIcobe_yrfilter(8:end-7),'color',[118 240 240]/255,'linewidth',1.5);
% plot([16+1:161],TPI2h,'color',[251 93 5]/255,'linewidth',2.5); % Hadisst
plot([16+1:154],TPI2h(8:end-8),'k','linewidth',2.5);
plot(zeros(1,154),'k')
set(gca,'XLim',[1,154]);
set(gca,'XTick',[10:10:154]);
set(gca,'XTickLabel',[1870:10:2020],'FontSize',14);
set(gca,'YLim',[-1,1]);
set(gca,'YTick',[-1:0.2:1],'FontSize',14);
legend('ERSST(NOAA)','HADISST1.1(NOAA)','COBE(NOAA)','HADISSTcal','Location','south','Orientation','horizontal')
legend('boxoff');
title('TPI','fontsize',14);
ylabel('Temperature Anomaly (^oC)','fontsize',14);
% print(Fig,['D:\figures\IAP\Yearly\TPIs.png'],'-dpng','-r300')
%% AMO index (hadley sst)
%
[AMO_zs AMO] = areamean(Tsdhadi,110:181,20:90,latDataH); % 70W-0,0-70N
dT = 1; % interval
cf = 1/8;
AMOf8_hadley = lanczosfilter(AMO,dT,cf,[],'low'); % 8 year filtered AMO index
save('MatFile/AMOf8_hadley.mat','AMOf8_hadley');

%%
close all;
Fig = figure('position',[2700 100 800 400])
plot(AMO,'--','color',[0,46,166]/255,'LineWidth',1.5)
hold on
plot(AMOf8,'color',[0,46,166]/255,'LineWidth',1.5)
plot(zeros(1,154),'k','LineWidth',1)
set(gca,'XLim',[1,154]);
set(gca,'XTick',[1:10:154]);
set(gca,'XTickLabel',[1870:10:2020],'FontSize',14);
set(gca,'YLim',[-0.6,0.6]);
set(gca,'YTick',[-0.6:0.2:0.6]);
xlabel('YEAR','FontSize',12),ylabel('AMO index','FontSize',14);
legend('unfiltered','filtered','Location','north','Orientation','horizontal')
legend('boxoff');
% print(Fig,['D:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\AMOindex.png'],'-dpng','-r300')
%%
figure(2)
plot(AMOf8)
hold on
plot(zeros(1,81))
% plot([11:82],TPI2ersst,'k')
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',12);
set(gca,'YLim',[-0.6,0.6]);
set(gca,'YTick',[-0.6:0.2:0.6]);
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
print(Fig,['D:\博后入站\AMOregSST.png'],'-dpng','-r300')
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
load('D:\data\NCEP\indices\AMO(unsmoothed).txt');
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
plot([12*70+1:1824],AMOf8,'k','linewidth',2.5);
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
% print(Fig,['D:\博后入站\AMOs.png'],'-dpng','-r300')


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
