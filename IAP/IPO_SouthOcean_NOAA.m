% Trend 0-700m Tem
clc,clear,close all;
addpath G:\1_matlab\help;
datadir='G:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% anomaly 
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
Tempa = double(Temp - nanmean(Temp(:,:,:,42:71),4)); % remove climatology from 1981-2010 
sst_iap = permute(Tempa(:,:,1,:),[1 2 4 3]);
%%
lats = 1:55;
dweit = depthData(2:27)-depthData(1:26); % depth weight
Tsubraw = permute(nansum(Tempa(:,lats,1,:)+Tempa(:,lats,2:27,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m
startyr = 1940;
endyr = 2015;
var = permute(Tsubraw(:,:,startyr-1939:endyr-1939),[3 1 2]);
x = [1:size(var,1)]';
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.05); % t test trend
max(trd,[],'all')
min(trd,[],'all')
close all;
ftsz = 12; ticks = .2;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd*10,[-ticks*5:0.01:ticks*5],ticks,lonData,latData(lats),12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K decade^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots(h0,'k',lonData,latData(lats));
print(Fig,['G:\figures\IAP\Yearly\20230427_IPO_SouthernOcean_45_60S\Tem0-700m_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% DI raw
lats = [30:45];
Tsub_r_0 = cat(1,Tsubraw(150:300,:,:),Tsubraw(301:360,:,:),Tsubraw(1:149,:,:));
[spacz_0 spac_0] = areamean(Tsub_r_0,1:151,lats,latData); 
[siaz_0 sia_0] = areamean(Tsub_r_0,152:360,lats,latData); 
DI_0 = zscore(spac_0-sia_0);
%%
P1val = mean(DI_0(1979-1939:1993-1939));
P2val = mean(DI_0(1994-1939:2012-1939));
close all;
Fig = figure('position',[700 100 800 400]);
plot(DI_0(1979-1939:2012-1939),'-','color','k','LineWidth',2)
hold on
plot(zeros(1,34),'k','LineWidth',1)
plot([1:15],zeros(15)+P1val,'-','color','r','LineWidth',2)
plot([16:34],zeros(19)+P2val,'-','color','b','LineWidth',2)
set(gca,'XLim',[1,34]);
set(gca,'XTick',[2:10:34]);
set(gca,'XTickLabel',[1980:10:2010],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_45_60S\DI_raw_1979_2012.png'],'-dpng','-r300')
%%
lats = 1:55;
com1 = nanmean(Tsubraw(:,:,1979-1939:1993-1939),3);
com2 = nanmean(Tsubraw(:,:,1994-1939:2012-1939),3);
close all;
ftsz = 12; ticks = 0.25;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(com1,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_45_60S\Com_raw_1979_1993.png'],'-dpng','-r300')

Fig = figure('position',[610 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(com2,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k');
set_colorbar([0.83 0.08 0.03 0.88],'K',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_45_60S\Com_raw_1994_2012.png'],'-dpng','-r300')
%% load Temperature
% clc,clear,close all;
addpath G:\1_matlab\help;
addpath G:\1_matlab\help\seawater\;
load('MatFile/Tempadlf');
Tempadlf_long = Tempadlf; % 1940-2020
Tempadlf(:,:,:,1:4) = [];
Tempadlf(:,:,:,end-3:end) = []; % 1944-2016
load("MatFile\TPI_filtered.mat");
TPI_filtered(1:4) = []; TPI_filtered(end-3:end) = [];
TPIfz = zscore(TPI_filtered);
load("MatFile\AMOf8.mat");
AMOf8(1:4) = []; AMOf8(end-3:end) =[];
AMOfz = zscore(AMOf8);
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
[nlon nlat nlev nyr] = size(Tempadlf);
%%  load u v 
filename = 'G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\10mSI-MO\uwnd.10m.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
uwnd = ncread(filename,'uwnd');
uwnd = permute(mean(reshape(uwnd(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
ua = double(uwnd - mean(uwnd(:,:,42:71),3)); % remove climatology from 1981 to 2010
filename = 'G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\10mSI-MO\vwnd.10m.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
vwnd = ncread(filename,'vwnd');
vwnd = permute(mean(reshape(vwnd(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
va = double(vwnd - mean(vwnd(:,:,42:71),3)); % remove climatology from 1981 to 2010
filename = 'G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\miscSI-MO\prmsl.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
slp = ncread(filename,'prmsl');
slp = permute(mean(reshape(slp(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
slpa = double(slp - mean(slp(:,:,42:71),3)); % remove climatology from 1981 to 2010
uar = permute(ua,[3 1 2]);
var = permute(va,[3 1 2]);
slpar = permute(slpa,[3 1 2]);
dT = 1; cf = 1/8;
for i = 1:360
    for j = 1:181
    uad(:,i,j) = detrend(uar(:,i,j)); % linear detrend
    uadf(:,i,j) = lanczosfilter(uad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    vad(:,i,j) = detrend(var(:,i,j)); % linear detrend
    vadf(:,i,j) = lanczosfilter(vad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    slpad(:,i,j) = detrend(slpar(:,i,j)); % linear detrend
    slpadf(:,i,j) = lanczosfilter(slpad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
uadf = permute(uadf,[2 3 1]);
uadf(:,:,1:4) = []; uadf(:,181,:) = [];
vadf = permute(vadf,[2 3 1]);
vadf(:,:,1:4) = []; vadf(:,181,:) = [];
slpadf = permute(slpadf,[2 3 1]);
slpadf(:,:,1:4) = []; slpadf(:,181,:) = [];
%% wind stress
clear taux tauy
for s = 1:size(uadf,3)
    % lon*lat
    [taux(:,:,s) tauy(:,:,s)]= ra_windstr(uadf(:,:,s),vadf(:,:,s));
end
%  wind stress curl
clear curlz0
for s = 1:size(uadf,3)
    % curlz lat*lon
    curlz0(:,:,s) = ra_windstrcurl(latData,lonData,uadf(:,:,s)',vadf(:,:,s)',0);
end
curlz = permute(curlz0,[2 1 3]);
%% load sst
filename = 'G:\data\hadley\HadISST_sst.nc'; % 1870-2021
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
sst = ncread(filename,'sst');
sst(find(sst == -1000)) = nan;
lonDataH = ncread(filename,'longitude'); % -179.5-179.5
latDataH = ncread(filename,'latitude');  % 89.5- -89.5 
sst = permute(mean(reshape(sst(:,:,841:1812),[360 180 12 972/12]),3),[1 2 4 3]); % yearly 1940-2020
ssta = double(sst - mean(sst(:,:,42:71),3)); % remove climatology from 1981 to 2010
% 转换为相同lonData latData
ssta1 = fliplr(ssta);
ssta2 = cat(1,ssta1(182:end,:,:,:),ssta1(1:181,:,:,:));
sstar = permute(ssta2,[3 1 2]);
dT = 1; cf = 1/8;
clear sstad sstadf
for i = 1:360
    for j = 1:180
    sstad(:,i,j) = detrend(sstar(:,i,j)); % linear detrend
    sstadf(:,i,j) = lanczosfilter(sstad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
sstadf = permute(sstadf,[2 3 1]);
sstadf(:,:,1:4) = []; sstadf(:,:,end-3:end) = [];
%% load surface heat data
[lhtfla] = LoadNOAA('G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\sfcFlxSI-MO\lhtfl.mon.mean.nc','lhtfl');
[shtfla] = LoadNOAA('G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\sfcFlxSI-MO\shtfl.mon.mean.nc','shtfl');
[dlwrfa] = LoadNOAA('G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\sfcFlxSI-MO\dlwrf.sfc.mon.mean.nc','dlwrf');
[dswrfa] = LoadNOAA('G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\sfcFlxSI-MO\dswrf.sfc.mon.mean.nc','dswrf');
[ulwrfa] = LoadNOAA('G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\sfcFlxSI-MO\ulwrf.sfc.mon.mean.nc','ulwrf');
[uswrfa] = LoadNOAA('G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\sfcFlxSI-MO\uswrf.sfc.mon.mean.nc','uswrf');
% surface net heat flux (positive downward)
nhflxa = (dswrfa-uswrfa)-lhtfla-shtfla+dlwrfa-ulwrfa;
clear nhflxad nhflxadf
dT = 1; cf = 1/8;
for i = 1:360
    for j = 1:181
        nhflxad(:,i,j) = detrend(nhflxa(:,i,j)); % linear detrend
        nhflxadf(:,i,j) = lanczosfilter(nhflxad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
nhflxadf = permute(nhflxadf,[2 3 1]);
nhflxadf(:,:,1:4) = []; nhflxadf(:,181,:) = [];
%% load z200
filename = 'G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\hgt.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
hgt = ncread(filename,'hgt');
level = ncread(filename,'level');
z200raw = permute(hgt(:,:,19,:),[1 2 4 3]);
z200 = permute(mean(reshape(z200raw(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
z200za = z200-mean(z200,1); % zonal deviation
z200zaa = double(z200za - mean(z200za(:,:,42:71),3)); % remove climatology from 1981 to 2010
z200ar = permute(z200zaa,[3 1 2]);
dT = 1; cf = 1/8;
for i = 1:360
    for j = 1:181
    z200ad(:,i,j) = detrend(z200ar(:,i,j)); % linear detrend
    z200adf(:,i,j) = lanczosfilter(z200ad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
z200adf = permute(z200adf,[2 3 1]);
z200adf(:,:,1:4) = []; z200adf(:,181,:) = [];

%% EOF
lats = 1:90;
dweit = depthData(2:27)-depthData(1:26); % depth weight
Tsub_long = permute(nansum(Tempadlf_long(:,lats,1,:)+Tempadlf_long(:,lats,2:27,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m, 1940-2020
Tsub = permute(nansum(Tempadlf(:,lats,1,:)+Tempadlf(:,lats,2:27,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m
[eof_maps,pc,expvar]=eofnan(Tsub);  
EOF1z = eof_maps(:,:,1)*std(pc(1,:)); 
PC1z = pc(1,:)/std(pc(1,:)); 
%%
close all;
ftsz = 12; ticks = 0.15;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(EOF1z,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_45_60S\Tem0-700m_EOFmap_0_90S_zscore.png'],'-dpng','-r300')
%%  time series PC1 & TPI
pc1z = zscore(pc(1,:));
num_pc1p = find(pc1z>0.5);
num_pc1n = find(pc1z<-0.5);
close all;
Fig = figure('position',[2700 100 800 400])
plot(zscore(pc(1,:)/10),'-','color',[0,46,166]/255,'LineWidth',1.5)
hold on
plot(TPIfz,'color',[252,198,48]/255,'LineWidth',1.5)
plot(zeros(1,81),'k','LineWidth',1)
plot(zeros(1,81)-0.5,'k--','LineWidth',1)
plot(zeros(1,81)+0.5,'k--','LineWidth',1)
rr = corrcoef(pc1z,TPIfz);
text(5,-2.2,['Corr = ',num2str(roundn(rr(1,2),-2)),'**'],'fontsize',14)
scatter(num_pc1p,pc1z(num_pc1p),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
scatter(num_pc1n,pc1z(num_pc1n),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
set(gca,'XLim',[1,73]);
set(gca,'XTick',[7:10:73]);
set(gca,'XTickLabel',[1950:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
legend('PC1','TPI','Location','north','Orientation','horizontal')
legend('boxoff');
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\TimeSeries.png'],'-dpng','-r300')
%% PC1 time series  with CESM PAC pacemaker
pc1z = zscore(pc(1,:));
num_pc1p = find(pc1z>0.5);
num_pc1n = find(pc1z<-0.5);
close all;
Fig = figure('position',[2700 100 800 400])
maxval = max(PC1pac,[],1);
minval = min(PC1pac,[],1);
xval = [1:86];
fill([xval,fliplr(xval)],[minval fliplr(maxval)],[.8 .8 .8]);
hold on
plot(zeros(1,93),'k','LineWidth',1)
for s = 1:10
    plot([1:86],PC1pac(s,:),'LineWidth',0.5)
end
p1 = plot([21:93],PC1z,'-','color','k','LineWidth',3)
p2 = plot([1:86],PC1ensm(1,:),'color','r','LineWidth',3)
set(gca,'XLim',[1,93]);
set(gca,'XTick',[7:10:93]);
set(gca,'XTickLabel',[1930:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
legend([p1,p2],'IAP','ensmean','Location','south','Orientation','horizontal')
legend('boxoff');
% print(Fig,['G:\figures\CESM\Yearly\PAC-pacemaker\20230118_IPO_SouthernOcean\ensmean\PC1.png'],'-dpng','-r300')
%% IPO time series  with CESM PAC pacemaker
% pc1z = zscore(pc(1,:));
close all;
Fig = figure('position',[2700 100 800 400])
maxval = max(TPIpac,[],2);
minval = min(TPIpac,[],2);
xval = [1:86];
fill([xval,fliplr(xval)],[minval' fliplr(maxval')],[.8 .8 .8]);
hold on
plot(zeros(1,93),'k','LineWidth',1)
for s = 1:10
    plot([1:86],TPIpac(:,s),'LineWidth',0.5)
end
p1 = plot([21:93],TPIfz,'-','color','k','LineWidth',3)
p2 = plot([1:86],TPIensm,'color','r','LineWidth',3)
set(gca,'XLim',[1,93]);
set(gca,'XTick',[7:10:93]);
set(gca,'XTickLabel',[1930:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
legend([p1,p2],'IAP','ensmean','Location','north','Orientation','horizontal')
legend('boxoff');
print(Fig,['G:\figures\CESM\Yearly\PAC-pacemaker\20230118_IPO_SouthernOcean\ensmean\IPO-TPI.png'],'-dpng','-r300')
%% PC1 .reg. zonal mean 150E-60W
lonm1 = permute(nanmean(Tempadlf(150:300,lats,1:27,:),1),[2 3 4 1]);
% lonm1m = permute(nanmean(Temp(150:300,lats,1:27,:),1),[2 3 4 1]); % climotology
% corr with PC1
var = permute(lonm1,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = reg1_ttest(pc(1,:),var(:,i,j),0.05,1); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:27),par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% [ac ah] = contour(latData(lats),-depthData(1:27),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.PC1.zonal mean temperature (South Pac)')
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\Tem0-700m_EOFmap_0_90S_Pac.png'],'-dpng','-r300')
%% PC1 .reg. zonal mean 30W-150E 
Temp_r = cat(1,Tempadlf(300:360,:,:,:),Tempadlf(1:299,:,:,:));
lonm1 = permute(nanmean(Temp_r(1:210,lats,1:27,:),1),[2 3 4 1]);
% corr with PC1
var = permute(lonm1,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(pc(1,:),var(:,i,j),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:27),par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% [ac ah] = contour(latData,-depthData(1:27),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.PC1.zonal mean temperature (South IO+Atl)')
print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\Tem0-700m_EOFmap_0_90S_IO+Atl.png'],'-dpng','-r300')
%% PC1 .reg. zonal mean 0-360
lonm1 = permute(nanmean(Tempadlf(:,lats,1:27,:),1),[2 3 4 1]);
% lonm1m = permute(nanmean(Temp(150:300,lats,1:27,:),1),[2 3 4 1]); % climotology
% corr with PC1
var = permute(lonm1,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(pc(1,:),var(:,i,j),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
%
% close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:27),par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% [ac ah] = contour(latData(lats),-depthData(1:27),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.PC1.zonal mean temperature (Southern Ocean)')
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\Tem0-700m_EOFmap_0_90S_zm_0_360.png'],'-dpng','-r300')
%% PC1 .reg. meridional mean 
lats = [30:45]; names = '45S-60S'; index = DI;
latm = latmean(Tempadlf,lats,latData);
% latmc = latmean(Temp,lats,latData); % climotology
% corr with 
var = permute(latm,[3 1 2]);
clear par11 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i,j),0.1,1);
    end
end
max(par11,[],'all')
min(par11,[],'all')
%
close all;
ticks = 0.1;
par11_r = cat(1,par11(150:300,:),par11(301:360,:),par11(1:149,:));
lonData_r = [lonData(150:300);lonData(301:360);lonData(1:149)+360];
h05_r = cat(1,h05(150:300,:),h05(301:360,:),h05(1:149,:));

Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,par11_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([300 300],[0 -700],'linewidth',1.5,'color','k')
[xx,yy] = find(h05_r==1);
plot(lonData_r(xx),-depthData(yy),'.','color',[.5 .5 .5])
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','K','Units','normalized','position',[4 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[150,360+149]);
set(gca,'XTick',[150:60:360+149]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
% title(['PC1.reg.meridional mean temperature (',names,')'])
print(Fig,['G:\figures\IAP\Yearly\20230427_IPO_SouthernOcean_45_60S\DI.reg.Tem0-700m.png'],'-dpng','-r300')
%% time series South Pac, IO+Atl
lats = [30:45];
Tsub_r = cat(1,Tsub(150:300,:,:),Tsub(301:360,:,:),Tsub(1:149,:,:));
lonData_r = [lonData(150:300);lonData(301:360);lonData(1:149)+360];
[spacz spac] = areamean(Tsub_r,1:151,lats,latData); 
[siaz sia] = areamean(Tsub_r,152:360,lats,latData); 
DI = zscore(spac-sia);
% long 1940-2020
Tsub_longr = cat(1,Tsub_long(150:300,:,:),Tsub_long(301:360,:,:),Tsub_long(1:149,:,:));
[spacz_long spac_long] = areamean(Tsub_longr,1:151,lats,latData); 
[siaz_long sia_long] = areamean(Tsub_longr,152:360,lats,latData); 
DI_long = zscore(spac_long-sia_long);

%%
t1 = DI_long; t2 = PC1z; str1 = 'DI'; str2 = 'PC1'; pngname = 'DI_PC1';
% t1 = spacz_long; t2 = siaz_long; str1 = 'Pac'; str2 = 'IO+Atl'; pngname = 'Pac_IO+Atl';

[r,p,n_eff2] = corr_eff(t1(5:end-4),t2,0.05)  %  90% confidence
% [r,p,n_eff2] = corr_eff(t1(5:end-4),t2(5:end-4),0.05)  %  90% confidence
close all;
Fig = figure('position',[700 100 800 400]);
plot([5:77],t1(5:end-4),'-','color','r','LineWidth',1.5)
hold on
plot([5:77],t2(1:end),'color','b','LineWidth',1.5)
% plot([5:77],t2(5:end-4),'color','b','LineWidth',1.5)

plot([1:4],t1(1:4),'--','color','r','LineWidth',1.5)
plot([78:81],t1(end-3:end),'--','color','r','LineWidth',1.5)

% plot([1:4],t2(1:4),'--','color','b','LineWidth',1.5)
% plot([78:81],t2(end-3:end),'--','color','b','LineWidth',1.5)
plot(zeros(1,81),'k','LineWidth',1)
plot(zeros(1,81)-0.5,'k--','LineWidth',1)
plot(zeros(1,81)+0.5,'k--','LineWidth',1)
text(5,-2.2,['Corr = ',num2str(roundn(r,-2)),'**'],'fontsize',14)
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
legend(str1,str2,'Location','north','Orientation','horizontal')
legend('boxoff');
print(Fig,['G:\figures\IAP\Yearly\20230427_IPO_SouthernOcean_45_60S\TimeSeries_',pngname,'.png'],'-dpng','-r300')
%% PC1/DI .reg. taux tauy slp curlz
varu = permute(taux,[3 1 2]);
varv = permute(tauy,[3 1 2]);
varslp = permute(slpadf,[3 1 2]);
varcurlz = permute(curlz,[3 1 2]);  
index = DI(1:72); pngname = 'DI';
clear par1 h1 par2 h2 par3 h3 par4 h4
for i = 1:size(varu,2);
    for j = 1:size(varu,3);
        [par1(i,j),h1(i,j),t] = reg1_ttest(index,varu(:,i,j),0.05,0);
        [par2(i,j),h2(i,j),t] = reg1_ttest(index,varv(:,i,j),0.05,0);
        [par3(i,j),h3(i,j),t] = reg1_ttest(index,varslp(:,i,j),0.05,0);
        [par4(i,j),h4(i,j),t] = reg1_ttest(index,varcurlz(:,i,j),0.05,0);
    end
end
par1(:,1:10) = nan; par1(:,170:180)=nan;
par2(:,1:10) = nan; par2(:,170:180)=nan;
%% slp + tau
close all;
ftsz = 12; ticks1 = 60; ticks2 = 0.08;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par3(:,1:180),[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(150,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'Pa',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
d = 6; dd = 0.03*10^-3;
[x,y] = meshgrid(lonData,latData);
umap = par1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par2(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1.5,'MaxHeadSize',5) 
% text(1.2,0.1,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
%     'margin',5,'position',[-3.073,1.38,-1])
% m_quiver(10,-60,0.5./dd*cosd(87),0,0,'color','r','linewidth',1,'MaxHeadSize',5)
print(Fig,['G:\figures\IAP\Yearly\20230220_IPO_SouthernOcean_46_61S\',pngname,'.reg.windstress&slp.png'],'-dpng','-r300')
%% curlz & wind
% close all;
ftsz = 12; ticks1 = 1; ticks2 = 0.08;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par4*10^9,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
d = 6; dd = 0.03*10^-3;
[x,y] = meshgrid(lonData,latData);
umap = par1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par2(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1.5,'MaxHeadSize',5) 
% text(1.2,0.1,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
%     'margin',5,'position',[-3.073,1.38,-1])
% m_quiver(10,-60,0.5./dd*cosd(87),0,0,'color','r','linewidth',1,'MaxHeadSize',5)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(150,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'10^-^9 Pa m^-^1',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['G:\figures\IAP\Yearly\20230220_IPO_SouthernOcean_46_61S\',pngname,'.reg.windstress&curlz.png'],'-dpng','-r300')
%% taux
ftsz = 12; ticks1 = 6*10^-4;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par1,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(150,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'10^-^4 Pa',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1]*10^4)
print(Fig,['G:\figures\IAP\Yearly\20230220_IPO_SouthernOcean_46_61S\',pngname,'.reg.Taux.png'],'-dpng','-r300')
%% tauy
T17 = permute(nansum(Tempadlf(:,1:90,12:27,:).*(permute(dweit(11:26),[3 2 1])),3)/600,[1 2 4 3]); % 100-700m
varv = permute(tauy,[3 1 2]);
index = areamean(T17,30:90,30:50,latData);
clear  par2 h2 
for i = 1:size(varv,2);
    for j = 1:size(varv,3);
        [par2(i,j),h2(i,j),t] = reg1_ttest(index(1:72),varv(:,i,j),0.05,0);
    end
end

ftsz = 12; ticks1 = 4*10^-4;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par2,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'10^-^4 Pa',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1]*10^4)
print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\30_90E.100-700mT.reg.Tauy.png'],'-dpng','-r300')
%% PC1 reg meridional mean (v wind stress) 45-60S 
lats = [30:45]; names = '45S-60S'; index = PC1z(1:72);
latm = latmean(tauy,lats,latData);
% latmc = latmean(Temp,lats,latData); % climotology
% corr with TPI
var = latm';
clear par11 h05
for i = 1:size(var,2);
    [par11(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i),0.05,1);
end
max(par11,[],'all')
min(par11,[],'all')
%%
close all;
ticks = 0.1;
par11_r = cat(1,par11(150:300,:),par11(301:360,:),par11(1:149,:));
lonData_r = [lonData(150:300);lonData(301:360);lonData(1:149)+360];
Fig = figure('position',[100 100 800 400]);
plot(lonData_r,par11_r,'g','linewidth',1.5)
hold on
plot(lonData_r,zeros(360),'k','linewidth',1.5)
set(gca,'XLim',[150,360+149+1]);
set(gca,'XTick',[150:60:360+149+1]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-1.8,1.8]*10^-4);
set(gca,'YTick',[-1.8:0.6:1.8]*10^-4,'FontSize',14);
ylabel('Meridional wind stress (Pa)')
print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\PC1.reg.meridionalmean.Vwindstress.50_60S.png'],'-dpng','-r300')
%% PC1/DI reg zonal mean (SST)
reg_zonalmean(sstadf,PC1z,[-0.1,0.15],[-0.1:0.05:0.15],'SST (K)','PC1.reg.zonalmean_SST',latData)
% reg_zonalmean(sstadf,DI,[-0.1,0.15],[-0.1:0.05:0.15],'SST (K)','DI.reg.zonalmean_SST',latData)
%% PC1/DI reg zonal mean (net heat flux)
reg_zonalmean(nhflxadf,PC1z(1:72),[-1.5,1.5],[-1.5:0.5:1.5],'Net Heat Flux (W m^-^2)','PC1.reg.zonalmean_NetHeatFlux',latData)
reg_zonalmean(nhflxadf,DI(1:72),[-1.5,1.5],[-1.5:0.5:1.5],'Net Heat Flux (W m^-^2)','DI.reg.zonalmean_NetHeatFlux',latData)
%% PC1/DI reg zonal mean (zonal wind)
reg_zonalmean(taux,PC1z(1:72),[-4.5,4]*10^-4,[-4:1:4]*10^-4,'Zonal Wind Stress (Pa)','PC1.reg.zonalmean_Taux',latData)
reg_zonalmean(taux,DI(1:72),[-4.5,4]*10^-4,[-4:1:4]*10^-4,'Zonal Wind Stress (Pa)','DI.reg.zonalmean_Taux',latData)
%% PC1/DI reg zonal mean (curlz)
reg_zonalmean(curlz,PC1z(1:72),[-5,10]*10^-10,[-5:5:10]*10^-10,'Wind Stress Curl (Pa/m)','PC1.reg.zonalmean_curlz',latData)
reg_zonalmean(curlz,DI(1:72),[-5,10]*10^-10,[-5:5:10]*10^-10,'Wind Stress Curl (Pa/m)','DI.reg.zonalmean_curlz',latData)
%% PC1/DI reg zonal mean (slp)
reg_zonalmean(slpadf,PC1z(1:72),[-50,50],[-50:10:50],'Sea Level Pressure (Pa)','PC1.reg.zonalmean_SLP',latData)
reg_zonalmean(slpadf,DI(1:72),[-50,50],[-50:10:50],'Sea Level Pressure (Pa)','DI.reg.zonalmean_SLP',latData)

%% u v wind
close all;
ticks = 0.35; ftsz = 12;
Fig = figure('position',[100 100 800 400]);
% m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
hold on
[c h] = m_contourf(lonData,latData,par1(:,1:180)',[-ticks*5:ticks/5:ticks*5],'linestyle','none');
caxis([-ticks,ticks])
hold on
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks],'position',[0.78 0.09 0.03 0.85],'fontsize',12);
set(ch.Label,'String','m s^-^1 STD ^-^1','Units','normalized','position',[4 0.5 0],'Rotation',-90,'fontsize',ftsz);
m_coast('linewidth',1,'color','k');
m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
d = 6; dd = 0.1;
[x,y] = meshgrid(lonData,latData);
umap = par1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par2(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\PC1_corr_UV.png'],'-dpng','-r300')
%% scatter plot
close all;
num_pc1p = find(pc1z>0);
num_pc1n = find(pc1z<0);
figure('position',[2700 100 500 500])
pc1p1 = find(pc1z >= 0.5 & pc1z < 0.9); pc1p2 = find(pc1z >= 0.9 & pc1z < 1.3); 
pc1p3 = find(pc1z >= 1.3 & pc1z < 1.7); pc1p4 = find(pc1z >= 1.7); 
numhb = [num_pc1n,num_pc1p];
ch = scatter(pc1z(numhb),TPIfz(numhb),[],pc1z(numhb),'filled')
% hold on
% scatter(AMOfz(find(pc1z>0 & pc1z<0.5)),TPIfz(find(pc1z>0 & pc1z<0.5)),'ro')
% scatter(AMOfz(find(pc1z>-0.5 & pc1z<0)),TPIfz(find(pc1z>-0.5 & pc1z<0)),'bo')
set(gca,'XLim',[-2,2]);
set(gca,'XTick',[-2:1:2]);
set(gca,'YLim',[-2,2]);
set(gca,'YTick',[-2:1:2]);
text(2.15,-0.2,'PC1','FontSize',12)
text(-0.4,2.25,'IPO','FontSize',12)
%
new_fig_handle = convert_to_std_coordinate_system(gca,0);
caxis([-2,2]);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
cb = colorbar;
set(cb,'Ticks',[-2:0.4:2],'TickLabels',[-2:0.4:2],'Position',[0.88,0.1,0.03,0.8],'FontSize',12);
set(cb.Label,'String','PC1','Rotation',270,'position',[4 0 0],'FontSize',12)
% saveas(new_fig_handle,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\scatter_PC1_IPO.png'],'png')
%% PC1/DI.reg.SST
index = DI; pngname = 'DI';
var = permute(sstadf,[3 1 2]);
clear par5 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par5(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i,j),0.05,1);
    end
end
max(par5,[],'all')
min(par5,[],'all')
%
ftsz = 12; ticks1 = 0.2;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par5,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(150,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['G:\figures\IAP\Yearly\20230220_IPO_SouthernOcean_46_61S\',pngname,'.reg.sst.png'],'-dpng','-r300')
%% PC1/DI.reg. net heat flux
index = DI(1:72); pngname = 'DI';
var = permute(nhflxadf,[3 1 2]);
clear par6 h06
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par6(i,j),h06(i,j),t] = reg1_ttest(index,var(:,i,j),0.05,1);
    end
end
max(par6,[],'all')
min(par6,[],'all')
%
ftsz = 12; ticks1 = 5;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par6,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(150,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['G:\figures\IAP\Yearly\20230220_IPO_SouthernOcean_46_61S\',pngname,'.reg.nhflx.png'],'-dpng','-r300')
%% PC1 1std composite 
num_PC1p = find(PC1z > 1);
num_PC1n = find(PC1z < -1);
clear evs_PC1p evs_PC1n
[evs_PC1p] = event_index(num_PC1p,5) % consecutive 5 days
[evs_PC1n] = event_index(num_PC1n,5)
clear map1
for s = 1:length(evs_PC1p)
    map1(:,:,s) = mean(Tsub(:,:,evs_PC1p{s}),3);
end
close all;
[Fig] = com_map(mean(map1,3),0.3,lonData,latData(1:90))
% print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\PC1p_0.5std.com.Tem700m.png'],'-dpng','-r300')
clear map2
for s = 1:length(evs_PC1n)
    map2(:,:,s) = mean(Tsub(:,:,evs_PC1n{s}),3);
end
[Fig] = com_map(mean(map2,3),0.3,lonData,latData(1:90))
% print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\PC1n_0.5std.com.Tem700m.png'],'-dpng','-r300')
%% PC1 1 std composite (meridional mean Tem)
clear map1
for s = 1:length(evs_PC1n)
    map1(:,:,s) = mean(latm(:,:,evs_PC1n{s}),3);
end
close all;
ticks = 0.2;
map = mean(map1,3);
map_r = cat(1,map(150:300,:),map(301:360,:),map(1:149,:));
lonData_r = [lonData(150:300);lonData(301:360);lonData(1:149)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,map_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([300 300],[0 -700],'linewidth',1.5,'color','k')
line([150 509],[-700 -700],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','K','Units','normalized','position',[4 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[150,360+149]);
set(gca,'XTick',[150:60:360+149]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\PC1n_1std.com.meridionalmeanT.png'],'-dpng','-r300')
%% PC1 corr other time series (scatter)
lats = [29:44];
var = nhflxadf; varname = 'NHFLX'; titlename = 'Net Heat Flux';
var_r = cat(1,var(150:300,:,:),var(301:360,:,:),var(1:149,:,:));
lonData_r = [lonData(150:300);lonData(301:360);lonData(1:149)+360];
[tz1 t1] = areamean(var_r,1:151,lats,latData); % Pac
[tz2 t2] = areamean(var_r,152:360,lats,latData); % Atl & IO
%
% tx1 = spacz(1:72); tx2 = tz1; pocname = 'Pac';
% tx1 = siaz(1:72); tx2 = tz2; pocname = 'IOAtl';
% tx1 = PC1z(1:72); tx2 = zscore(t1-t2); pocname = 'PC1';
% tx1 = DI(1:72); tx2 = zscore(t1-t2); pocname = 'DI';
tx1 = TPIfz(1:72); tx2 = zscore(t1-t2); pocname = 'IPO';
%
xlimv = [-3,3];
ylimv = [-4,3];
[r,p,n_eff2] = corr_eff(tx1(1:72),tx2,0.05)  %  90% confidence  
close all;
Fig = figure('position',[700 100 400 400]);
scatter(tx1,tx2,'*','Color',[0,46,166]/255,'LineWidth',2);
hold on
yb = polyfit(tx1,tx2,1);
plot(tx1,polyval(yb,tx1),'k','linewidth',3)
text(-2.5,-2.8,['Corr = ',num2str(roundn(r,-2))],'fontsize',12)
text(-2.5,-3.4,['Slope = ',num2str(roundn(yb(1),-2))],'fontsize',12)
set(gca,'XLim',xlimv);
set(gca,'FontSize',12);
set(gca,'YLim',ylimv);
xlabel('Temperature','FontSize',12),ylabel(titlename,'FontSize',12);
print(Fig,['G:\figures\IAP\Yearly\20230220_IPO_SouthernOcean_46_61S\',pocname,'.regwith.',varname,'.scatter.png'],'-dpng','-r300')

%% MOHT calculate
% ro*cp*OHC*Vy
% Vy = My/ro = -taux/f/ro 
datadir='G:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
f = sw_f(latData);
clear tauxraw tauyraw
for s = 1:size(vwnd,3)
    % lon*lat
    [tauxraw(:,:,s) tauyraw(:,:,s)]= ra_windstr(uwnd(:,1:180,s),vwnd(:,1:180,s));
end
My = -tauxraw./f'; % Ekman 质量输运
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
Vy = My/ro;
Vy(:,85:96,:) = nan; % 赤道不适用ekman输送
londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
levdist(1) = 1; levdist(2:27) = dweit; % z distance
clear OHC OHCT MOHT
OHC = permute(ro*cp*sum((Temp(:,1:180,1:27,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3])/700;
OHCT = OHC(:,:,1:76).*Vy;
OHCTw = OHCT.*londist';
OHCTw_r = cat(1,OHCTw(150:300,:,:),OHCTw(301:360,:,:),OHCTw(1:149,:,:));
MOHT1 = permute(nansum(OHCTw_r(1:151,:,:),1),[2 3 1]); % Pac unit:W
MOHT2 = permute(nansum(OHCTw_r(152:end,:,:),1),[2 3 1]); % IO+Atl unit:W
%%
close all;
MOHT = permute(nansum(OHCTw_r,1),[2 3 1]); % unit:W
Fig = figure('position',[100 100 600 300])
map = mean(MOHT,2)/10^15;
% map(85:96,:) = nan;
plot(latData,map,'linewidth',1.5)
set(gca,'XLim',[-90,90],'YLim',[-2,2],'XGrid','on','YGrid','on');
set(gca,'XTick',[-80:20:80],'YTick',[-2:0.5:2]);
ylabel('MOHT (PW)');xlabel('Latitude')
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_45_60S\MOHT_climotology.png'],'-dpng','-r300')
%% MOHT anomaly,detrend, and 8-yr filter
MOHT1a = MOHT1-mean(MOHT1,2);
MOHT2a = MOHT2-mean(MOHT2,2);
dT = 1; cf = 1/8;
for i = 1:180
    MOHT1ad(:,i) = detrend(MOHT1a(i,:)'); % linear detrend
    MOHT1adf(:,i) = lanczosfilter(MOHT1ad(:,i),dT,cf,[],'low'); % 8 year filtered
    MOHT2ad(:,i) = detrend(MOHT2a(i,:)'); % linear detrend
    MOHT2adf(:,i) = lanczosfilter(MOHT2ad(:,i),dT,cf,[],'low'); % 8 year filtered
end
MOHT1adf(1:4,:) = []; MOHT2adf(1:4,:) = [];
%%
index = TPIfz(1:72); pngname = 'IPO';
clear par7 h07 par8 h08 
for i = 1:size(MOHT1adf,2);
    [par7(i),h07(i),t] = reg1_ttest(index,MOHT1adf(:,i),0.05,1);
    [par8(i),h08(i),t] = reg1_ttest(index,MOHT2adf(:,i),0.05,1);
end
%
close all;
Fig = figure('position',[100 100 800 400]);
plot(latData,par7/10^15,'r','linewidth',1.5)
hold on
plot(latData,par8/10^15,'b','linewidth',1.5)
plot(latData,zeros(1,length(latData)),'k','linewidth',1.5)
set(gca,'XLim',[-70,-30]);
set(gca,'XTick',[-70:10:-30]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',[-.01,.01],'YTick',[-.01:.005:.01]);
set(gca,'YTick',yticks,'FontSize',14);
ylabel('MOHT (PW)')
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_45_60S\',pngname,'.reg.zonalmean_MOHT.png'],'-dpng','-r300')

%% IPO+
num_ipop = find(TPIfz>0.5);
num_ipon = find(TPIfz<-0.5);
var = permute(Tsub(:,1:55,:),[3 1 2]);
index = TPIfz;
clear par11 h01 par12 h02
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h01(i,j),t] = reg1_ttest(index,var(:,i,j),0.05,1);
    end
end
%
lats = 1:55;
close all;
ftsz = 12; ticks = 0.15;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par11,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(150,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['G:\figures\IAP\Yearly\20230220_IPO_SouthernOcean_46_61S\IPO_reg_Tem0-700m.png'],'-dpng','-r300')
%% zonal mean Pac 150E-60W
Temcom1z = permute(nanmean(Temcom1(150:300,lats,:),1),[2 3 4 1]);
close all;
ticks = 0.1;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData,Temcom1z',[-1:0.02:1],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' IPO+ zonal mean Tem (South Pac)'])
print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\IPOp_zonalmean_Tem_150E_60W.png'],'-dpng','-r300')
%% zonal mean IO+Atl 60W-150E
lats = 1:45;
Temcom1_r = cat(1,Temcom1(300:360,:,:,:),Temcom1(1:299,:,:,:));
Temcom1z = permute(nanmean(Temcom1_r(1:210,lats,1:27,:),1),[2 3 4 1]);
% Temcom1z = permute(nanmean(Temcom1(150:300,lats,:),1),[2 3 4 1]);
close all;
ticks = 0.1;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:27),Temcom1z',[-1:0.02:1],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' IPO+ zonal mean Tem (South IO+Atl)'])
print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\IPOp_zonalmean_Tem_60W_150E.png'],'-dpng','-r300')
%% meridional mean
lats = [29:44]; names = '46S-61S'; index = TPIfz;
latm = latmean(Tempadlf,lats,latData);
% latmc = latmean(Temp,lats,latData); % climotology
% corr with TPI
var = permute(latm,[3 1 2]);
clear par11 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i,j),0.05,1);
    end
end
max(par11,[],'all')
min(par11,[],'all')
%%
close all;
ticks = 0.15;
par11_r = cat(1,par11(150:300,:),par11(301:360,:),par11(1:149,:));
lonData_r = [lonData(150:300);lonData(301:360);lonData(1:149)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,par11_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([300 300],[0 -700],'linewidth',1.5,'color','k')
line([150 509],[-700 -700],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
set(ch.Label,'String','K','Units','normalized','position',[4 0.5 0],'Rotation',-90,'fontsize',12);
set(gca,'XLim',[150,360+149]);
set(gca,'XTick',[150:60:360+149]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title(['Reg.IPO.meridional mean temperature (',names,')'])
print(Fig,['G:\figures\IAP\Yearly\20230109_IPO_SouthernOcean\IPO_reg_meridionalmean_Tem_45-60S.png'],'-dpng','-r300')
%% IPO reg taux tauy slp curlz 
varu = permute(taux,[3 1 2]);
varv = permute(tauy,[3 1 2]);
varslp = permute(slpadf,[3 1 2]);
varcurlz = permute(curlz,[3 1 2]);  
index = TPIfz(1:72);
clear par1 h1 par2 h2 par3 h3 par4 h4
for i = 1:size(varu,2);
    for j = 1:size(varu,3);
        [par1(i,j),h1(i,j),t] = reg1_ttest(index,varu(:,i,j),0.05,0);
        [par2(i,j),h2(i,j),t] = reg1_ttest(index,varv(:,i,j),0.05,0);
        [par3(i,j),h3(i,j),t] = reg1_ttest(index,varslp(:,i,j),0.05,0);
        [par4(i,j),h4(i,j),t] = reg1_ttest(index,varcurlz(:,i,j),0.05,0);
    end
end
par1(:,1:10) = nan; par1(:,170:180)=nan;
%% slp + tau
close all;
ftsz = 12; ticks1 = 60; ticks2 = 0.08;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par3(:,1:180),[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'Pa',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
d = 6; dd = 0.03*10^-3;
[x,y] = meshgrid(lonData,latData);
umap = par1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par2(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1.5,'MaxHeadSize',5) 
print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\IPO.reg.windstress&slp.png'],'-dpng','-r300')
%% curlz & tau
% close all;
ftsz = 12; ticks1 = 1; ticks2 = 0.08;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par4*10^9,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'10^-^9 Pa m^-^1',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
d = 6; dd = 0.03*10^-3;
[x,y] = meshgrid(lonData,latData);
umap = par1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par2(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1.5,'MaxHeadSize',5) 
print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\IPO.reg.windstress&curlz.png'],'-dpng','-r300')
%% taux
ftsz = 12; ticks1 = 6*10^-4;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par1,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'10^-^4 Pa',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1]*10^4)
print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\IPO.reg.Taux.png'],'-dpng','-r300')
%% IPO.reg.SST
index = TPIfz;
var = permute(sstadf,[3 1 2]);
clear par5 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par5(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i,j),0.05,1);
    end
end
max(par5,[],'all')
min(par5,[],'all')
%
ftsz = 12; ticks1 = 0.2;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par5,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['G:\figures\IAP\Yearly\20230109_IPO_SouthernOcean\IPO.reg.sst.png'],'-dpng','-r300')
%% IPO.reg. net heat flux
index = TPIfz(1:72);
var = permute(nhflxadf,[3 1 2]);
clear par6 h06
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par6(i,j),h06(i,j),t] = reg1_ttest(index,var(:,i,j),0.05,1);
    end
end
max(par6,[],'all')
min(par6,[],'all')
%%
close all;
ftsz = 12; ticks1 = 5;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par6,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\IPO.reg.nhflx.png'],'-dpng','-r300')
%% IPO reg zonal mean (SST)
reg_zonalmean(sstadf,TPIfz,[-0.1,0.15],[-0.1:0.05:0.15],'SST (K)','IPO.reg.zonalmean_SST',latData)
%% IPO reg zonal mean (net heat flux)
reg_zonalmean(nhflxadf,TPIfz(1:72),[-1.5,1.5],[-1.5:0.5:1.5],'Net Heat Flux (W m^-^2)','IPO.reg.zonalmean_NetHeatFlux',latData)
%% IPO reg zonal mean (zonal wind)
reg_zonalmean(taux,TPIfz(1:72),[-4,4]*10^-4,[-4:1:4]*10^-4,'Zonal Wind Stress (Pa)','IPO.reg.zonalmean_Taux',latData)
%% IPO reg zonal mean (curlz)
reg_zonalmean(curlz,TPIfz(1:72),[-0.5,0.9]*10^-9,[-0.8:0.4:0.8]*10^-9,'Wind Stress Curl (Pa/m)','IPO.reg.zonalmean_curlz',latData)
%% IPO reg zonal mean (slp)
reg_zonalmean(slpadf,TPIfz(1:72),[-50,50],[-50:10:50],'Sea Level Pressure (Pa)','IPO.reg.zonalmean_SLP',latData)

%% (teleconnection) IPO reg Z200 zonal deviation
varz = permute(z200adf,[3 1 2]);
index = TPIfz(1:72); pngname = 'IPO';
clear par1 h1 par2 h2 par3 h3 par4 h4
for i = 1:size(varz,2);
    for j = 1:size(varz,3);
        [par1(i,j),h1(i,j),t] = reg1_ttest(index,varz(:,i,j),0.05,0);
    end
end
par1(:,1:10) = nan; par1(:,170:180)=nan;
%%
close all;
ticks1 = 10; ftsz = 12;
Fig = figure('position',[100 100 650 300])
contourVARra(par1,[-12:1:12],ticks1,0,360,-90,90,lonData,latData)
set(gca,'position',[0.08,0.11,0.775,0.815]);
set_colorbar([0.87 0.095 0.026 0.832],'Z200 (gpm)',3.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_45_60S\IPO.reg.Z200za.png'],'-dpng','-r300')

%%
close all;
ftsz = 12; ticks1 = 4*10^-3;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par5,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'10^-^4 Pa',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1]*10^4)















function [varar] = LoadNOAA(filename,varname)
% filename is where the file locates
% varname is name of the variable in raw file
% filename = 'G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\miscSI-MO\prmsl.mon.mean.nc'; % 1836-2015
    ncid=netcdf.open(filename,'NOWRITE');
    ncdisp(filename);
    var = ncread(filename,varname);
    var = permute(mean(reshape(var(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
    vara = double(var - mean(var(:,:,:),3)); % remove climatology from 1981 to 2010
%42:71
    varar = permute(vara,[3 1 2]);
end
function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function contourfSPolar(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end
function contourVARra(var,val,ctr,lon1,lon2,lat1,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[lat1 lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
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
function [Fig] = diffdepth(Temcom,lonData,latData,lats,val,ticks)
    Fig = figure('position',[100 100 600 300]);
%     contourVARra(Temcom(:,:,nlevel),[-val:ticks/10:val],1,0,360,-90,90,lonData,latData);
    contourfSPolar(Temcom,[-val:ticks/10:val],ticks,lonData,latData(lats),12)
    caxis([-ticks,ticks]);
    ch = colorbar;
    set(ch,'Ticks',[-ticks:ticks/5:ticks],'position',[0.78 0.09 0.03 0.85]);
    % testdots(h05,'k',lonData,latData)
end
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end
function [res] = gradient_lat(lat,var)
% 此函数用于计算沿经向的梯度
% Input:
% var is lat*lon 
% lat is 单调递增，自南向北
% Output:
% res 即为梯度结果。
% degree to radian
rad=pi/180;
% computation of curl
[lt, ln]=size(var);
a=diff(lat);
aa=NaN*ones(length(a)-1,1);
for ii=1:length(a)-1
    if (a(ii) == a(ii+1))
        aa(ii)=a(ii);
    else
        error('Latitude difference is not consistance')
    end % endif
    dlat=mean(aa);
end % endfor
clear ii
deltay=dlat*111176;
res=NaN(lt, ln);
% Centeral difference method in x and y 
for ii=2:lt-1
    for jj=1:ln
        res(ii, jj)= (var(ii+1, jj)-var(ii-1, jj))/(2*deltay) ;
    end % endfor
end % endfor
clear ii jj
% Forward difference method in x and y 
for jj=1:ln
    res(1, jj)= (var(2, jj)-var(1, jj))/deltay ; 
end 
clear ii jj
% Backward difference method in x and y
for jj=1:ln
    res(lt, jj)= (var(lt, jj)-var(lt-1, jj))/deltay ;
end
clear ii jj
end
function reg_zonalmean(nhflx,index,ylimval,yticks,titlename,pngname,latData)
lats = [1:90];
lonmh1 = permute(nanmean(nhflx(150:300,lats,:),1),[2 3 1]);
nhflx_r = cat(1,nhflx(300:360,:,:),nhflx(1:299,:,:));
lonmh2 = permute(nanmean(nhflx_r(1:210,lats,:),1),[2 3 1]);
% corr with PC1
var1 = lonmh1';
var2 = lonmh2';
clear par11 h01 par12 h02
for i = 1:size(var1,2);
    [par11(i),h01(i),t] = reg1_ttest(index,var1(:,i),0.05,1); % filtered
    [par12(i),h02(i),t] = reg1_ttest(index,var2(:,i),0.05,1); % filtered
end
max(par11,[],'all')
min(par11,[],'all')
%
close all;
Fig = figure('position',[100 100 800 400]);
plot(latData(lats),par11,'r','linewidth',1.5)
hold on
plot(latData(lats),par12,'b','linewidth',1.5)
plot(latData(lats),zeros(1,length(lats)),'k','linewidth',1.5)
set(gca,'XLim',[-70,-30]);
set(gca,'XTick',[-70:10:-30]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',ylimval);
set(gca,'YTick',yticks,'FontSize',14);
ylabel(titlename)
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_45_60S\',pngname,'.png'],'-dpng','-r300')
end

function [eyr] = event_index(var,days)
% This function is used to find events which is defined by consecutive 
% days/years.
% Input var is all days/years that satisfied the condition. e.g. 1std
% Input days is the number of defined consecutive days/years
% Output eyr is the choosed events
enum = 1; k = 0;
for s = 2:length(var);
    diff = var(s)-var(s-1);
    if diff == 1 
        k = k + 1;
        if k >= days-1
            eyr{enum} = var(s-k:s);
        end
    else 
        if k >= days-1
            eyr{enum} = var(s-1-k:s-1);
            enum = enum + 1;
        end
        k = 0;
    end
end
end

function [Fig] = com_map(var,ticks1,lonData,latData)
ftsz = 12;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(var,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-60,'linewidth',2,'color','k'); 
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
end
function testdots(h,clor,lonData,latData)
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',4,'color',clor); % F-test dots
end
function h0 = trendtest(var,alp)
% note F95val need to search based F-test and sample number
% var is time*lon*lat
d1 = size(var,1);
d2 = size(var,2);
d3 = size(var,3);
h0 = zeros(d2,d3);
x = [1:d1]';
for i=1:d2;
    for j=1:d3;
        m3=var(:,i,j);
        par(i,j,:)=polyfit(x,m3,1); % regression parameters
        y1=polyval(permute(par(i,j,:),[3 2 1]),x);
        Q = sum((m3-y1).^2); % t test
        c = 1/sum((x-mean(x)).^2);
        T = par(i,j,1)/sqrt(c)/sqrt(Q/(d1-2));
        if abs(T) > tinv(1-alp/2,d1-2); % 双侧 t test
            h0(i,j) = 1;
        end
    end
end
end