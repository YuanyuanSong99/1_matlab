% Trend 0-700m Tem150
clc,clear,close all;
addpath(genpath('E:\1_matlab\help'));
filename1 = 'G:\CMIP6\historical\ptem\ensmean.thetao.1850-2014.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
Tem = ncread(filename1,'thetao');
Temp = Tem(:,:,:,111:end); % 1960-2014
clear Tem
%%
latData = ncread(filename1,'lat');
lonData = ncread(filename1,'lon');
depthData = ncread(filename1,'lev');
%%
Tempa = double(Temp - mean(Temp(:,:,:,22:51),4,'omitmissing')); % remove climatology from 1981-2010 
%% meridional mean
lats = 35:55; % 90S-35S
latstr = '35S_55S';
[Tempa_mm] = latmean(Tempa,lats,latData);
% ts is lon*depth*time
startyr = 1960;
endyr = 2014;
yrstr = '1960-2014';
var = permute(Tempa_mm(:,:,startyr-1959:endyr-1959),[3 1 2]);
x = [1:size(var,1)]';
clear trd
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.05); % t test trend
max(trd,[],'all')
min(trd,[],'all')
%% 150 longitude
close all;
ftsz = 12;
ticks = 0.25;
map = trd*10;
map_r = cat(1,map(150:360,:),map(1:150,:));
lonData_r = [lonData(150:360);lonData(1:150)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,map_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([300 300],[0 -2000],'linewidth',1.5,'color','k')
% line([150 509],[-700 -700],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('E:/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(flipud(bl_re4));
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','K decade^-^1','Units','normalized','position',[4.5 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[150,360+150]);
set(gca,'XTick',[150:60:360+150]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:400:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-400:0],'FontSize',14);
ylabel('Depth (m)')
% print(Fig,['E:\figures\IAP\Yearly\20240416_SO_warming\meridionalmean_',latstr,'_trend_',yrstr,'_150E.png'],'-dpng','-r300')
%% OHC
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthData(2:53)-depthData(1:52); % depth weight
OHCsub = cp*ro*permute(nansum(cat(3,Temp(:,:,1,:),Temp(:,:,2:53,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]); % J/m2

OHCw = OHCsub/365/24/60/60; % W/m2
%% OHC trend
startyr = 1960;
endyr = 2014;
var = permute(OHCw(:,:,startyr-1959:endyr-1959),[3 1 2]);
x = [1:size(var,1)]';
clear trd 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.05); % t test trend
h0(find(isnan(Temp(:,:,1,1)) == 1)) = nan;
max(trd,[],'all')
min(trd,[],'all')
%% SO plot
close all;
ftsz = 12; ticks = 1.5;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd,[-ticks*5:ticks/10:ticks*5],ticks,lonData,latData,12)
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-60,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots0(h0,[.4 .4 .4],lonData,latData);
print(Fig,['D:\figures\CMIP6\SOtrend\ensmean\20240423_SO_warming\OHC',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% original OHC time series ZJ
lats = [35:55]; % 55S-35S
levs = [1:53]; % 700 m
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
Temp = double(Temp);
Tempz = cat(3,Temp(:,:,1,:),Temp(:,:,2:53,:).*(permute(dweit,[3 2 1])));
VVz = cat(3,VV(:,:,1),VV(:,:,2:53).*(permute(dweit,[3 2 1])));
Tempzy = Tempz*111*1000; % 111 km / latitude
VVzy = VVz*111*1000;
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(150:360,:,:,:),Tempzyx(1:149,:,:,:)); 
VVzyx(find(isnan(Tempzyx(:,:,:,1)) == 1)) = nan;
Vzyxsub_r_0 = cat(1,VVzyx(150:360,:,:),VVzyx(1:149,:,:)); 
spac_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),2),3));
spacV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,lats,levs),1),2),3));
spac0V = spac_0./spacV;
sia_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),2),3));
siaV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,lats,levs),1),2),3));
sia0V = sia_0./siaV;
sdiff = spac0V-sia0V;
max(spac0V,[],'all')
%% OHC/m3  time series
close all;
Fig = figure('position',[700 100 800 400]);
p1 = bar([1:55],sdiff,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7],'BarWidth',1);
hold on
set(gca,'YLim',[-0,1],'YTick',[-0:0.2:1],'YColor','r');
ylabel('Difference  (W m^-^3)','FontSize',14);
yb = polyfit([1:55],sdiff,1);
plot([1:55],polyval(yb,[1:55]),'--','Color',[0.3,0.3,0.3],'linewidth',2)
text(55,0.2,[num2str(roundn(yb(1),-4)*10),' W/m^3 per decade'],'fontsize',14,'Color',[0.3,0.3,0.3])

yyaxis right
ymax = 4; ymin = 3; yint = 0.2;
p2 = plot(sia0V,'-','color',[.64 .08 .18],'LineWidth',3);
hold on
p3 = plot(spac0V,'-','color',[.04 .35 .56],'LineWidth',3);
plot(zeros(1,55),'k','LineWidth',1)
yb = polyfit([1:55],sia0V,1);
plot([1:55],polyval(yb,[1:55]),'--','Color',[.64 .08 .18],'linewidth',2)
text(55,3.5,[num2str(roundn(yb(1),-4)*10),' W/m^3 per decade'],'fontsize',14,'Color',[.64 .08 .18])
yb = polyfit([1:55],spac0V,1);
plot([1:55],polyval(yb,[1:55]),'--','Color',[.04 .35 .56],'linewidth',2)
text(55,3.8,[num2str(roundn(yb(1),-4)*10),' W/m^3 per decade'],'fontsize',14,'Color',[.04 .35 .56])
set(gca,'XLim',[1,55]);
set(gca,'XTick',[1:10:55]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel(' (W m^-^3)','FontSize',14);
legend([p3, p2, p1],'Pacific','Atlantic-Indian Ocean','Difference','Location','northwest')
legend('boxoff')
% print(Fig,['E:\figures\IAP\Yearly\20240416_SO_warming\TimeSeries_OHC_Wm3_',depthstr,'_35S_90S.png'],'-dpng','-r300')
%% OHC W  time series （anomaly）
sdiff = (spac_0-sia_0)/10^21;
bsln = 'baseline_1960-1970';
ln1 = (spac_0-mean(spac_0(1:11),1))/10^21; % baseline 1960-1970
ln2 = (sia_0-mean(sia_0(1:11),1))/10^21;
close all;
Fig = figure('position',[700 100 800 400]);
p3 = plot(ln1,'-','color',[.04 .35 .56],'LineWidth',3);
hold on
set(gca,'YLim',[-50,50],'YTick',[-50:10:20],'YColor','r');
ylabel('Pacific  (ZJ)','FontSize',14);
yb = polyfit([1:55],spac_0/10^15,1);
plot([1:55],polyval(yb,[1:55]),'--','Color',[.04 .35 .56],'linewidth',2)
text(55,221,[num2str(roundn(yb(1),-4)*10),' PW per decade'],'fontsize',14,'Color',[.04 .35 .56])

yyaxis right
ymax = 50; ymin = -50; yint = 10;
p2 = plot(ln2,'-','color',[.64 .08 .18],'LineWidth',3);
hold on
plot(zeros(1,55),'k','LineWidth',1)
yb = polyfit([1:55],sia_0/10^15,1);
plot([1:55],polyval(yb,[1:55]),'--','Color',[.64 .08 .18],'linewidth',2)
text(55,390,[num2str(roundn(yb(1),-4)*10),' PW per decade'],'fontsize',14,'Color',[.64 .08 .18])
set(gca,'XLim',[1,55]);
set(gca,'XTick',[1:10:55]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel('Atlantic-Indian Ocean (ZJ)','FontSize',14);
legend([p3, p2],'Pacific','Atlantic-Indian Ocean','Location','northwest')
legend('boxoff')
% print(Fig,['E:\figures\IAP\Yearly\20240416_SO_warming\TimeSeries_OHC_PW_',depthstr,'m_35S_90S_',bsln,'.png'],'-dpng','-r300')

%%  load u v ；taux tauy curlz
datadir='E:\data\ERA5_1x1\10m_wind\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE');
lonDataE = ncread(filetemp,'longitude'); % 0-360
lonEmap = lonDataE; lonEmap(361) = lonEmap(1);
latDataE = ncread(filetemp,'latitude'); % 90 - -90
%
k=length(filelist);
for s=1:k-2;  % 1940-2020
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    u10 = ncread(filename1,'u10');     uwnd(:,:,s) = mean(u10,3);   % yearly mean
    v10 = ncread(filename1,'v10');     vwnd(:,:,s) = mean(v10,3);   % yearly mean
    netcdf.close(ncid);  %关闭文件    
    [taux_raw(:,:,s) tauy_raw(:,:,s)]= ra_windstr(uwnd(:,:,s),vwnd(:,:,s)); % lon*lat
    curlz0(:,:,s) = ra_windstrcurl(latDataE,lonDataE,uwnd(:,:,s)',vwnd(:,:,s)',0); % lat*lon
end
curlz_raw = permute(curlz0,[2 1 3]);
%%  load slp
datadir='E:\data\ERA5_1x1\sea_level_pressure\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE');
k=length(filelist);
for s=1:k-2;  % 1940-2020
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    slp0 = ncread(filename1,'msl');     slp(:,:,s) = mean(slp0,3);   % yearly mean
    netcdf.close(ncid);  %关闭文件    
end
%% load sst and sic
filename = 'E:\data\hadley\HadISST_sst.nc'; % 1870-2021
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
sst = ncread(filename,'sst');
sst(find(sst == -1000)) = nan;
lonDataH = ncread(filename,'longitude'); % -179.5-179.5
latDataH = ncread(filename,'latitude');  % 89.5- -89.5 
sst = permute(mean(reshape(sst(:,:,841:1812),[360 180 12 972/12]),3),[1 2 4 3]); % yearly 1940-2020
% 转换为与ERA5相同lonData latData
clear sst1 sst_raw
sst_raw = cat(1,sst(181:end,:,:),sst(1:180,:,:));
sst_raw(:,181,:) = sst_raw(:,180,:);

filename = 'E:\data\hadley\HadISST_ice.nc'; % 1870-2021
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
sic = ncread(filename,'sic');
% sst(find(sst == -1000)) = nan;
sic = permute(mean(reshape(sic(:,:,841:1812),[360 180 12 972/12]),3),[1 2 4 3]); % yearly 1940-2020
% 转换为与ERA5相同lonData latData
sic_raw = cat(1,sic(181:end,:,:),sic(1:180,:,:));
sic_raw(:,181,:) = sic_raw(:,180,:);

%% load surface heat data
datadir='E:\data\ERA5_1x1\heat\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
k=length(filelist);
for s=1:k-2;  % 1940-2020
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    slhf = ncread(filename1,'slhf');     lhf(:,:,s) = mean(slhf,3);   % latent heat flux  upward  J m-2 
    ssr = ncread(filename1,'ssr');       swf(:,:,s) = mean(ssr,3);    % short wave flux   downward
    str = ncread(filename1,'str');       lwf(:,:,s) = mean(str,3);    % long wave flux    upward
    sshf = ncread(filename1,'sshf');     shf(:,:,s) = mean(sshf,3);   % sensible heat flux upward    
    netcdf.close(ncid);  %关闭文件    
end
%% climo rectangle
nhflx = (swf+lhf+lwf+shf)/3600; % positive downward W m-2
mapcli = mean(nhflx,3,'omitmissing');
max(mapcli,[],'all')
min(mapcli,[],'all')
close all;
ftsz = 12; ticks = 100;
Fig = figure('position',[10 50 850 400]);
ax = axes('Position',[0.05 0.07 0.85 0.87],'fontsize',ftsz,'box','on');
contourVARra(mapcli,[-ticks*5:ticks/10:ticks*5],ticks,0,360,-90,90,lonDataE,latDataE)
% title(['Climotological net heat flux (downward positive)'])
set_colorbar([0.9 0.065 0.02 0.875],'Pa m^-^1',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% print(Fig,['E:\figures\IAP\Yearly\climo_Curlz_ERA5.png'],'-dpng','-r300')
%% climo polar
close all;
% mapcli = nanmean(nhflx(:,:,1960-1939:2020-1939),3);
% ticks = 100;

% mapcli = nanmean(sst(:,:,1960-1939:2020-1939),3);
% ticks = 10;

% mapcli = nanmean(taux_raw(:,:,1960-1939:2020-1939),3);
% ticks = 0.08;

mapcli = nanmean(curlz_raw(:,:,1960-1939:2020-1939),3);
ticks = 1*10^-7;

Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(mapcli,[-ticks*5:ticks/10:ticks*5],ticks,lonDataE,latDataE,12)
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-60,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],'K',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['E:\figures\IAP\Yearly\20240416_SO_warming\CURLZ_climo_1960_2020_Polar_ERA5.png'],'-dpng','-r300')
%% trend
startyr = 1960;
endyr = 2020;
var = permute(sic_raw(:,:,startyr-1939:endyr-1939),[3 1 2]);
trd = trend_cal_3D(var);
h0 = trendtest(var,0.05); % t test trend
h0(find(isnan(Temp(:,:,1,1)) == 1)) = nan;
max(trd,[],'all'), min(trd,[],'all')
%%
% units = '10^-^9 Pa m^-^1'; name = 'CURLZ'; 
% map = trd*10*10^9; ticks = 3;

% units = 'K'; name = 'SST'; 
% map = trd*10; ticks = .2;

units = '%'; name = 'SIC'; 
map = trd*10*100; ticks = 6;
map(find(map == 0)) = nan;

% units = '10^-^3 Pa'; name = 'TAUX'; 
% map = trd*10*10^3; ticks = 3;

% units = 'W m^-^2'; name = 'sensitiveheat'; 
% map = trd*10; ticks = 4;

close all;
ftsz = 12; 
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(map,[-ticks*5:ticks/10:ticks*5],ticks,lonDataE,latDataE,12)
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-60,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],[units,' decade^-^1'],5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots0(h0,[.4 .4 .4],lonDataE,latDataE);
print(Fig,['E:\figures\IAP\Yearly\20240416_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_3.png'],'-dpng','-r300')
%% area mean Time series (anomaly)
lats = [126:146]; latstr = '30-55S';
% lats = [146:171]; latstr = '55-80S';

% var_r_0 = cat(1,nhflx(150:360,:,:),nhflx(1:149,:,:)); 
% str1 = 'Net Heat Flux (W m^-^2)';
% str2 = 'NHFLX';
% [spac_0] = areamean(var_r_0(:,:,21:81),1:151,lats,latDataE); % 1960-1970
% [sat_0] = areamean(var_r_0(:,:,21:81),152:231,lats,latDataE); 
% [sio_0] = areamean(var_r_0(:,:,21:81),232:360,lats,latDataE); 
% ln1 = spac_0;  max(ln1)
% ln2 = sat_0;  max(ln2)
% ln3 = sio_0; min(ln3)
% % ymax = 25; ymin = -5; yint = 5;
% % tt01 = 14; tt02 = 21;
% % tt11 = 27; tt12 = 21;
% % tt21 = 42; tt22 = 21;
% ymax = 20; ymin = -10; yint = 5;
% tt01 = 14; tt02 = 16;
% tt11 = 27; tt12 = 16;
% tt21 = 42; tt22 = 16;

% var_r_0 = cat(1,curlz_raw(150:360,:,:),curlz_raw(1:149,:,:)); 
% str1 = 'Wind Stress Curl (10 ^-^9 Pa m^-^1)';
% str2 = 'CURLZ';
% [spac_0] = areamean(var_r_0(:,:,21:81),1:151,lats,latDataE); % 1960-1970
% [sat_0] = areamean(var_r_0(:,:,21:81),152:231,lats,latDataE); 
% [sio_0] = areamean(var_r_0(:,:,21:81),232:360,lats,latDataE); 
% ln1 = spac_0*10^9;  max(ln1)
% ln2 = sat_0*10^9;  max(ln2)
% ln3 = sio_0*10^9; min(ln3)
% % ymax = 36; ymin = 6; yint = 6;
% % tt01 = 14; tt02 = 32;
% % tt11 = 27; tt12 = 32;
% % tt21 = 42; tt22 = 32;
% ymax = -6; ymin = -36; yint = 6;
% tt01 = 14; tt02 = -10;
% tt11 = 27; tt12 = -10;
% tt21 = 42; tt22 = -10;

var_r_0 = cat(1,taux_raw(150:360,:,:),taux_raw(1:149,:,:)); 
str1 = 'Zonal Wind Stress (10 ^-^3 Pa)';
str2 = 'TAUX';
[spac_0] = areamean(var_r_0(:,:,21:81),1:151,lats,latDataE); % 1960-1970
[sat_0] = areamean(var_r_0(:,:,21:81),152:231,lats,latDataE); 
[sio_0] = areamean(var_r_0(:,:,21:81),232:360,lats,latDataE); 
ln1 = spac_0*10^3;  max(ln1)
ln2 = sat_0*10^3; min(ln2)
ln3 = sio_0*10^3; min(ln3)
ymax = 51; ymin = 15; yint = 6;
tt01 = 14; tt02 = 46;
tt11 = 27; tt12 = 46;
tt21 = 42; tt22 = 46;
% ymax = 30; ymin = -6; yint = 6;
% tt01 = 14; tt02 = 25;
% tt11 = 27; tt12 = 25;
% tt21 = 42; tt22 = 25;

lonDara_r = cat(1,lonData(150:360),lonData(1:149));
close all;
bsln = '1960-1970';
sdiff = smooth(ln2-ln1);
Fig = figure('position',[700 50 600 400]);
% p3 = bar([1:61],sdiff,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7],'BarWidth',1);
% hold on
% set(gca,'YLim',[-18,18],'YTick',[ymin:yint:ymax],'YColor','k');
% ylabel('Difference  (K)','FontSize',14);
% yb = polyfit([1:61],sdiff,1);
% plot([1:61],polyval(yb,[1:61]),'--','Color',[0.3,0.3,0.3],'linewidth',2)
% text(58,0.3,[num2str(roundn(yb(1),-4)*10),' K per decade'],'fontsize',14,'Color',[0.3,0.3,0.3])
% 
p1 = plot(smooth(ln1),'-','color',[.04 .35 .56],'LineWidth',2);
hold on
p2 = plot(smooth(ln2),'-','color',[.64 .08 .18],'LineWidth',2);
p3 = plot(smooth(ln3),'-','color',[0.93,0.69,0.1],'LineWidth',2);
plot(zeros(1,61),'k','LineWidth',1)
yb = polyfit([1:61],ln1,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[.04 .35 .56],'linewidth',2)
text(tt01,tt02,[num2str(roundn(yb(1),-4))],'fontsize',12,'Color',[.04 .35 .56])
yb = polyfit([1:61],ln2,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[.64 .08 .18],'linewidth',2)
text(tt11,tt12,[num2str(roundn(yb(1),-4))],'fontsize',12,'Color',[.64 .08 .18])
yb = polyfit([1:61],ln3,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[0.93,0.69,0.13],'linewidth',2)
text(tt21,tt22,[num2str(roundn(yb(1),-4))],'fontsize',12,'Color',[0.93,0.69,0.13])
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel(str1,'FontSize',14);
% text(tt21,tt22,['Baseline ',bsln],'fontsize',12,'Color','k')
% legend([p2, p1, p3],'Atlantic-Indian Ocean','Pacific','Difference','Location','north','Orientation','horizontal')
legend([p1,p2,p3],'Pacific','Atlantic','Indian Ocean','Location','north','Orientation','horizontal')
legend('boxoff')
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',str2,'_TimeSeries_',latstr,'_',bsln,'.png'],'-dpng','-r300')
%% area mean Time series (raw)
% lats = [126:146]; % 55S-35S
% str3 = '35S_55S';
lats = [146:171]; % 80S-55S
str3 = '55S_80S';

var_r_0 = cat(1,nhflx(150:360,:,:),nhflx(1:149,:,:)); 
str1 = 'Net Heat Flux (W m^-^2)';
str2 = 'NHFLX';
[spacz_0 spac_0] = areamean(var_r_0(:,:,21:81),1:151,lats,latDataE); % 1960-1970
[siaz_0 sia_0] = areamean(var_r_0(:,:,21:81),152:360,lats,latDataE); 
ln1 = spac_0;  max(ln1)
ln2 = sia_0; min(ln2)
% ymax1 = 18; ymin1 = -3; yint1 = 3;
% ymax2 = 18; ymin2 = -3; yint2 = 3;
% tt01 = 35; tt02 = 15;
% tt11 = 15; tt12 = 15;
% tt21 = 50; tt22 = 15;
ymax1 = 21; ymin1 = -9; yint1 = 6;
ymax2 = 12; ymin2 = -18; yint2 = 6;
tt01 = 35; tt02 = 8;
tt11 = 15; tt12 = 8;
tt21 = 50; tt22 = 17;

% var_r_0 = cat(1,curlz_raw(150:360,:,:),curlz_raw(1:149,:,:)); 
% str1 = 'Wind Stress Curl (10 ^-^9 Pa m^-^1)';
% str2 = 'CURLZ';
% [spacz_0 spac_0] = areamean(var_r_0(:,:,21:81),1:151,lats,latDataE); % 1960-1970
% [siaz_0 sia_0] = areamean(var_r_0(:,:,21:81),152:360,lats,latDataE); 
% ln1 = spac_0*10^9;  max(ln1)
% ln2 = sia_0*10^9; min(ln2)
% ymax1 = 12; ymin1 = -9; yint1 = 3;
% ymax2 = 30; ymin2 = 9; yint2 = 3;
% tt01 = 35; tt02 = 27;
% tt11 = 15; tt12 = 27;
% tt21 = 50; tt22 = 9;
% ymax1 = 12; ymin1 = -18; yint1 = 6;
% ymax2 = -6; ymin2 = -36; yint2 = 6;
% tt01 = 35; tt02 = -10;
% tt11 = 15; tt12 = -10;
% tt21 = 50; tt22 = 8;

% var_r_0 = cat(1,taux_raw(150:360,:,:),taux_raw(1:149,:,:)); 
% str1 = 'Zonal Wind Stress (10 ^-^3 Pa)';
% str2 = 'TAUX';
% [spacz_0 spac_0] = areamean(var_r_0(:,:,21:81),1:151,lats,latDataE); % 1960-1970
% [siaz_0 sia_0] = areamean(var_r_0(:,:,21:81),152:360,lats,latDataE); 
% ln1 = spac_0*10^3;  max(ln1)
% ln2 = sia_0*10^3; min(ln2)
% ymax1 = 30; ymin1 = -10; yint1 = 5;
% ymax2 = 50; ymin2 = 10; yint2 = 5;
% tt01 = 35; tt02 = 44;
% tt11 = 15; tt12 = 44;
% tt21 = 50; tt22 = 24;
% ymax1 = 15; ymin1 = -25; yint1 = 5;
% ymax2 = 35; ymin2 = -5; yint2 = 5;
% tt01 = 35; tt02 = 29.5;
% tt11 = 15; tt12 = 29.5;
% tt21 = 50; tt22 = 9.5;

lonDara_r = cat(1,lonData(150:360),lonData(1:149));
close all;
bsln = '1960-1970';
sdiff = smooth(ln2-ln1);
Fig = figure('position',[700 50 600 400]);
p3 = bar([1:61],sdiff,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7],'BarWidth',1);
hold on
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint1:ymax1],'YColor','k');
ylabel('Difference','FontSize',14);
yb = polyfit([1:61],sdiff,1);
% plot([1:61],polyval(yb,[1:61]),'--','Color',[0.3,0.3,0.3],'linewidth',2)
text(tt21,tt22,[num2str(roundn(yb(1),-4))],'fontsize',12,'Color',[0.3,0.3,0.3])

yyaxis right
p1 = plot(smooth(ln1),'-','color',[.04 .35 .56],'LineWidth',2);
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:yint2:ymax2],'YColor','k');
hold on
p2 = plot(smooth(ln2),'-','color',[.64 .08 .18],'LineWidth',2);
plot(zeros(1,61),'k','LineWidth',1)
yb = polyfit([1:61],ln1,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[.04 .35 .56],'linewidth',2)
text(tt01,tt02,[num2str(roundn(yb(1),-4))],'fontsize',12,'Color',[.04 .35 .56])
yb = polyfit([1:61],ln2,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[.64 .08 .18],'linewidth',2)
text(tt11,tt12,[num2str(roundn(yb(1),-4))],'fontsize',12,'Color',[.64 .08 .18])
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel(str1,'FontSize',14);
legend([p2, p1, p3],'Atlantic-Indian Ocean','Pacific','Difference','Location','north','Orientation','horizontal')
legend('boxoff')
print(Fig,['E:\figures\IAP\Yearly\20240416_SO_warming\',str2,'_TimeSeries_',str3,'_raw.png'],'-dpng','-r300')


%% zonal mean
% name = 'NHFLX_zonalm';
% var_raw = nhflx;
% name = 'CURLZ_zonalm';
% var_raw = curlz_raw;
name = 'TAUX_zonalm';
var_raw = taux_raw;
var_r = cat(1,var_raw(151:360,:,:),var_raw(1:150,:,:));
varpac = permute(mean(var_r(1:151,:,:,:),1,'omitmissing'),[2 3 4 1]); % Pac unit:W
varatl = permute(mean(var_r(152:231,:,:,:),1,'omitmissing'),[2 3 4 1]); % Atl unit:W
vario = permute(mean(var_r(232:end,:,:,:),1,'omitmissing'),[2 3 4 1]); % IO unit:W
varall = permute(mean(var_r,1,'omitmissing'),[2 3 4 1]); % all
% trend
startyr = 1960;
endyr = 2020;
var1 = permute(varpac(:,startyr-1939:endyr-1939),[2 1]);
var2 = permute(varatl(:,startyr-1939:endyr-1939),[2 1]);
var3 = permute(vario(:,startyr-1939:endyr-1939),[2 1]);
x = [1:size(var1,1)]';
clear trd1
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd1(i) = par1(1); 
        par2=polyfit(x,var2(:,i),1); % regression parameters
        trd2(i) = par2(1); 
        par3=polyfit(x,var3(:,i),1); % regression parameters
        trd3(i) = par3(1); 
end
% h0 = trendtest(var1,0.05); % t test trend
max(trd2,[],'all')
min(trd2,[],'all')
%
climap1 = mean(varpac,2,'omitmissing');
climap2 = mean(varatl,2,'omitmissing');
climap3 = mean(vario,2,'omitmissing');
max(climap2,[],'all')
min(climap2,[],'all')
close all;
plot(latDataE,trd1,'r')
hold on
plot(latDataE,trd2,'b')
plot(latDataE,trd3,'y')
set(gca,'XLim',[-90,90],'xtick',[-90:10:90],'fontsize',12,'XGrid','on','YGrid','on');
%% nhflx plot
Fig = plotzm_nhflx(climap1,trd1,'Pacific',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_Pac_150E.png'],'-dpng','-r300')
%
Fig = plotzm_nhflx(climap2,trd2,'Atlantic',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_Atl_150E.png'],'-dpng','-r300')

Fig = plotzm_nhflx(climap3,trd3,'Indian Ocean',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_IO_150E.png'],'-dpng','-r300')
%% curlz plot
Fig = plotzm_curlz(climap1*10^7,trd1*10^9,'Pacific',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_Pac_150E.png'],'-dpng','-r300')
%
Fig = plotzm_curlz(climap2*10^7,trd2*10^9,'Atlantic',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_Atl_150E.png'],'-dpng','-r300')
%
Fig = plotzm_curlz(climap3*10^7,trd3*10^9,'Indian Ocean',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_IO_150E.png'],'-dpng','-r300')
%% taux plot
Fig = plotzm_taux(climap1*10^2,trd1*10^3,'Pacific',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_Pac_150E.png'],'-dpng','-r300')
%
Fig = plotzm_taux(climap2*10^2,trd2*10^3,'Atlantic',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_Atl_150E.png'],'-dpng','-r300')
%
Fig = plotzm_taux(climap3*10^2,trd3*10^3,'Indian Ocean',latDataE)
print(Fig,['E:\figures\IAP\Yearly\20240410_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_IO_150E.png'],'-dpng','-r300')

%%
unit1 = 'W';
unit2 = 'W decade^-^1';
close all;
Fig = figure('position',[700 100 600 400]);
ymax = 50; ymin = -40; yint = 10;
p1 = plot(latDataE(121:172),climap1(121:172),'--','color',[.04 .35 .56]','LineWidth',2);
hold on
p2 = plot(latDataE(121:172),climap2(121:172),'--','color',[.64 .08 .18]','LineWidth',2);
p3 = plot(latDataE(121:172),climap3(121:172),'--','color',[0.93,0.69,0.1],'LineWidth',2);
plot(latDataE(121:172),zeros(1,52),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12),ylabel(['Climo (',unit1,')'],'FontSize',14);

yyaxis right
ymax = 5; ymin = -4; yint = 1;
p4 = plot(latDataE(121:172),trd1(121:172)*10,'-','color',[.04 .35 .56],'LineWidth',2);
p5 = plot(latDataE(121:172),trd2(121:172)*10,'-','color',[.64 .08 .18],'LineWidth',2);
p6 = plot(latDataE(121:172),trd3(121:172)*10,'-','color',[0.93,0.69,0.1],'LineWidth',2);
plot(latDataE(126)*ones(10),[ymin:(ymax-ymin)/9:ymax],'-')
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
hold on
set(gca,'XLim',[-80,-30],'xtick',[-80:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o'});
ylabel(['Trend  (',unit2,')'],'FontSize',14);
% legend([p1, p2],'Climo (Atlantic-Indian Ocean)','Trend (Atlantic-Indian Ocean)','Location','northwest')
legend([p1, p2],'Climo (Pacific)','Trend (Pacific)','Location','northwest')
legend('boxoff')
% print(Fig,['E:\figures\IAP\Yearly\20240416_SO_warming\',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_Pac_150E.png'],'-dpng','-r300')
%% geostrophic MOHT
datadir='E:\data\IAP\geostrophic current\IAP_Ocean_geostrophic_current_0_2000m\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; % 查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE'); 
lonDatav = ncread(filetemp,'lon');
latDatav = ncread(filetemp,'lat');
depthDatav = ncread(filetemp,'lat');
k=length(filelist);
l = 0; 
clear vgeo0 vgeo
for s=1:k;  % 1940-2020
    
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');    
    if mod(s,12) ~= 0
        l=l+1
        vgeo0(:,:,:,l) = ncread(filename1,'v'); 
    else
        s/12
        vgeo0(:,:,:,l+1) = ncread(filename1,'v'); 
        vgeo(:,:,:,s/12) = nanmean(double(vgeo0),4); % yearly mean
        l = 0;
        clear vgeo0
        netcdf.close(ncid);  %关闭文件
    end
end
vgeo(:,:,:,82) = [];
vgeo(:,:,180,:) = vgeo(:,:,179,:);
Vvel = permute(vgeo,[2 3 1 4]);
%% MOHT
f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:end)-depthData(1:end-1);
levdist(1) = 4; levdist(2:41) = dweit; % z distance
TemVel = Temp(:,1:90,:,:).*Vvel(:,1:90,:,:);
OHCz = permute(ro*cp*sum((TemVel(:,1:90,1:41,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z
OHCzy = OHCz.*londist';
OHCzy_r = cat(1,OHCzy(180:300,:,:),OHCzy(301:360,:,:),OHCzy(1:179,:,:));
MOHTpac = permute(nansum(OHCzy_r(1:121,:,:,:),1),[2 3 4 1]); % Pac unit:W
MOHTia = permute(nansum(OHCzy_r(122:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
MOHTall = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
%% climo
close all;
map = nanmean(MOHTall,2);
plot(map)
%% trend
startyr = 1960;
endyr = 2020;
var1 = permute(MOHTpac(:,startyr-1939:endyr-1939),[2 1]);
var2 = permute(MOHTia(:,startyr-1939:endyr-1939),[2 1]);
x = [1:size(var1,1)]';
clear trd1
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd1(i) = par1(1); 
        par2=polyfit(x,var2(:,i),1); % regression parameters
        trd2(i) = par2(1); 
end
h0 = trendtest(var1,0.05); % t test trend
max(trd1,[],'all')
min(trd1,[],'all')
%%
climap1 = nanmean(MOHTpac,2);
climap2 = nanmean(MOHTia,2);
close all;
plot(trd1,'r')
hold on
plot(trd2,'b')
close all;
Fig = figure('position',[700 100 800 400]);
ymax = 8; ymin = -1; yint = 0.5;
p3 = plot(climap2/10^15,'-','color','k','LineWidth',3);
hold on
% p4 = plot(climap1,'-','color','k','LineWidth',3);
plot(zeros(1,81),'k','LineWidth',1)
%
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12),ylabel('Climo (10^1^5 W)','FontSize',14);


yyaxis right
ymax = 8; ymin = -5; yint = 0.5;
p1 = plot(trd2/10^12,'-','color',[.64 .08 .18],'LineWidth',3);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
hold on
% p2 = plot(trd1/10^12,'-','color',[.04 .35 .56],'LineWidth',3);
% set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.04 .35 .56]);
set(gca,'XLim',[10,70],'xtick',[10:10:90],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o'});
ylabel('Trend  (10^1^2 W)','FontSize',14);
% legend([p2, p3, p1],'Atlantic-Indian Ocean','Pacific','Difference','Location','northwest')
% legend('boxoff')
% print(Fig,['E:\figures\CESM\Yearly\LE\0_LEM\TimeSeries_MOHT_AtlIO',depthstr,'m_trend_1940_2005.png'],'-dpng','-r300')

%% Ekman MOHT calculate
% ro*cp*OHC*Vy
% Vy = My/ro = -taux/f/ro 
datadir='E:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
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
tauxraw = cat(1,fliplr(taux_raw(2:360,1:end-1,:)),fliplr(taux_raw(1,1:end-1,:)));
My = -tauxraw./f'; % Ekman 质量输运
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
Vy = My/ro;
Vy(:,85:96,:) = nan; % 赤道不适用ekman输送
londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
levdist(1) = 1; levdist(2:27) = dweit; % z distance
clear OHC OHCT MOHT
OHC = permute(ro*cp*sum((Temp(:,:,1:27,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3])/700;
OHCT = OHC.*Vy;
OHCTw = OHCT.*londist';
OHCTw_r = cat(1,OHCTw(160:300,:,:),OHCTw(301:360,:,:),OHCTw(1:159,:,:));
MOHT1 = permute(nansum(OHCTw_r(1:141,:,:),1),[2 3 1]); % Pac unit:W
MOHT2 = permute(nansum(OHCTw_r(142:360,:,:),1),[2 3 1]); % IO+Atl unit:W
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
% print(Fig,['E:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\MOHT_climotology.png'],'-dpng','-r300')
%% MOHT anomaly,detrend, and 8-yr filter
MOHT1a = MOHT1-mean(MOHT1,2);
MOHT2a = MOHT2-mean(MOHT2,2);
dT = 1; cf = 1/8;
clear MOHT1ad MOHT1adf MOHT2ad MOHT2adf
for i = 1:180
    MOHT1ad(:,i) = detrend(MOHT1a(i,:)'); % linear detrend
    MOHT1adf(:,i) = lanczosfilter(MOHT1ad(:,i),dT,cf,[],'low'); % 8 year filtered
    MOHT2ad(:,i) = detrend(MOHT2a(i,:)'); % linear detrend
    MOHT2adf(:,i) = lanczosfilter(MOHT2ad(:,i),dT,cf,[],'low'); % 8 year filtered
end
MOHT1adf(1:4,:) = []; MOHT1adf(end-3:end,:) = [];
MOHT2adf(1:4,:) = []; MOHT2adf(end-3:end,:) = [];
%%
taux_r = cat(1,taux(160:300,:,:),taux(301:360,:,:),taux(1:159,:,:));
lonmh1 = permute(nanmean(taux_r(1:141,:,:),1),[2 3 1]);
lonmh2 = permute(nanmean(taux_r(142:360,:,:),1),[2 3 1]);
index = DI; pngname = 'DI';
clear par7 h07 par8 h08 par7t h07t par8t h08t
for i = 1:size(MOHT1adf,2); 
    [par7(i),h07(i),t] = reg1_ttest(index,MOHT1adf(:,i),0.05,0);
    [par8(i),h08(i),t] = reg1_ttest(index,MOHT2adf(:,i),0.05,0);
    [par7t(i),h07t(i),t] = reg1_ttest(index,lonmh1(i,:)',0.05,0);
    [par8t(i),h08t(i),t] = reg1_ttest(index,lonmh2(i,:)',0.05,0);
end
%%
close all;
Fig = figure('position',[100 100 800 400]);
p1 = plot(latData,par7/10^15,'r-','linewidth',1.5);
hold on
p2 = plot(latData,par8/10^15,'b-','linewidth',1.5);
plot(-61*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
plot(-46*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
h07(1:30) = 0; h08(1:30) = 0;
scatter(latData(find(h07 == 1)),par7(find(h07 == 1))/10^15,'o','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h08 == 1)),par8(find(h08 == 1))/10^15,'o','MarkerEdgeColor','b','linewidth',1)
plot(latData,zeros(1,length(latData)),'k','linewidth',1.5)
% set(gca,'YLim',[-.025,.032],'YTick',[-.03:.01:.03]); % IPO
set(gca,'YLim',[-.02,.02],'YTick',[-.03:.01:.03]); % DI
set(gca,'YTick',yticks,'FontSize',14);
ylabel('MOHT (PW)')
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
% print(Fig,['E:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.zonalmean_MOHT_Ekman.png'],'-dpng','-r300')

yyaxis right
p3 = plot(latDataE(1:180),par7t,'-','color',[0.90,0.51,0.31],'linewidth',1.5); % ERA lat
p4 = plot(latDataE(1:180),par8t,'-','color',[0.17,0.73,0.69],'linewidth',1.5); % ERA lat
scatter(latDataE(find(h07t == 1)),par7t(find(h07t == 1)),'o','MarkerEdgeColor',[0.90,0.51,0.31],'linewidth',1)
scatter(latDataE(find(h08t == 1)),par8t(find(h08t == 1)),'o','MarkerEdgeColor',[0.17,0.73,0.69],'linewidth',1)
set(gca,'YLim',[-.0023,.0023],'YTick',[-.002:.001:.002],'ycolor','k');
set(gca,'YTick',yticks,'FontSize',14);
ylabel('Zonal Wind Stress (Pa)','color','k','rotation',270,'position',[-31.8 0 -1]);
print(Fig,['E:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.zonalmean_MOHT&Taux.png'],'-dpng','-r300')
%%
Fig = figure('position',[100 100 800 400]);
p1 = plot(latData,par7/10^15,'r-','linewidth',1.5);
hold on
p2 = plot(latData,par8/10^15,'b-','linewidth',1.5);
set(gca,'YLim',[-0.4,0.4],'FontSize',14);
legend([p1,p2],'Pacific MOHT','Atlantic & Indian Ocean MOHT','Location','north','Orientation','vertical','position',[0.49,0.75,0.35,0.12])
legend('boxoff')
print(Fig,['E:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\legend2.png'],'-dpng','-r300')

Fig = figure('position',[100 100 800 400]);
p3 = plot(latDataE(1:180),par7t,'-','color',[0.90,0.51,0.31],'linewidth',1.5); % ERA lat
hold on
p4 = plot(latDataE(1:180),par8t,'-','color',[0.17,0.73,0.69],'linewidth',1.5); % ERA lat
set(gca,'YLim',[-.003 .003],'FontSize',14);
legend([p3,p4],'Pacific Zonal Wind Stress','Atlantic & Indian Ocean Zonal Wind Stress','Location','north','Orientation','vertical','position',[0.49,0.75,0.35,0.12])
legend('boxoff')
print(Fig,['E:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\legend3.png'],'-dpng','-r300')
















function [adf] = caladf(var)
% var lon*lat*time
vara = double(var - mean(var(:,:,42:71),3)); % remove climatology from 1981 to 2010
varr = permute(vara,[3 1 2]);
dT = 1; cf = 1/8;
for i = 1:size(varr,2);
    for j = 1:size(varr,3);
        ad(:,i,j) = detrend(varr(:,i,j)); % linear detrend
        adf(:,i,j) = lanczosfilter(ad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
adf = permute(adf,[2 3 1]);
end
function [ad] = detrend3d(var)
    varr = permute(var,[3 1 2]);
    for i = 1:size(varr,2);
        for j = 1:size(varr,3);
            ad(:,i,j) = detrend(varr(:,i,j)); % linear detrend
        end
    end
    ad = permute(ad,[2 3 1]);
end
function [af] = filter3d(var,dT,cf)
% var lon*lat*time
varr = permute(var,[3 1 2]);
for i = 1:size(varr,2);
    for j = 1:size(varr,3);
        af(:,i,j) = lanczosfilter(varr(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
af = permute(af,[2 3 1]);
end
function [varar] = LoadNOAA(filename,varname)
% filename is where the file locates
% varname is name of the variable in raw file
% filename = 'E:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\miscSI-MO\prmsl.mon.mean.nc'; % 1836-2015
    ncid=netcdf.open(filename,'NOWRITE');
    ncdisp(filename);
    var = ncread(filename,varname);
    var = permute(mean(reshape(var(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
    vara = double(var - mean(var(:,:,:),3)); % remove climatology from 1981 to 2010
%42:71
    varar = permute(vara,[3 1 2]);
end
function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('E:/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function contourfSPolar(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',80,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('E:/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
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
    load('E:/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
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
function [ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(sum(sum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1,'omitmissing'),2,'omitmissing')/sum(cos(latData(lats)'/180*pi).*var2,'all','omitmissing'),size(var1,3),1);
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
lonmh1 = permute(nanmean(nhflx(160:300,:,:),1),[2 3 1]);
nhflx_r = cat(1,nhflx(300:360,:,:),nhflx(1:299,:,:));
lonmh2 = permute(nanmean(nhflx_r(1:220,:,:),1),[2 3 1]);
% reg with PC1
var1 = lonmh1';
var2 = lonmh2';
clear par11 h01 par12 h02
for i = 1:size(var1,2);
    [par11(i),h01(i),t] = reg1_ttest(index,var1(:,i),0.1,1); % filtered
    [par12(i),h02(i),t] = reg1_ttest(index,var2(:,i),0.1,1); % filtered
end
max(par11,[],'all')
min(par11,[],'all')
%
Fig = figure('position',[100 100 800 400]);
plot(latData,par11,'r','linewidth',1.5)
hold on
plot(latData,par12,'b','linewidth',1.5)
plot(latData,zeros(1,length(latData)),'k','linewidth',1.5)
plot(-61*ones(20),[ylimval(1):(ylimval(2)-ylimval(1))/19:ylimval(2)],'k--','linewidth',1.5);
plot(-46*ones(20),[ylimval(1):(ylimval(2)-ylimval(1))/19:ylimval(2)],'k--','linewidth',1.5);
scatter(latData(find(h01 == 1)),par11(find(h01 == 1)),'o','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h02 == 1)),par12(find(h02 == 1)),'o','MarkerEdgeColor','b','linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
set(gca,'YLim',ylimval);
set(gca,'YTick',yticks,'FontSize',14);
ylabel(titlename)
% legend('Pacific','Atlantic & Indian Ocean','Location','north','Orientation','vertical','position',[0.61,0.75,0.28,0.17],'edgecolor','w')
saveas(Fig,['E:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.png'])
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
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
end
function testdots1(h,clor,lonData,latData)
% used to dot significant area
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',2,'color',clor); % F-test dots
end
function testdots0(h,clor,lonData,latData)
% used to dot insignificant area
    dots = mod(find(h == 0),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 0)/length(lonData))),'.','markersize',2,'color',clor); % F-test dots
end
function h0 = trendtest(var,alp)
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
function trd = trend_cal_3D(var)
% var is time*d1*d2
    x = [1:size(var,1)]';
    clear trd
    for i = 1:size(var,2);
        for j = 1:size(var,3);
            par=polyfit(x,var(:,i,j),1); % regression parameters
            trd(i,j) = par(1);
        end
    end
end
function Fig = plotzm_nhflx(climap,trd,str,latDataE)
unit1 = 'W m^-^2 decade^-^1';
unit2 = 'W m^-^2';
close all;
Fig = figure('position',[700 100 450 400]);
ymax = 5; ymin = -4; yint = 1;
p1 = plot(latDataE(121:172),trd(121:172)*10,'-','color',[.64 .08 .18],'LineWidth',2);
hold on
plot(latDataE(126)*ones(10),[ymin:(ymax-ymin)/9:ymax],'k--','LineWidth',2)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
ylabel(['Trend  (',unit2,')'],'FontSize',14);

yyaxis right
ymax = 50; ymin = -40; yint = 10;
p2 = plot(latDataE(121:172),climap(121:172),'-','color','k','LineWidth',2);
hold on
plot(latDataE(121:172),zeros(1,52),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','k');
ylabel(['Climo (',unit1,')'],'FontSize',14);

yyaxis left
set(gca,'YColor',[.64 .08 .18]);
ylabel(['Trend  (',unit2,')'],'FontSize',14);


set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12);
set(gca,'XLim',[-80,-30],'xtick',[-80:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o'});
legend([p1, p2],['Trend (',str,')'],['Climo (',str,')'],'Location','northwest')
legend('boxoff')
end
function Fig = plotzm_curlz(climap,trd,str,latDataE)
unit1 = '10^-^9 Pa m^-^1 decade^-^1';
unit2 = '10^-^7 Pa m^-^1';
close all;
Fig = figure('position',[700 100 450 400]);
ymax = 6; ymin = -6; yint = 2;
p1 = plot(latDataE(121:172),trd(121:172)*10,'-','color',[.64 .08 .18],'LineWidth',2);
hold on
plot(latDataE(126)*ones(10),[ymin:(ymax-ymin)/9:ymax],'k--','LineWidth',2)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
ylabel(['Trend  (',unit2,')'],'FontSize',14);
yyaxis right
ymax = 1.2; ymin = -1.2; yint = .4;
p2 = plot(latDataE(121:172),climap(121:172),'-','color','k','LineWidth',2);
hold on
plot(latDataE(121:172),zeros(1,52),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','k');
ylabel(['Climo (',unit1,')'],'FontSize',14);

yyaxis left
set(gca,'YColor',[.64 .08 .18]);
ylabel(['Trend  (',unit2,')'],'FontSize',14);


set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12);
set(gca,'XLim',[-80,-30],'xtick',[-80:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o'});
legend([p1, p2],['Trend (',str,')'],['Climo (',str,')'],'Location','northwest')
legend('boxoff')
end
function Fig = plotzm_taux(climap,trd,str,latDataE)
unit1 = '10^-^3 Pa decade^-^1';
unit2 = '10^-^2 Pa';
close all;
Fig = figure('position',[700 100 450 400]);
ymax = 3; ymin = -1.8; yint = 0.6;
p1 = plot(latDataE(121:172),trd(121:172)*10,'-','color',[.64 .08 .18],'LineWidth',2);
hold on
plot(latDataE(126)*ones(10),[ymin:(ymax-ymin)/9:ymax],'k--','LineWidth',2)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
ylabel(['Trend  (',unit2,')'],'FontSize',14);
yyaxis right
ymax = 9; ymin = -5.4; yint = 1.8;
p2 = plot(latDataE(121:172),climap(121:172),'-','color','k','LineWidth',2);
hold on
plot(latDataE(121:172),zeros(1,52),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','k');
ylabel(['Climo (',unit1,')'],'FontSize',14);

yyaxis left
set(gca,'YColor',[.64 .08 .18]);
ylabel(['Trend  (',unit2,')'],'FontSize',14);


set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12);
set(gca,'XLim',[-80,-30],'xtick',[-80:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o'});
legend([p1, p2],['Trend (',str,')'],['Climo (',str,')'],'Location','northwest')
legend('boxoff')
end