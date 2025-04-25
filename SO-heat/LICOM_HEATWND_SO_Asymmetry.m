% era5 forcing Data 
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
load("/Volumes/Togo4T/1_matlab/MatData/LICOM/lonData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/LICOM/latData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/LICOM/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
filename1 = '/Volumes/Togo4T/data/LICOMpost/Tem_2000m.era5.HEAT.grid.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
lonData = ncread(filename1,'lon');
latData = ncread(filename1,'lat');
depthData = ncread(filename1,'lev');
Tem1 = ncread(filename1,'ts');
filename2 = '/Volumes/Togo4T/data/LICOMpost/Tem_2000m.era5.WND.grid.nc';
Tem2 = ncread(filename2,'ts');
Tem = Tem1+Tem2;
%% ------------------------- 0-700m ----------------------------------------
depthstr = '0-700';
dweit = -(depthData(2:21)-depthData(1:20));
lats = 1:60;
T100 = Tem(:,:,21,:)+(Tem(:,:,22,:)-Tem(:,:,21,:))*(700-depthData(21))/(depthData(22)-depthData(21));
Temsub = permute(nansum(cat(3,Tem(:,:,1,:)*5,Tem(:,:,2:21,:).*permute(dweit,[3 2 1 4]),T100*79),3)/700,[1 2 4 3]);
%% ------------------------- 0-2000m ----------------------------------------
depthstr = '0-2000';
dweit = -(depthData(2:25)-depthData(1:24));
lats = 1:60;
Temsub = permute(nansum(cat(3,Tem(:,:,1,:)*5,Tem(:,:,2:25,:).*permute(dweit,[3 2 1 4])),3)/700,[1 2 4 3]);
%% meridional mean
lats = 35:55; % 55S-35S
latstr = '35S_55S';
[Tempa_mm] = latmean(Tem,lats,latData);
% ts is lon*depth*time
startyr = 1960;
endyr = 2020;
yrstr = '1960-2020';
var = permute(Tempa_mm(:,:,startyr-1940:endyr-1940),[3 1 2]);
x = [1:size(var,1)]';
clear trd
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.01); % t test trend
max(trd,[],'all')
min(trd,[],'all')
%% 150 longitude
close all;
ftsz = 20;
ticks = 0.1;
map = trd*10;
map_r = cat(1,map(150:360,:),map(1:150,:));
h0_r = cat(1,h0(150:360,:),h0(1:150,:));
lonData_r = [lonData(150:360);lonData(1:150)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,depthData,map_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([290 290],[0 -2000],'linewidth',1.5,'color','k')
line([150 509],[-700 -700],'linewidth',1.5,'color','k','linestyle','--')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(flipud(bl_re4));
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','K decade^-^1','Units','normalized','position',[6.5 0.5 0],'Rotation',-90,'fontsize',ftsz);
[i,j] = find(h0_r == 0); % significance test
plot(lonData_r(i(1:3:end)),depthData(j(1:3:end)),'.','MarkerSize',5,'color',[.4 .4 .4])
set(gca,'XLim',[150,360+150]);
set(gca,'XTick',[150:60:360+150]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',ftsz);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:500:0],'FontSize',ftsz);
set(gca,'YTickLabel',[2000:-500:0],'FontSize',ftsz);
ylabel('Depth (m)')
print(Fig,['/Users/yysong/Desktop/figures/LICOM/Yearly/20240705_SO_warming/WNDHEAT_Tem.meridionalmean_',latstr,'_trend_',yrstr,'_150E.png'],'-dpng','-r300')
%% OHC
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
%------------------------- 0-2000m -----------------------------------------
% depthstr = '0-2000m';
% dweit = depthData(2:41)-depthData(1:40); % depth weight
% OHCsub = cp*ro*permute(nansum(cat(3,Temp(:,:,1,:),Temp(:,:,2:41,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]); % J/m2
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
depthData_r = -depthData;
dweit = depthData_r(2:21)-depthData_r(1:20); % depth weight
Tem700 = Tem(:,:,21,:)+(Tem(:,:,22,:)-Tem(:,:,21,:))*(700-depthData_r(21))/(depthData_r(22)-depthData_r(21));
OHCsub = cp*ro*permute(nansum(cat(3,Tem(:,:,1,:)*5,Tem(:,:,2:21,:).*(permute(dweit,[3 2 1])),Tem700*(700-depthData_r(21))),3),[1 2 4 3]); % J/m2

OHCw = OHCsub/365/24/60/60; % W/m2
%% OHC trend
startyr = 1960;
endyr = 2020;
var = permute(OHCw(:,:,startyr-1940:endyr-1940),[3 1 2]);
x = [1:size(var,1)]';
clear trd 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.01); % t test trend
h0(find(isnan(Tem(:,:,1,1)) == 1)) = nan;
max(trd,[],'all')
min(trd,[],'all')
%% SO plot
close all;
mapr = cat(1,trd,trd(1,:));
lonr = [lonData;lonData(1)];
ftsz = 20; ticks = 1.5;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(mapr,[-ticks*5:ticks/10:ticks*5],ticks,lonr,latData,ftsz)
m_coast('patch','w');
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-60,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2',5.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots0(h0(1:2:end,1:2:end),[.4 .4 .4],lonData(1:2:end),latData(1:2:end));
print(Fig,['/Users/yysong/Desktop/figures/LICOM/Yearly/20240705_SO_warming/WNDHEAT_OHC',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.3.png'],'-dpng','-r300')
%% original OHC time series
lats = [35:55]; % 55S-35S
levs = [1:22]; % 700 m
lons = 150; % 
VV = ones(length(lonData),length(latData),length(depthData_r)); % calculate volume
Tem = double(Tem);
Temp700 = Tem(:,:,21,:)+(Tem(:,:,22,:)-Tem(:,:,21,:))*(700-depthData_r(21))/(depthData_r(22)-depthData_r(21));
Tempz = cat(3,Tem(:,:,1,:)*5,Tem(:,:,2:21,:).*(permute(dweit,[3 2 1])),Temp700*(700-depthData_r(21)));
VV700 = VV(:,:,21,:)+(VV(:,:,22,:)-VV(:,:,21,:))*(700-depthData_r(21))/(depthData_r(22)-depthData_r(21));
VVz = cat(3,VV(:,:,1)*5,VV(:,:,2:21).*(permute(dweit,[3 2 1])),VV700*(700-depthData_r(21)));
Tempzy = Tempz*111*1000; % 111 km / latitude
VVzy = VVz*111*1000;
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(lons+1:360,:,:,:),Tempzyx(1:lons,:,:,:)); 
VVzyx(find(isnan(Tempzyx(:,:,:,1)) == 1)) = nan;
Vzyxsub_r_0 = cat(1,VVzyx(lons+1:360,:,:),VVzyx(1:lons,:,:)); 
spac_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:300-lons,lats,levs,:),1),2),3));
spacS = squeeze(nansum(nansum(Vzyxsub_r_0(1:300-lons,lats,1),1),2)); % Square
spac0S = spac_0./spacS;
sia_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(300-lons+1:360,lats,levs,:),1),2),3));
siaS = squeeze(nansum(nansum(Vzyxsub_r_0(300-lons+1:360,lats,1),1),2)); % Square
sia0S = sia_0./siaS;
%%
% IAP square
datadir='/Volumes/Togo4T/data/IAP/temperature/Yearly/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
Temp_i = h5read(filetemp,'/temp in degC');
aa = load("/Volumes/Togo4T/1_matlab/MatData/IAP/lonData.mat"); lonData_i = aa.lonData;
aa = load("/Volumes/Togo4T/1_matlab/MatData/IAP/latData.mat"); latData_i = aa.latData;
aa = load("/Volumes/Togo4T/1_matlab/MatData/IAP/depthData.mat"); depthData_i = aa.depthData;
dweit_i = depthData_i(2:27)-depthData_i(1:26);
VV = ones(length(lonData_i),length(latData_i),length(depthData_i)); % calculate volume
VVz = cat(3,VV(:,:,1),VV(:,:,2:27).*(permute(dweit_i,[3 2 1])));
VVzy = VVz*111*1000; % 111km /latitude
clear dx Tempzyx VVzyx
for j = 1:length(latData_i);
    dx(j) = sw_dist([latData_i(j),latData_i(j)],[lonData_i(1),lonData_i(2)],'km'); % 每个纬度上，经度之间距离不一样
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
VVzyx(find(isnan(Temp_i(:,:,1:27,end)) == 1)) = nan;
Vzyxsub_r_0 = cat(1,VVzyx(lons+1:360,:,:),VVzyx(1:lons,:,:)); 
spacS_i = squeeze(nansum(nansum(Vzyxsub_r_0(1:300-lons,lats,1),1),2)); % Square
siaS_i = squeeze(nansum(nansum(Vzyxsub_r_0(300-lons+1:360,lats,1),1),2)); % Square

%% 
close all;
ln1 = (spac_0(20:80)-mean(spac_0(20:30)))/spacS_i/10^9; % 1960-2020
ln2 = (sia_0(20:80)-mean(sia_0(20:30)))/siaS_i/10^9; % 1960-2020
sdiff = ln2-ln1;
ftsz=20;
Fig = figure('position',[700 100 800 400]);
p1 = bar([1:61],sdiff,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7],'BarWidth',1);
hold on
yb = polyfit([1:61],sdiff,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[0.3,0.3,0.3],'linewidth',2)
text(43,1.25,[num2str(roundn(yb(1),-4))],'fontsize',ftsz,'Color',[0.3,0.3,0.3])

ymax = 1.5; ymin = -.3; yint = 0.3;
p2 = plot(ln2,'-','color',[.64 .08 .18],'LineWidth',3);
hold on
p3 = plot(ln1,'-','color',[.04 .35 .56],'LineWidth',3);
plot(zeros(1,81),'k','LineWidth',1)
yb = polyfit([1:61],ln2,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[.64 .08 .18],'linewidth',2)
text(20,1.25,[num2str(roundn(yb(1),-4))],'fontsize',ftsz,'Color',[.64 .08 .18])
yb = polyfit([1:61],ln1,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[.04 .35 .56],'linewidth',2)
text(33,1.25,[num2str(roundn(yb(1),-4))],'fontsize',ftsz,'Color',[.04 .35 .56])
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',ftsz);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'ycolor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',ftsz),ylabel(' (GJ m^-^2)','FontSize',ftsz);
legend([p2, p3, p1],'Atlantic-Indian Ocean','Pacific','Difference','Location','north','Orientation','horizontal')
legend('boxoff')
% print(Fig,['/Users/yysong/Desktop/figures/LICOM/Yearly/20240615_SO_warming/CTR_TimeSeries_OHC_GJm2_',depthstr,'_35S_55S_160W.png'],'-dpng','-r300')
%% net heat flux
startyr = 1960; endyr = 2020;
[nhflx lonData latData] = read_1lev_1nc('/Volumes/Togo4T/data/LICOMpost/CTLyrmean/net1_yrmean.nc','net1');
var = nhflx; varname = 'NHFLX'; 
zlim = 4; unit = 'W m^-^2 decade^-^1';
[trd h0] = lntrend3d(var(:,:,startyr-1940:endyr-1940),0.05);
max(trd,[],'all')
min(trd,[],'all')
[Fig] = pltrdpolar(trd*10,zlim,lonData,latData,20,unit)
print(Fig,['/Users/yysong/Desktop/figures/LICOM/Yearly/20240615_SO_warming/CTR_',varname,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
% [Fig] = pltrdpolar(mean(var(:,:,startyr-1940:endyr-1940),3),100,lonData,latData,12,unit)
% print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_climo_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% latent heat flux
startyr = 1960; endyr = 2020;
[lthf lonData latData] = read_1lev_1nc('/Volumes/Togo4T/data/LICOMpost/CTLyrmean/lthf_yrmean.nc','lthf');
var = lthf; varname = 'LTHF'; 
zlim = 4; unit = 'W m^-^2 decade^-^1';
[trd h0] = lntrend3d(var(:,:,startyr-1940:endyr-1940),0.05);
max(trd,[],'all')
min(trd,[],'all')
[Fig] = pltrdpolar(trd*10,zlim,lonData,latData,12,unit)
print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
% [Fig] = pltrdpolar(mean(var(:,:,startyr-1940:endyr-1940),3),100,lonData,latData,12,unit)
% print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_climo_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% sensible heat flux
startyr = 1960; endyr = 2020;
[sshf lonData latData] = read_1lev_1nc('/Volumes/Togo4T/data/LICOMpost/CTLyrmean/sshf_yrmean.nc','sshf');
var = sshf; varname = 'SSHF'; 
zlim = 4; unit = 'W m^-^2 decade^-^1';
[trd h0] = lntrend3d(var(:,:,startyr-1940:endyr-1940),0.05);
max(trd,[],'all')
min(trd,[],'all')
[Fig] = pltrdpolar(trd*10,zlim,lonData,latData,12,unit)
print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
% [Fig] = pltrdpolar(mean(var(:,:,startyr-1940:endyr-1940),3),100,lonData,latData,12,unit)
% print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_climo_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% long wave
startyr = 1960; endyr = 2020;
[lwv lonData latData] = read_1lev_1nc('/Volumes/Togo4T/data/LICOMpost/CTLyrmean/lwv_yrmean.nc','lwv');
var = lwv; varname = 'LWV'; 
zlim = 4; unit = 'W m^-^2 decade^-^1';
[trd h0] = lntrend3d(var(:,:,startyr-1940:endyr-1940),0.05);
max(trd,[],'all')
min(trd,[],'all')
[Fig] = pltrdpolar(trd*10,zlim,lonData,latData,12,unit)
print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
% [Fig] = pltrdpolar(mean(var(:,:,startyr-1940:endyr-1940),3),100,lonData,latData,12,unit)
% print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_climo_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% short wave
startyr = 1960; endyr = 2020;
[swv lonData latData] = read_1lev_1nc('/Volumes/Togo4T/data/LICOMpost/CTLyrmean/swv_yrmean.nc','swv');
var = swv; varname = 'SWV'; 
zlim = 4; unit = 'W m^-^2 decade^-^1';
[trd h0] = lntrend3d(var(:,:,startyr-1940:endyr-1940),0.05);
max(trd,[],'all')
min(trd,[],'all')
[Fig] = pltrdpolar(trd*10,zlim,lonData,latData,12,unit)
print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
% [Fig] = pltrdpolar(mean(var(:,:,startyr-1940:endyr-1940),3),100,lonData,latData,12,unit)
% print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_climo_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% sum of 4 heat
var = swv+lwv+lthf+sshf; varname = 'SUMheat'; 
zlim = 4; unit = 'W m^-^2 decade^-^1';
[trd h0] = lntrend3d(var(:,:,startyr-1940:endyr-1940),0.05);
max(trd,[],'all')
min(trd,[],'all')
[Fig] = pltrdpolar(trd*10,zlim,lonData,latData,12,unit)
% print(Fig,['D:/figures/LICOM/Yearly/20240425_SO_warming/era5.CTR.',varname,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% zonal plot  consistent with CESM LE
%%%%%%%%    net heat flux
lonstr = 151;
startyr = 1960;
endyr = 2020;
var_r = cat(1,nhflx(lonstr:360,:,:),nhflx(1:lonstr,:,:));
var1 = squeeze(nanmean(var_r(1:140,1:90,startyr-1939:endyr-1939),1))';
var2 = squeeze(nanmean(var_r(141:360,1:90,startyr-1939:endyr-1939),1))';
cli1 = nanmean(var1,1);
cli2 = nanmean(var2,1);
x = [1:size(var1,1)]';
clear trd* 
for i = 1:size(var1,2);
    par=polyfit(x,var1(:,i),1); % regression parameters
    trd1(i) = par(1);
    par=polyfit(x,var2(:,i),1); % regression parameters
    trd2(i) = par(1);
end
%%
close all;
ftsz = 20;
Fig = figure('position',[700 100 800 400]);
ymin1 = -.24; ymax1 = .32; yint1 = .08;
ymin2 = -39; ymax2 = 52; yint2 = 13;
p1 = plot(latData(20:60),trd1(20:60),'-','color',[.04 .35 .56],'LineWidth',3);
hold on
p2 = plot(latData(20:60),trd2(20:60),'-','color',[.64 .08 .18],'LineWidth',3);
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint1:ymax1],'YColor','k');
set(gca,'XLim',[-70.5,-30.5],'xtick',[-70.5:10:-30.5],'fontsize',ftsz,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('Trend  (W/yr)','FontSize',ftsz);
plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
p3 = plot(latData(20:60),cli1(20:60),'--','color',[.04 .35 .56],'LineWidth',1.5);
hold on
p4 = plot(latData(20:60),cli2(20:60),'--','color',[.64 .08 .18],'LineWidth',1.5);
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:yint2:ymax2],'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',ftsz),ylabel('Climo (W)','FontSize',ftsz);
legend([p1, p3, p2, p4],['Pacific trend'],['Pacific climo'],['Atlantic-Indian Ocean trend'],...
    ['Atlantic-Indian Ocean climo'],'Location','northwest','numcolumns',2,'edgecolor','w')
% 手动保存 20240615
%% MOHT
% cal vt
[ts] = read_nlev_1nc('/Volumes/Togo4T/data/LICOMpost/Tem_2000m.era5.grid.nc','ts');
[vs] = read_nlev_1nc('/Volumes/Togo4T/data/LICOMpost/CTLyrmean/vs_yrmean.nc','vs');
vs = vs(:,:,:,1:82);

inlev = 10; depthstr = '0-100m';
splon = 150; lonstr = '150E';
f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dweit = -(depthData(2:inlev)-depthData(1:inlev-1));
clear levdist
levdist(1) = 5; levdist(2:inlev) = dweit; levdist(inlev+1) = 100-(-depthData(inlev));% z distance
T100 = ts(:,:,10,:)+(ts(:,:,11,:)-ts(:,:,10,:))*(levdist(11))/(-(depthData(11)-depthData(10)));
Temp = cat(3,ts(:,:,1:10,:),T100);
v100 = vs(:,:,10,:)+(vs(:,:,11,:)-vs(:,:,10,:))*(levdist(11))/(-(depthData(11)-depthData(10)));
Vvel = cat(3,vs(:,:,1:10,:),v100);
%% 0-100m MOHT
TemVel = Temp(:,1:90,:,:).*Vvel(:,1:90,:,:); % 1960-
OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z
OHCzy = OHCz.*londist';
OHCzy_r = cat(1,OHCzy(splon:360,:,:,:),OHCzy(1:splon-1,:,:,:));
MOHTp = permute(nansum(OHCzy_r(1:290-splon,:,:,:),1),[2 3 4 1]); % Pac unit:W 150E-70w
MOHTia = permute(nansum(OHCzy_r(291-splon:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
MOHTi = permute(nansum(OHCzy_r(382-splon:end,:,:,:),1),[2 3 4 1]); % io 20E-150E unit:W
MOHTa = permute(nansum(OHCzy_r(291-splon:381-splon,:,:,:),1),[2 3 4 1]); % Atl 60W-20E unit:W

MOHTcp = nanmean(MOHTp,2); % climo
MOHTcia = nanmean(MOHTia,2);
MOHTci = nanmean(MOHTi,2);
MOHTca = nanmean(MOHTa,2);
startyr = 1960;
endyr = 2020;
var1 = permute(MOHTp(:,startyr-1939:endyr-1939),[2 1]);
var2 = permute(MOHTia(:,startyr-1939:endyr-1939),[2 1]);
cli1 = nanmean(MOHTp,2);
cli2 = nanmean(MOHTia,2);
x = [1:size(var1,1)]';
clear trd*
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd1(i) = par1(1); 
        par2=polyfit(x,var2(:,i),1); % regression parameters
        trd2(i) = par2(1); 
end

%% MOHT/dlon -> W/m-1

elon = 1000*sw_dist([45,45],[1,2],'km');
map1 = trd1/(140*elon); 
map2 = trd2/(210*elon); 
map1c = cli1/(140*elon);
map2c = cli2/(210*elon);
% map = trd3/(90*elon); regstr = 'IO_30E-120E'
close all;
ftsz = 20;
Fig = figure('position',[700 100 800 400]);
ymax1 = 200; ymin1 = -50; yint1 = 25; 
ymax2 = 100; ymin2 = -50; yint2 = 25; 
p1 = plot(latData(20:60),map1(20:60)/10^3,'-','color',[.04 .35 .56],'LineWidth',3);
hold on
p2 = plot(latData(20:60),map2(20:60)/10^3,'-','color',[.64 .08 .18],'LineWidth',3);
plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint1:ymax1]);
ylabel('Trend  (kW m^-^1 yr^-^1)','FontSize',ftsz);

yyaxis right
p3 = plot(latData(20:60),map1c(20:60)/10^6,'--','color',[.04 .35 .56],'LineWidth',1.5);
hold on
p4 = plot(latData(20:60),map2c(20:60)/10^6,'--','color',[.64 .08 .18],'LineWidth',1.5);
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:yint2:ymax2],'ycolor','k');
ylabel('Climo  (10^3 kW m^-^1)','FontSize',ftsz,'Color','k');
set(gca,'XLim',[-70,-30],'xtick',[-70:10:-30],'fontsize',ftsz,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
set(gca,'XGrid','on','YGrid','on');
legend([p1, p3, p2, p4],['Pacific trend'],['Pacific climo'],['Atlantic-Indian Ocean trend'],...
    ['Atlantic-Indian Ocean climo'],'Location','northwest','numcolumns',2,'edgecolor','w')

% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240503_SO_Warming/MOHT_VRES_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
% 手动哦保存
%% dMOHT/dy m-1
ymin1 = -3.9; ymax1 = 3.9;  ymin2 = -900; ymax2 = 900;
%% dMOHT/dy m-1

dMdyp = dMOHTdy(MOHTp(lats,startyr-1939:endyr-1939));
dMdyia = dMOHTdy(MOHTia(lats,startyr-1939:endyr-1939));
cli1 = nanmean(dMdyp,2);
cli2 = nanmean(dMdyia,2);
%
var1 = dMdyp';
var2 = dMdyia';
x = [1:size(var1,1)]';
clear trd1 trd2
for i = 1:size(var1,2);
    par1=polyfit(x,var1(:,i),1); % regression parameters
    trd1(i) = par1(1);
    par2=polyfit(x,var2(:,i),1); % regression parameters
    trd2(i) = par2(1);
end
%
close all;
figure(1)
plot(latData(lats),trd1,'r')
hold on
plot(latData(lats),trd2,'b')
legend('pac','atl','IO')
%% MOHT/dlon -> W/m-1
elon = 1000*sw_dist([45,45],[1,2],'km');

% Pacific
map1 = trd1/(140*elon); 
map2 = cli1/(140*elon);

% Atlantic-Indian Ocean
map1 = trd2/(210*elon); 
map2 = cli2/(210*elon);

ymin1 = -.24; ymax1 = .24; yint = .08
ymin2 = -90; ymax2 = 90; yint2 = 30
latd = latData(lats);
close all;
ftsz = 20;
Fig = figure('position',[100 100 800 400])
nump = find(map1 > 0); numn = find(map1 < 0);
h1 = bar(latd(nump),map1(nump),'FaceColor',[.64 .08 .18],'EdgeColor',[.64 .08 .18],'barwidth',1)
hold on
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint:ymax1]);
h2 = bar(latd(numn),map1(numn),'FaceColor',[.04 .35 .56],'EdgeColor',[.04 .35 .56],'barwidth',1)
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
plot(-35.5*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55.5*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
set(gca,'XLim',[-70.5,-30.5],'XTick',[-70.5:10:-30.5],'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',ftsz,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',ftsz),ylabel('Trend (W m^-^2 yr^-^1)','FontSize',ftsz);

yyaxis right
numpc = find(map2 > 0); numnc = find(map2 < 0);
h3 = bar(latd(numpc),map2(numpc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
h4 = bar(latd(numnc),map2(numnc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:yint2:ymax2],'ycolor','k');
ylabel('Climo (W m^-^2)','FontSize',ftsz);
legend([h1, h2],'Divergence trend','Convergence trend','Location','northwest','edgecolor','w')

% 手动保存 MOHTdydx_0-100m_Pac_1960-2020.png


function [dMdy] = dMOHTdy(var)
[d1 d2] = size(var);
dy = sw_dist([-90,-89],[0,0],'km')*10^3; % m
dMdy(1,:) = nan;
dMdy(d1,:) = nan;
for t = 1:d2
    for i = 2:d1-1
        dMOHT(i,t) = var(i+1,t)-var(i-1,t);
        dMdy(i,t) = dMOHT(i,t)/(2*dy);

    end
end
end
function [Fig] = pltrdpolar(mval,zlim,lon,lat,ftsz,clbarstr)
close all;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(mval,[-zlim*5:zlim/10:zlim*5],zlim,lon,lat,ftsz)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-55,'linewidth',2,'color','k'); 
m_line(-60,[-55:1:-35],'linewidth',2,'color','k'); 
m_line(150,[-55:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],clbarstr,5.5,ftsz,[-zlim:zlim/5:zlim],[-zlim:zlim/5:zlim])
% testdots(h0,[.4 .4 .4],lonData,latData);
end
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
% filename = 'E:/data/NOAA/20CRv3/20thC_ReanV3/Monthlies/miscSI-MO/prmsl.mon.mean.nc'; % 1836-2015
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
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','none','xaxisloc','top');
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
saveas(Fig,['E:/figures/IAP/Yearly/20230606_IPO_SouthernOcean_46_61S/',pngname,'.png'])
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
function h0 = trendtest(var,alp)
% var is time*lon*lat
d1 = size(var,1);
d2 = size(var,2);
d3 = size(var,3);
h0 = nan(d2,d3);
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
        elseif abs(T) < tinv(1-alp/2,d1-2); % 双侧 t test
            h0(i,j) = 0;
        end
    end
end
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
    m_plot(lonData(dots),latData(ceil(find(h == 0)/length(lonData))),'.','markersize',5,'color',clor); % F-test dots
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