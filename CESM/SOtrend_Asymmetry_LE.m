% Data 
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
load("MatFile/lonData.mat");
load("MatFile/latData.mat");
load("MatFile/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
%% LE ensmean (EX)
filename1 = '/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_1920-2005.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
TemMMM1 = ncread(filename1,'TEMP');
depthData = ncread(filename1,'z_t')/100; % cm -> m
filename2 = '/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_2006-2022.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
TemMMM2 = ncread(filename2,'TEMP');
TemMMMall = cat(4,TemMMM1,TemMMM2(:,:,1:47,:));
Temp = TemMMMall;
clear TemMMM1 TemMMM2
%%  Tem
% ------------------------- 0-700m ----------------------------------------
depthstr = '0-700';
dweit = depthData(2:37)-depthData(1:36);
lats = 1:60;
TemMMMsub = permute(nansum(cat(3,TemMMMall(:,:,1,:)*5,TemMMMall(:,:,2:37,:).*permute(dweit,[3 2 1 4])),3)/707,[1 2 4 3]);
% ------------------------- 0-2000m ----------------------------------------
% depthstr = '0-2000';
% dweit = depthData(2:47)-depthData(1:46);
% lats = 1:60;
% TemMMMsub = permute(nansum(cat(3,TemMMMall(:,:,1,:)*5,TemMMMall(:,:,2:47,:).*permute(dweit,[3 2 1 4])),3)/2000,[1 2 4 3]);

Tsub_r = cat(1,TemMMMsub(182:360,lats,:),TemMMMsub(1:181,lats,:));
lonData_r = [lonData(182:360);lonData(1:181)+360];
%% Tem trend 2d
startyr = 2005;
endyr = 2020;
var = permute(Tsub_r(:,:,startyr-1919:endyr-1919),[3 1 2]);
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
%% Tem plot
close all;
ftsz = 12; ticks = 0.25;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd*10,[-ticks*5:ticks/10:ticks*5],ticks,lonData_r,latData(lats),12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-50,'linewidth',2,'color','k'); 
m_line(-60,[-50:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-50:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K decade^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots(h0,[.4 .4 .4],lonData_r,latData(lats));
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/Tem',depthstr,'m_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% Time series
lats = [40:55]; % latitude
[spacz_long spac_long] = areamean(Tsub_r,1:120,lats,latData); 
[siaz_long sia_long] = areamean(Tsub_r,121:360,lats,latData); 
DI_rawLE = sia_long-spac_long;
%% meridional mean
lats = 35:55; % 55S-35S
latstr = '35S_55S';
[Tempa_mm] = latmean(Temp,lats,latData);
% ts is lon*depth*time
startyr = 1960;
endyr = 2020;
yrstr = '1960-2020';
var = permute(Tempa_mm(:,:,startyr-1939:endyr-1939),[3 1 2]);
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
%% 180 longitude
close all;
ftsz = 12;
ticks = 0.1;
map = trd*10;
map_r = cat(1,map(180:360,:),map(1:180,:));
lonData_r = [lonData(180:360);lonData(1:180)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,map_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([300 300],[0 -2000],'linewidth',1.5,'color','k')
% line([150 509],[-700 -700],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(flipud(bl_re4));
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','K decade^-^1','Units','normalized','position',[4.5 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[180,360+180]);
set(gca,'XTick',[180:60:360+180]);
set(gca,'XTickLabel',{'180^oE','120^oW','60^oW','0^o','60^oE','120^oE','180^oE'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:400:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-400:0],'FontSize',14);
ylabel('Depth (m)')
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/Tem_meridionalmean_',latstr,'_trend_',yrstr,'.png'],'-dpng','-r300')
%% 150 longitude
close all;
ftsz = 12;
ticks = 0.1;
map = trd*10;
map_r = cat(1,map(150:360,:),map(1:150,:));
lonData_r = [lonData(150:360);lonData(1:150)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,map_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([290 290],[0 -2000],'linewidth',1.5,'color','k')
% line([150 509],[-700 -700],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
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
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/meridionalmean_',latstr,'_trend_',yrstr,'_150E.png'],'-dpng','-r300')
%% OHC
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
%------------------------- 0-2000m -----------------------------------------
% depthstr = '0-2000m';
% dweit = depthData(2:46)-depthData(1:45); % depth weight
% Temp2000 = Temp(:,:,46,:)+(Temp(:,:,47,:)-Temp(:,:,46,:))*(2000-depthData(46))/(depthData(47)-depthData(46));
% OHCsub = cp*ro*permute(nansum(cat(3,TemMMMall(:,:,1,:)*5,TemMMMall(:,:,2:46,:).*(permute(dweit,[3 2 1])),Temp2000*(2000-depthData(46))),3),[1 2 4 3]); % J/m2
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthData(2:36)-depthData(1:35); % depth weight
Temp700 = Temp(:,:,36,:)+(Temp(:,:,37,:)-Temp(:,:,36,:))*(700-depthData(36))/(depthData(37)-depthData(36));
OHCsub = cp*ro*permute(nansum(cat(3,TemMMMall(:,:,1,:)*5,TemMMMall(:,:,2:36,:).*(permute(dweit,[3 2 1])),Temp700*(700-depthData(36))),3),[1 2 4 3]); % J/m2
OHCw = OHCsub/365/24/60/60; % W/m2
%% Delta OHC
yrspan = 15;
yrs1 = 1970; 
yrs2 = 2005; 

map = nanmean(OHCw(:,:,yrs2-1919:yrs2-1919+yrspan-1),3)-nanmean(OHCw(:,:,yrs1-1919:yrs1-1919+yrspan-1),3);
max(map,[],'all')
min(map,[],'all')
% SO asymmetry
close all;
ftsz = 12; ticks = 60;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS((map),[-ticks*5:ticks/10:ticks*5],ticks,lonData,latData,12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
% m_line([0:1:360],-50,'linewidth',2,'color','k'); 
m_line(-60,[-90:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-90:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'/DeltaOHC  (W m^-^2)',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots(h0,[.4 .4 .4],lonData,latData);
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/Delta_OHC_W_raw',depthstr,'_',num2str(yrs1),'_',num2str(yrs2),'_',num2str(yrspan),'yrs.png'],'-dpng','-r300')
%% OHC trend
startyr = 1980;
endyr = 2020;
var = permute(OHCw(:,:,startyr-1919:endyr-1919),[3 1 2]);
x = [1:size(var,1)]';
clear trd 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
%
h0 = trendtest(var,0.05); % t test trend
h0(find(isnan(TemMMMall(:,:,1,1)) == 1)) = nan;
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
m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2 yr^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots0(h0,[.4 .4 .4],lonData,latData);
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/OHC',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.3.png'],'-dpng','-r300')
%% original OHC time series
lats = [35:55]; % 55S-35S
intlev = 37; intdep = 700; % 700m
dweit = depthData(2:intlev-1)-depthData(1:intlev-2);
levs = [1:intlev];
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
Tempint = Temp(:,:,intlev-1,:)+(Temp(:,:,intlev,:)-Temp(:,:,intlev-1,:))*(intdep-depthData(intlev-1))/(depthData(intlev)-depthData(intlev-1));
Tempz = cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:intlev-1,:).*(permute(dweit,[3 2 1])),Tempint*(intdep-depthData(intlev-1)));
VVint = VV(:,:,intlev-1,:)+(VV(:,:,intlev,:)-VV(:,:,intlev-1,:))*(intdep-depthData(intlev-1))/(depthData(intlev)-depthData(intlev-1));
VVz = cat(3,VV(:,:,1)*5,VV(:,:,2:intlev-1).*(permute(dweit,[3 2 1])),VVint*(intdep-depthData(intlev-1)));
Tempzy = Tempz*111*1000; % 111 km / latitude
VVzy = VVz*111*1000;
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(151:360,:,:,:),Tempzyx(1:150,:,:,:)); 
VVzyx(find(isnan(Tempzyx(:,:,:,1)) == 1)) = nan;
Vzyxsub_r_0 = cat(1,VVzyx(151:360,:,:),VVzyx(1:150,:,:)); 
spac_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:150,lats,levs,:),1),2),3));
spacV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:150,lats,levs),1),2),3));
spac0V = spac_0./spacV;
sia_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(151:360,lats,levs,:),1),2),3));
siaV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(151:360,lats,levs),1),2),3));
sia0V = sia_0./siaV;
sdiff = spac0V-sia0V;
sat_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(151:230,lats,levs,:),1),2),3));
satV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(151:230,lats,levs),1),2),3));
sat0V = sat_0./satV;
sio_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(231:360,lats,levs,:),1),2),3));
sioV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(231:360,lats,levs),1),2),3));
sio0V = sio_0./sioV;

zpac = cp*ro*squeeze(nansum(nansum(Tzyxsub_r_0(1:150,1:90,levs,:),1),3)); % zonal sum
pacS = cp*ro*squeeze(nansum(nansum(Vzyxsub_r_0(1:150,1:90,levs),1),3));
zpacS = zpac./pacS';
zatl = cp*ro*squeeze(nansum(nansum(Tzyxsub_r_0(151:230,1:90,levs,:),1),3)); % zonal sum
atlS = cp*ro*squeeze(nansum(nansum(Vzyxsub_r_0(151:230,1:90,levs),1),3));
zatlS = zatl./atlS';
zio = cp*ro*squeeze(nansum(nansum(Tzyxsub_r_0(231:360,1:90,levs,:),1),3)); % zonal sum
ioS = cp*ro*squeeze(nansum(nansum(Vzyxsub_r_0(231:360,1:90,levs),1),3));
zioS = zio./ioS';
%% OHC/m3  time series
sdiff = sdiff(21:101); % 1920-
sia0V = sia0V(21:101);
spac0V = spac0V(21:101);
close all;
Fig = figure('position',[700 100 800 400]);
p1 = bar([1:81],sdiff,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7],'BarWidth',1);
hold on
set(gca,'YLim',[-0.2,0.8],'YTick',[-0.2:0.2:0.8],'YColor','r');
ylabel('Difference  (W m^-^3)','FontSize',14);
yb = polyfit([1:81],sdiff,1);
plot([1:81],polyval(yb,[1:81]),'--','Color',[0.3,0.3,0.3],'linewidth',2)
text(55,0.2,[num2str(roundn(yb(1),-4)*10),' J/m^3 per decade'],'fontsize',14,'Color',[0.3,0.3,0.3])

yyaxis right
ymax = 4; ymin = 3; yint = 0.2;
p2 = plot(sia0V,'-','color',[.64 .08 .18],'LineWidth',3);
hold on
p3 = plot(spac0V,'-','color',[.04 .35 .56],'LineWidth',3);
plot(zeros(1,81),'k','LineWidth',1)
yb = polyfit([1:81],sia0V,1);
plot([1:81],polyval(yb,[1:81]),'--','Color',[.64 .08 .18],'linewidth',2)
text(55,3.5,[num2str(roundn(yb(1),-4)*10),' J/m^3 per decade'],'fontsize',14,'Color',[.64 .08 .18])
yb = polyfit([1:81],spac0V,1);
plot([1:81],polyval(yb,[1:81]),'--','Color',[.04 .35 .56],'linewidth',2)
text(55,3.8,[num2str(roundn(yb(1),-4)*10),' J/m^3 per decade'],'fontsize',14,'Color',[.04 .35 .56])
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel(' (J m^-^3)','FontSize',14);
legend([p3, p2, p1],'Pacific','Atlantic-Indian Ocean','Difference','Location','northwest')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/TimeSeries_OHC_Wm3_',depthstr,'_35S_90S.png'],'-dpng','-r300')
%% OHC J  time series
sdiff = (spac_0-sia_0)/10^15;
sdiff = sdiff(21:101); % 1920-
sia_0 = sia_0(21:101);
spac_0 = spac_0(21:101);
%%
close all;
Fig = figure('position',[700 100 800 400]);
p3 = plot(spac_0/10^15,'-','color',[.04 .35 .56],'LineWidth',3);
hold on
set(gca,'YLim',[195,245],'YTick',[195:10:245],'YColor','r');
ylabel('Pacific  (PW)','FontSize',14);
yb = polyfit([1:81],spac_0/10^15,1);
plot([1:81],polyval(yb,[1:81]),'--','Color',[.04 .35 .56],'linewidth',2)
text(55,221,[num2str(roundn(yb(1),-4)*10),' PJ per decade'],'fontsize',14,'Color',[.04 .35 .56])

yyaxis right
ymax = 400; ymin = 350; yint = 10;
p2 = plot(sia_0/10^15,'-','color',[.64 .08 .18],'LineWidth',3);
hold on
plot(zeros(1,81),'k','LineWidth',1)
yb = polyfit([1:81],sia_0/10^15,1);
plot([1:81],polyval(yb,[1:81]),'--','Color',[.64 .08 .18],'linewidth',2)
text(55,390,[num2str(roundn(yb(1),-4)*10),' PJ per decade'],'fontsize',14,'Color',[.64 .08 .18])
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel('Atlantic-Indian Ocean (PW)','FontSize',14);
legend([p3, p2],'Pacific','Atlantic-Indian Ocean','Location','northwest')
legend('boxoff')
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/TimeSeries_OHC_PW_',depthstr,'m_35S_90S.png'],'-dpng','-r300')
%% zonal OHC
lonstr = 150;
startyr = 1960;
endyr = 2020;
regstr1 = 'Pacific'; regstr2 = 'Pac'; clor = [.04 .35 .56];
var = zpac(:,startyr-1919:endyr-1919)';
% regstr1 = 'Atlantic'; regstr2 = 'Atl'; clor = [.64 .08 .18];
% var = zatl(:,startyr-1919:endyr-1919)';
% regstr1 = 'Indian Ocean'; regstr2 = 'IO'; clor = [.64 .08 .18];
% var = zio(:,startyr-1919:endyr-1919)';
x = [1:size(var,1)]';
clear trd 
for i = 1:size(var,2);
    par=polyfit(x,var(:,i),1); % regression parameters
    trd(i) = par(1);
end
cli = nanmean(var,1);
map1 = cli; map2 = trd; 
close all;
Fig = figure('position',[700 100 800 400]);
ymin = 0; ymax = .05; yint = .01;
p2 = plot(latData(20:60),map2(20:60)/10^21,'-','color',clor,'LineWidth',3);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);
hold on
set(gca,'XLim',[-70,-30],'xtick',[-70:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('Trend  (ZJ yr^-^1)','FontSize',14);
plot(-35*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
p1 = plot(latData(20:60),map1(20:60)/10^21,'-','color','k','LineWidth',3);
hold on
% p4 = plot(climap1,'-','color','k','LineWidth',3);
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax]*1000,'YTick',[ymin:yint:ymax]*1000,'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12),ylabel('Climo (ZJ)','FontSize',14);
legend([p1, p2],['Climo (',regstr1,')'],['Trend (',regstr1,')'],'Location','northwest')
legend('boxoff')

yyaxis left
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);
print(Fig,['D:/figures/CESM/Yearly/LE/0_LEM/20240425_SO_Warming/OHCzonal_',regstr2,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% zonal OHC /m2
lonstr = 150;
startyr = 1960;
endyr = 2020;
% regstr1 = 'Pacific'; regstr2 = 'Pac'; clor = [.04 .35 .56];
% var = zpacS(:,startyr-1939:endyr-1939)';
% regstr1 = 'Atlantic'; regstr2 = 'Atl'; clor = [.64 .08 .18];
% var = zatlS(:,startyr-1939:endyr-1939)';
regstr1 = 'Indian Ocean'; regstr2 = 'IO'; clor = [.64 .08 .18];
var = zioS(:,startyr-1939:endyr-1939)';

x = [1:size(var,1)]';
clear trd 
for i = 1:size(var,2);
    par=polyfit(x,var(:,i),1); % regression parameters
    trd(i) = par(1);
end
cli = nanmean(var,1);
map1 = cli; map2 = trd; 
close all;
Fig = figure('position',[700 100 800 400]);
ymin = 0; ymax = 16; yint = 4;
p2 = plot(latData(20:60),map2(20:60)*1000,'-','color',clor,'LineWidth',3);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);
hold on
set(gca,'XLim',[-70,-30],'xtick',[-70:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('Trend  (ZJ m^-^2 yr^-^1)','FontSize',14);
plot(-35*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
p1 = plot(latData(20:60),map1(20:60),'-','color','k','LineWidth',3);
hold on
% p4 = plot(climap1,'-','color','k','LineWidth',3);
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12),ylabel('Climo (ZJ m^-^2)','FontSize',14);
legend([p1, p2],['Climo (',regstr1,')'],['Trend (',regstr1,')'],'Location','northwest')
legend('boxoff')

yyaxis left
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);
print(Fig,['D:/figures/CESM/Yearly/LE/0_LEM/20240425_SO_Warming/OHCm2_zonal_',regstr2,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')

%% Tendency of OHC (minus 1960)
% sat_0 = sat_0(21:101);
% sio_0 = sio_0(21:101);
%%
close all;
Fig = figure('position',[700 100 800 400]);
ln1 = (sio_0-sio_0(1))/10^21;
p3 = plot(ln1,'-','color',[.04 .35 .56],'LineWidth',3);
hold on
set(gca,'YLim',[0,60],'YTick',[195:10:245]);
ylabel('Pacific  (PW)','FontSize',14);
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel('Atlantic-Indian Ocean (PW)','FontSize',14);
% legend([p3, p2],'Pacific','Atlantic-Indian Ocean','Location','northwest')
% legend('boxoff')
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/TimeSeries_OHC_PW_',depthstr,'m_35S_90S.png'],'-dpng','-r300')

%%
ncid=netcdf.open('/Volumes/CESM-post2/CESM-post/LE/Temperature/his_sst_ensmean.nc','NOWRITE');
ncdisp('/Volumes/CESM-post2/CESM-post/LE/Temperature/his_sst_ensmean.nc');
tsMMM = ncread('/Volumes/CESM-post2/CESM-post/LE/Temperature/his_sst_ensmean.nc','TEMP'); 
sstMMM = permute(tsMMM(:,:,1,:),[1 2 4 3]);
datadir='/Volumes/CESM-post2/CESM-post/LE/Temperature/his_sst_yr1x1/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); 
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
%
for s = 1
    s
    filename=[datadir,filelist(s).name];    
    ncid=netcdf.open(filename,'NC_NOWRITE');
    Tem = ncread(filename,'TEMP');
% sst minus MMM & anomaly (without filter)
    sst0 = permute(Tem(:,:,1,:),[1,2,4,3]);
    sst1 = sst0 - sstMMM; % minus MMM (trend)
    sstd(:,:,:,s) = sst1 - mean(sst1,3); % anomaly

    load(['MatFile/0_700m_South/Temsubadf',num2str(s),'.mat']);
    Temsubadf(:,:,:,1:4) = []; Temsubadf(:,:,:,end-3:end) = [];
% sst minus MMM & anomaly & (8-yr filter) 
    load(['MatFile/sst_Global/sstadf',num2str(s),'.mat']);
    sstadf(:,:,1:4) = []; sstadf(:,:,end-3:end) = [];
    sstdf(:,:,:,s) = sstadf; 
% 0-707m mean
    Tsubm(:,:,:,s) = permute(nansum(Temsubadf(:,:,1,:)*5+Temsubadf(:,:,2:37,:).*permute(dweit,[3 2 1 4]),3)/707,[1 2 4 3]);
% % meridional mean
%     names = '46S-61S';
%     Tem1m(:,:,:,s) = latmean(Temsubadf(:,:,1:37,:),[29:44],latData);
% end
% %
% for s = 1:k
%     load(['MatFile/nhflx_Global/nhflxadf',num2str(s),'.mat']);
%     nhflxadf(:,:,1:4) = []; nhflxadf(:,:,end-3:end) = [];
%     nhflxdf(:,:,:,s) = nhflxadf;
% end
% for s = 1:40
%     load(['MatFile/slp_Global/slpadf',num2str(s),'.mat']);
%     slpadf(:,:,1:4) = []; slpadf(:,:,end-3:end) = [];
%     slp(:,:,:,s) = double(slpadf);
end
nyr = size(Tsubm,3);
%% surface fields 
% CURLZ
datadirx1='/Volumes/CESM-post2/CESM-post/LE/TAUX/historical/taux_yr_1x1/'; %指定批量数据所在的文件夹
filelistx1=dir([datadirx1,'*.nc']); %指定批量数据的类型
ncdisp([datadirx1,filelistx1(1).name]);
datadirx2='/Volumes/CESM-post2/CESM-post/LE/TAUX/'; %指定批量数据所在的文件夹
filelistx2=dir([datadirx2,'*.nc']); %指定批量数据的类型

datadiry1='/Volumes/CESM-post2/CESM-post/LE/TAUY/historical/tauy_yr_1x1/'; %指定批量数据所在的文件夹
filelisty1=dir([datadiry1,'*.nc']); %指定批量数据的类型
ncdisp([datadiry1,filelisty1(1).name]);
datadiry2='/Volumes/CESM-post2/CESM-post/LE/TAUY/'; %指定批量数据所在的文件夹
filelisty2=dir([datadiry2,'*.nc']); %指定批量数据的类型

k=length(filelisty1);
clear Curlz
for s=1:k
    s 
    filenamex1=[datadirx1,filelistx1(s).name];
    ncid=netcdf.open(filenamex1,'NC_NOWRITE');
    Taux1 = ncread(filenamex1,'TAUX');
    filenamex2=[datadirx2,filelistx2(s).name];
    ncid=netcdf.open(filenamex2,'NC_NOWRITE');
    Taux2 = ncread(filenamex2,'TAUX');
    Taux = cat(3,Taux1(:,:,1:86),Taux2(:,:,1:75));

    filenamey1=[datadiry1,filelisty1(s).name];
    ncid=netcdf.open(filenamey1,'NC_NOWRITE');
    Tauy1 = ncread(filenamey1,'TAUY');
    filenamey2=[datadiry2,filelisty2(s).name];
    ncid=netcdf.open(filenamey2,'NC_NOWRITE');
    Tauy2 = ncread(filenamey2,'TAUY');
    Tauy = cat(3,Tauy1(:,:,1:86),Tauy2(:,:,1:75));

    for i = 1:size(Taux,3);
        Curlz(:,:,i,s) = ra_windstrcurl(latData,lonData,Taux(:,:,i)',Tauy(:,:,i)',1);  % curlz lat*lon
    end
end
[d1 d2 d3 d4] = size(Curlz);
curlz = permute(Curlz,[2 1 3 4]);
curlzMMM = nanmean(curlz,4);
% TAUX
load('MatFile/lonData.mat'); load('MatFile/latData.mat');
filename1 = '/Volumes/CESM-post2/CESM-post/LE/TAUX/taux_ensmean.nc';
filename2 = '/Volumes/CESM-post2/CESM-post/LE/TAUX/TAUXensmean.nc';
tauxMMM = loadsurface(filename1,filename2,'TAUX');
%% NHFLX
filename1 = '/Volumes/CESM-post2/CESM-post/LE/FLNS/FLNSensmean_1920-2005.nc';
filename2 = '/Volumes/CESM-post2/CESM-post/LE/FLNS/FLNSensmean_2006-2080.nc';
flnsMMM = loadsurface(filename1,filename2,'FLNS');
filename1 = '/Volumes/CESM-post2/CESM-post/LE/FSNS/FSNSensmean_1920-2005.nc';
filename2 = '/Volumes/CESM-post2/CESM-post/LE/FSNS/FSNSensmean_2006-2080.nc';
fsnsMMM = loadsurface(filename1,filename2,'FSNS');
filename1 = '/Volumes/CESM-post2/CESM-post/LE/LHFLX/LHFLXensmean_1920-2005.nc';
filename2 = '/Volumes/CESM-post2/CESM-post/LE/LHFLX/LHFLXensmean_2006-2080.nc';
lhflxMMM = loadsurface(filename1,filename2,'LHFLX');
filename1 = '/Volumes/CESM-post2/CESM-post/LE/SHFLX/SHFLXensmean_1920-2005.nc';
filename2 = '/Volumes/CESM-post2/CESM-post/LE/SHFLX/SHFLXensmean_2006-2080.nc';
shflxMMM = loadsurface(filename1,filename2,'SHFLX');
nhflxMMM = fsnsMMM-shflxMMM-lhflxMMM-flnsMMM;
%% climo
startyr = 1960; endyr = 2020;

units = 'W m^-^2'; name = 'NHFLX'; 
map = mean(mean(nhflxMMM(:,:,1960-1919:2020-1919),3),4);
ticks = 100;

% units = 'Pa'; name = 'TAUX'; 
% map = -mean(mean(tauxMMM(:,:,1960-1919:2020-1919),3),4);
% ticks = 0.5;

% units = 'Pa m^-^1'; name = 'CURLZ'; 
% map = -mean(mean(curlzMMM(:,:,1960-1919:2020-1919),3),4);
% ticks = 3*10^-7;

close all;
ftsz = 12; 
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(map,[-ticks*5:ticks/10:ticks*5],ticks,lonData,latData,12)
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],[units],5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/',name,'_climo_',num2str(startyr),'_',num2str(endyr),'_2.png'],'-dpng','-r300')
%% zonal plot
% net heat flux
lonstr = 150;
startyr = 1960;
endyr = 2020;
var_r = cat(1,nhflxMMM(1:lonstr,:,:),nhflxMMM(lonstr+1:end,:,:));
var1 = squeeze(nanmean(nhflxMMM(1:140,1:90,startyr-1919:endyr-1919),1))';
var2 = squeeze(nanmean(nhflxMMM(141:360,1:90,startyr-1919:endyr-1919),1))';
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
Fig = figure('position',[700 100 800 400]);
ymin1 = -.21; ymax1 = .28; yint1 = .07;
ymin2 = -33; ymax2 = 44; yint2 = 11;
p1 = plot(latData(20:60),trd1(20:60),'-','color',[.04 .35 .56],'LineWidth',3);
hold on
p2 = plot(latData(20:60),trd2(20:60),'-','color',[.64 .08 .18],'LineWidth',3);
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint1:ymax1],'YColor',clor);
set(gca,'XLim',[-70.5,-30.5],'xtick',[-70.5:10:-30.5],'fontsize',14,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('Trend  (W/yr)','FontSize',14);
plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
p3 = plot(latData(20:60),cli1(20:60),'--','color',[.04 .35 .56],'LineWidth',1.5);
hold on
p4 = plot(latData(20:60),cli2(20:60),'--','color',[.64 .08 .18],'LineWidth',1.5);
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:yint2:ymax2],'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',14),ylabel('Climo (W)','FontSize',14);
legend([p1, p3, p2, p4],['Pacific trend'],['Pacific climo'],['Atlantic-Indian Ocean trend'],...
    ['Atlantic-Indian Ocean climo'],'Location','northwest','numcolumns',2,'edgecolor','w')
% 手动保存
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_Warming/NHFLXzonal_',regstr2,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%%
startyr = 1960;  endyr = 2020;
var = permute(curlzMMM(:,:,startyr-1919:endyr-1919),[3 1 2]);
trd = trend_cal_3D(var);
max(trd,[],'all'), min(trd,[],'all')

units = '10^-^8 Pa m^-^1 decade^-^1'; name = 'CURLZ'; 
map = -trd*10*10^8; ticks = 1;

% units = 'K'; name = 'SST'; 
% map = trd*10; ticks = .2;
% 
% units = '10^-^3 Pa decade^-^1'; name = 'TAUX'; 
% map = -trd*10*10^3; ticks = 8;

% units = 'W m^-^2'; name = 'NHFLX'; 
% map = trd*10; ticks = 4;

close all;
ftsz = 12; 
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(map,[-ticks*5:ticks/10:ticks*5],ticks,lonData,latData,12)
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],[units,' decade^-^1'],5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots0(h0,[.4 .4 .4],lonDataE,latDataE);
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_warming/',name,'_trend_',num2str(startyr),'_',num2str(endyr),'_2.png'],'-dpng','-r300')
%% MOC
datadirx1='/Volumes/CESM-post2/CESM-post/LE/MOC/historical/'; %指定批量数据所在的文件夹
filelistx1=dir([datadirx1,'*.nc']); %指定批量数据的类型
ncdisp([datadirx1,filelistx1(1).name]);
datadirx2='/Volumes/CESM-post2/CESM-post/LE/MOC/'; %指定批量数据所在的文件夹
filelistx2=dir([datadirx2,'*2006*.nc']); %指定批量数据的类型
k=length(filelistx1);
clear Curlz
for s=1
    s 
    filenamex1=[datadirx1,filelistx1(s).name];
    ncid=netcdf.open(filenamex1,'NC_NOWRITE');
    MOC01 = ncread(filenamex1,'MOC');
    filenamex2=[datadirx2,filelistx2(s).name];
    ncid=netcdf.open(filenamex2,'NC_NOWRITE');
    MOC02 = ncread(filenamex2,'MOC');
    MOCm1 = cat(5,MOC01(:,:,:,:,111:156),MOC02(:,:,:,:,1:15)); % 1960-2020
    lat_aux_grid = ncread(filenamex1,'lat_aux_grid');
    moc_z = ncread(filenamex1,'moc_z');
    moc_comp = ncread(filenamex1,'moc_comp');
    transport_reg = ncread(filenamex1,'transport_reg');
end
for s=2:k
    s 
    filenamex1=[datadirx1,filelistx1(s).name];
    ncid=netcdf.open(filenamex1,'NC_NOWRITE');
    MOC01 = ncread(filenamex1,'MOC');
    filenamex2=[datadirx2,filelistx2(s).name];
    ncid=netcdf.open(filenamex2,'NC_NOWRITE');
    MOC02 = ncread(filenamex2,'MOC');
    MOCmo(:,:,:,:,:,s) = cat(5,MOC01(:,:,:,:,41:86),MOC02(:,:,:,:,1:15)); % 1960-2020
end
MOC = cat(6,MOCm1,MOCmo);
clear MOCm1 MOCmo
MOCcli = nanmean(nanmean(MOC,6),5);
%% plot
close all;
map = MOCcli(:,:,1,2);
max(map,[],'all')
min(map,[],'all')
Fig = figure('position',[10 50 600 500]);
ax = axes('Position',[0.1 0.09 0.8 0.89],'fontsize',12,'box','on');
contourf(lat_aux_grid,-moc_z,map',[-50:5:50])
caxis([-50,50]);
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(flipud(bl_re4));
set_colorbar([0.92 0.09 0.02 0.89],['Sv'],4.5,12,[-50:10:50],[-50:10:50])
set(gca,'YLim',[-550000,0],'YTick',[-500000:100000:0],'YTickLabel',[5000:-1000:0])
ylabel('Depth (m)','fontsize',12); xlabel('Latitude','fontsize',12)
set(gca,'XLim',[-90,90],'XTick',[-90:30:90],'XTickLabel',[-90:30:90])
print(Fig,['/Users/yysong/Desktop/figures/CESM/climo/MOC.Atlantic.EularianMean.climo.1960-2020.png'],'-dpng','-r300')









function varMMM = loadsurface(filename1,filename2,varstr)
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
varMMM1 = ncread(filename1,varstr);
varMMM1(:,:,87) = [];
lonData = ncread(filename1,'lon');
latData = ncread(filename1,'lat');

ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
varMMM2 = ncread(filename2,varstr);

varMMM = double(cat(3,varMMM1(:,:,1:86),varMMM2(:,:,1:75)));
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
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
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
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end
function testdots(h,clor,lonData,latData)
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',2,'color',clor); % F-test dots
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
function [ts] = latmean(var,lats,latData)
% average along longitude (all latitudes)
% var is lon*lat*depth*time
% ts is lon*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,2)./nansum(cos(latData(lats)'/180*pi).*var2,2),[1 3 4 2]);
end
