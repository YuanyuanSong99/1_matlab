clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
datadir='/Volumes/Togo4T/data/EN4/annual/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'EN4*.h5']); %指定批量数据的类型
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
Tempa = double(Temp - nanmean(Temp(:,:,:,32:61),4)); % remove climatology from 1981-2010 
%% Tem
sst_iap = permute(Tempa(:,:,1,:),[1 2 4 3]);
lats = 1:60;
%------------------------- 0-2000m -----------------------------------------
% depthstr = '0-2000';
% dweit = depthData(2:31)-depthData(1:30); % depth weight
% Tsubraw = permute(nansum(cat(3,Temp(:,lats,1,:),Temp(:,lats,2:31,:).*(permute(dweit,[3 2 1]))),3)/2000,[1 2 4 3]); 
%% meridional mean
lats = 29:49; % 55S-35S
latstr = '35S_55S';
[Tempa_mm] = latmean(Tempa,lats,latData);
% ts is lon*depth*time
startyr = 1960;
endyr = 2020;
yrstr = '1960-2020';
var = permute(Tempa_mm(:,:,startyr-1949:endyr-1949),[3 1 2]);
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
set(gca,'YLim',[-5000,0]);
set(gca,'YTick',[-2000:400:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-400:0],'FontSize',14);
ylabel('Depth (m)')
% print(Fig,['/Volumes/Togo4T/figures/EN4.2.1/Yearly/20240101_SO_warming/meridionalmean_',latstr,'_trend_',yrstr,'.png'],'-dpng','-r300')
%% 150 longitude
close all;
ftsz = 20;
ticks = 0.1;
map = trd*10;
map_r = cat(1,map(150:360,:),map(1:150,:));
h0_r = cat(1,h0(150:360,:),h0(1:150,:));
lonData_r = [lonData(150:360);lonData(1:150)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,map_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
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
set(ch.Label,'String','K decade^-^1','Units','normalized','position',[6 0.5 0],'Rotation',-90,'fontsize',ftsz);
[i,j] = find(h0_r == 0); % significance test
plot(lonData_r(i(1:3:end)),-depthData(j(1:3:end)),'.','MarkerSize',5,'color',[.4 .4 .4])
set(gca,'XLim',[150,360+150]);
set(gca,'XTick',[150:60:360+150]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',ftsz);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:500:0],'FontSize',ftsz);
set(gca,'YTickLabel',[2000:-500:0],'FontSize',ftsz);
ylabel('Depth (m)')
print(Fig,['/Users/yysong/Desktop/figures/EN4.2.1/Yearly/20240705_SO_warming/meridionalmean_',latstr,'_trend_',yrstr,'_150E.png'],'-dpng','-r300')
%% original Tem time series 
lats = [1:55]; % 90S-35S
Tzyxsub_r_0 = cat(1,Tsubraw(180:360,:,:),Tsubraw(1:179,:,:)); 
lonDara_r = cat(1,lonData(180:360),lonData(1:179))
[spacz_0 spac_0] = areamean(Tzyxsub_r_0,1:121,lats,latData); 
[siaz_0 sia_0] = areamean(Tzyxsub_r_0,122:360,lats,latData); 
sdiff = spac_0-sia_0;
close all;
Fig = figure('position',[700 100 800 400]);
ymax = 0.9; ymin = 0.1; yint = 0.2;
p1 = bar([1:71],sdiff,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7],'BarWidth',1);
hold on
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','r');
ylabel('Difference  (K)','FontSize',14);
yb = polyfit([1:71],sdiff,1);
plot([1:71],polyval(yb,[1:71]),'--','Color',[0.3,0.3,0.3],'linewidth',2)
text(58,0.3,[num2str(roundn(yb(1),-4)*10),' K per decade'],'fontsize',14,'Color',[0.3,0.3,0.3])

yyaxis right
ymax = 4.7; ymin = 3.9; yint = 0.2;
p2 = plot(sia_0,'-','color',[.64 .08 .18],'LineWidth',3);
hold on
p3 = plot(spac_0,'-','color',[.04 .35 .56],'LineWidth',3);
plot(zeros(1,71),'k','LineWidth',1)
yb = polyfit([1:71],sia_0,1);
plot([1:71],polyval(yb,[1:71]),'--','Color',[.64 .08 .18],'linewidth',2)
text(58,4.4,[num2str(roundn(yb(1),-4)*10),' K per decade'],'fontsize',14,'Color',[.64 .08 .18])
yb = polyfit([1:71],spac_0,1);
plot([1:71],polyval(yb,[1:71]),'--','Color',[.04 .35 .56],'linewidth',2)
text(58,4.6,[num2str(roundn(yb(1),-4)*10),' K per decade'],'fontsize',14,'Color',[.04 .35 .56])
set(gca,'XLim',[1,71]);
set(gca,'XTick',[1:10:71]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel('Original (K)','FontSize',14);
legend([p2, p3, p1],'Atlantic-Indian Ocean','Pacific','Difference','Location','northwest')
legend('boxoff')
% print(Fig,['/Volumes/Togo4T/figures/IAP/Yearly/20231104_SO_warming/TimeSeries_Tem_original_',depthstr,'m_35S_90S.png'],'-dpng','-r300')
%% OHC
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
%------------------------- 0-2000m -----------------------------------------
% depthstr = '0-2000m';
% dweit = depthData(2:30)-depthData(1:29); % depth weight
% OHCsub = cp*ro*permute(nansum(cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:31,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]); % J/m2
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthData(2:24)-depthData(1:23); % depth weight
Temp700 = Temp(:,:,24,:)+(Temp(:,:,25,:)-Temp(:,:,24,:))*(700-depthData(24))/(depthData(25)-depthData(24));
OHCsub = cp*ro*permute(nansum(cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:24,:).*(permute(dweit,[3 2 1])),Temp700*(700-depthData(24))),3),[1 2 4 3]); % J/m2

OHCw = OHCsub/365/24/60/60; % W/m2
% Delta OHC
yrspan = 15;
yrs1 = 1960; 
yrs2 = 2005; 

map = nanmean(OHCw(:,:,yrs2-1949:yrs2-1949+yrspan-1),3)-nanmean(OHCw(:,:,yrs1-1949:yrs1-1949+yrspan-1),3);
max(map,[],'all')
min(map,[],'all')
% SO asymmetry
close all;
ftsz = 12; ticks = 60;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS((map),[-ticks*5:ticks/5:ticks*5],ticks,lonData,latData,12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
% m_line([0:1:360],-50,'linewidth',2,'color','k'); 
m_line(-60,[-90:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-90:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'/DeltaOHC  (W m^-^2)',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots(h0,[.4 .4 .4],lonData,latData);
% print(Fig,['/Volumes/Togo4Tfigures/IAP/Yearly/20231107_SO_warming/Delta_OHC_W_raw',depthstr,'_',num2str(yrs1),'_',num2str(yrs2),'_',num2str(yrspan),'yrs.png'],'-dpng','-r300')
%% OHC trend
startyr = 1960;
endyr = 2020;
var = permute(OHCw(:,:,startyr-1949:endyr-1949),[3 1 2]);
x = [1:size(var,1)]';
clear trd 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.01); % t test trend
h0(find(isnan(Temp(:,:,1,1)) == 1)) = nan;
max(trd,[],'all')
min(trd,[],'all')
%% SO plot
close all;
ftsz = 20; ticks = 1.5;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd,[-ticks*5:ticks/10:ticks*5],ticks,lonData,latData,ftsz)
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2',5.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots0(h0(1:2:end,1:2:end),[.4 .4 .4],lonData(1:2:end,1:2:end),latData(1:2:end,1:2:end));
print(Fig,['/Users/yysong/Desktop/figures/EN4.2.1/Yearly/20240622_SO_warming/OHC',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.3.png'],'-dpng','-r300')
%% Global plot
close all;
ftsz = 12; ticks = 15;
Fig = figure('position',[10 50 850 400]);
ax = axes('Position',[0.05 0.07 0.85 0.87],'fontsize',ftsz,'box','on');
contourVARra(trd*10,[-ticks*5:ticks/5:ticks*5],ticks,0,360,-90,90,lonData,latData)
title([depthstr,' OHC trend (W m^-^2 decade^-^1)'])
set_colorbar([0.9 0.065 0.02 0.875],'W m^-^2 decade^-^1',4,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots0(h0,[.4 .4 .4],lonData,latData);
print(Fig,['/Volumes/Togo4T/figures/EN4.2.1/Yearly/20240101_SO_warming/OHC',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.2.png'],'-dpng','-r300')
%% original OHC time series
lats = [1:49]; % 90S-35S
levs = [1:31]; % 2000 m
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
Temp = double(Temp);
Temp2000 = Temp(:,:,30,:)+(Temp(:,:,31,:)-Temp(:,:,30,:))*(2000-depthData(30))/(depthData(31)-depthData(30));
Tempz = cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:30,:).*(permute(dweit,[3 2 1])),Temp2000*(2000-depthData(30)));
VV2000 = VV(:,:,30,:)+(VV(:,:,31,:)-VV(:,:,30,:))*(2000-depthData(30))/(depthData(31)-depthData(30));
VVz = cat(3,VV(:,:,1)*5,VV(:,:,2:30).*(permute(dweit,[3 2 1])),VV2000*(2000-depthData(30)));
Tempzy = Tempz*111*1000; % 111 km / latitude
VVzy = VVz*111*1000;
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(180:360,:,:,:),Tempzyx(1:179,:,:,:)); 
VVzyx(find(isnan(Tempzyx(:,:,:,1)) == 1)) = nan;
Vzyxsub_r_0 = cat(1,VVzyx(180:360,:,:),VVzyx(1:179,:,:)); 
spac_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:121,lats,levs,:),1),2),3));
spacV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:121,lats,levs),1),2),3));
spac0V = spac_0./spacV;
sia_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(122:360,lats,levs,:),1),2),3));
siaV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(122:360,lats,levs),1),2),3));
sia0V = sia_0./siaV;
sdiff = spac0V-sia0V;
max(spac0V,[],'all')
%% OHC/m3  time series 
close all;
Fig = figure('position',[700 100 800 400]);
p1 = bar([1:71],sdiff,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7],'BarWidth',1);
hold on
set(gca,'YLim',[-0,1],'YTick',[-0:0.2:1],'YColor','r');
ylabel('Difference  (W m^-^3)','FontSize',14);
yb = polyfit([1:71],sdiff,1);
plot([1:71],polyval(yb,[1:71]),'--','Color',[0.3,0.3,0.3],'linewidth',2)
text(55,0.2,[num2str(roundn(yb(1),-4)*10),' W/m^3 per decade'],'fontsize',14,'Color',[0.3,0.3,0.3])

yyaxis right
ymax = 4; ymin = 3; yint = 0.2;
p2 = plot(sia0V,'-','color',[.64 .08 .18],'LineWidth',3);
hold on
p3 = plot(spac0V,'-','color',[.04 .35 .56],'LineWidth',3);
plot(zeros(1,81),'k','LineWidth',1)
yb = polyfit([1:71],sia0V,1);
plot([1:71],polyval(yb,[1:71]),'--','Color',[.64 .08 .18],'linewidth',2)
text(55,3.5,[num2str(roundn(yb(1),-4)*10),' W/m^3 per decade'],'fontsize',14,'Color',[.64 .08 .18])
yb = polyfit([1:71],spac0V,1);
plot([1:71],polyval(yb,[1:71]),'--','Color',[.04 .35 .56],'linewidth',2)
text(55,3.8,[num2str(roundn(yb(1),-4)*10),' W/m^3 per decade'],'fontsize',14,'Color',[.04 .35 .56])
set(gca,'XLim',[0,71]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1950:10:2020],'FontSize',14);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel(' (W m^-^3)','FontSize',14);
legend([p3, p2, p1],'Pacific','Atlantic-Indian Ocean','Difference','Location','northwest')
legend('boxoff')
print(Fig,['/Volumes/Togo4T/figures/EN4.2.1/Yearly/20240101_SO_warming/TimeSeries_OHC_Wm3_',depthstr,'_35S_90S.png'],'-dpng','-r300')
%% OHC W  time series
sdiff = (spac_0-sia_0)/10^15;
close all;
Fig = figure('position',[700 100 800 400]);
p3 = plot(spac_0/10^15,'-','color',[.04 .35 .56],'LineWidth',3);
hold on
set(gca,'YLim',[215,265],'YTick',[215:10:265],'YColor','r');
ylabel('Pacific  (PW)','FontSize',14);
yb = polyfit([1:71],spac_0/10^15,1);
plot([1:71],polyval(yb,[1:71]),'--','Color',[.04 .35 .56],'linewidth',2)
text(50,218,[num2str(roundn(yb(1),-4)*10),' PW per decade'],'fontsize',14,'Color',[.04 .35 .56])

yyaxis right
ymax = 420; ymin = 370; yint = 10;
p2 = plot(sia_0/10^15,'-','color',[.64 .08 .18],'LineWidth',3);
hold on
plot(zeros(1,81),'k','LineWidth',1)
yb = polyfit([1:71],sia_0/10^15,1);
plot([1:71],polyval(yb,[1:71]),'--','Color',[.64 .08 .18],'linewidth',2)
text(50,390,[num2str(roundn(yb(1),-4)*10),' PW per decade'],'fontsize',14,'Color',[.64 .08 .18])
set(gca,'XLim',[1,71]);
set(gca,'XTick',[1:10:71]);
set(gca,'XTickLabel',[1950:10:2020],'FontSize',14);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',[.64 .08 .18]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel('Atlantic-Indian Ocean (PW)','FontSize',14);
legend([p3, p2],'Pacific','Atlantic-Indian Ocean','Location','northwest')
legend('boxoff')
print(Fig,['/Volumes/Togo4T/figures/EN4.2.1/Yearly/20240101_SO_warming/TimeSeries_OHC_PW_',depthstr,'m_35S_90S.png'],'-dpng','-r300')








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
% filename = '/Volumes/Togo4T/data/NOAA/20CRv3/20thC_ReanV3/Monthlies/miscSI-MO/prmsl.mon.mean.nc'; % 1836-2015
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
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function contourfSPolar(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',80,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
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
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
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
saveas(Fig,['/Volumes/Togo4T/figures/IAP/Yearly/20230606_IPO_SouthernOcean_46_61S/',pngname,'.png'])
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
    m_plot(lonData(dots),latData(ceil(find(h == 0)/length(lonData))),'.','markersize',5,'color',clor); % F-test dots
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
        elseif abs(T) < tinv(1-alp/2,d1-2);
            h0(i,j) = 0;
        end
    end
end
end
