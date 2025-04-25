clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
load("MatFile/lonData.mat");
load("MatFile/latData.mat");
load("MatFile/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
%% VVEL
addpath /Volumes/Togo4T/1_matlab/help;
VVEL1 = ncread('/Volumes/CESM-post2/CESM-post/LE/VVEL/VVELensmean_1920-2005.nc','VVEL');
VVEL2 = ncread('/Volumes/CESM-post2/CESM-post/LE/VVEL/VVELensmean_2006-2080.nc','VVEL');
VVEL = cat(4,VVEL1(:,:,1:47,41:end),VVEL2(:,:,1:47,1:15)); % 1960-2020
VISOP1 = ncread('/Volumes/CESM-post2/CESM-post/LE/VISOP/VISOPensmean_1920-2005.nc','VISOP');
VISOP2 = ncread('/Volumes/CESM-post2/CESM-post/LE/VISOP/VISOPensmean_2006-2080.nc','VISOP');
VISOP = cat(4,VISOP1(:,:,1:47,41:end),VISOP2(:,:,1:47,1:15));
VSUBM1 = ncread('/Volumes/CESM-post2/CESM-post/LE/VSUBM/VSUBMensmean_1920-2005.nc','VSUBM');
VSUBM2 = ncread('/Volumes/CESM-post2/CESM-post/LE/VSUBM/VSUBMensmean_2006-2080.nc','VSUBM');
VSUBM = cat(4,VSUBM1(:,:,1:47,41:end),VSUBM2(:,:,1:47,1:15));
VRES = (VVEL + VISOP + VSUBM);
vcli = nanmean(VRES,4);
clear VVEL* VISOP* VSUBM*
%% WVEL
WVEL1 = ncread('/Volumes/CESM-post2/CESM-post/LE/WVEL/WVELensmean_1920-2005.nc','WVEL');
WVEL2 = ncread('/Volumes/CESM-post2/CESM-post/LE/WVEL/WVELensmean_2006-2080.nc','WVEL');
WVEL = cat(4,WVEL1(:,:,1:47,41:end),WVEL2(:,:,1:47,1:15)); % 1960-2020
WISOP1 = ncread('/Volumes/CESM-post2/CESM-post/LE/WISOP/WISOPensmean_1920-2005.nc','WISOP');
WISOP2 = ncread('/Volumes/CESM-post2/CESM-post/LE/WISOP/WISOPensmean_2006-2080.nc','WISOP');
WISOP = cat(4,WISOP1(:,:,1:47,41:end),WISOP2(:,:,1:47,1:15));
WSUBM1 = ncread('/Volumes/CESM-post2/CESM-post/LE/WSUBM/WSUBMensmean_1920-2005.nc','WSUBM');
WSUBM2 = ncread('/Volumes/CESM-post2/CESM-post/LE/WSUBM/WSUBMensmean_2006-2080.nc','WSUBM');
WSUBM = cat(4,WSUBM1(:,:,1:47,41:end),WSUBM2(:,:,1:47,1:15));
WRES = (WVEL + WISOP + WSUBM);
wcli = nanmean(WRES,4);
clear WVEL* WISOP* WSUBM*
%% TEMP
TEMP1 = ncread('/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_1920-2005.nc','TEMP');
TEMP2 = ncread('/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_2006-2022.nc','TEMP');
TEMP = cat(4,TEMP1(:,:,1:47,41:end),TEMP2(:,:,1:47,1:15)); % 1960-2020
Tcli = nanmean(TEMP,4);
clear TEMP*
%% plot
lats = 1:60;
VVELcli_r = cat(1,vcli(150:360,lats,1:37,:),vcli(1:149,lats,1:37,:));
WVELcli_r = cat(1,wcli(150:360,lats,1:37,:),wcli(1:149,lats,1:37,:));
TEMPcli_r = cat(1,Tcli(150:360,lats,1:37,:),Tcli(1:149,lats,1:37,:));
vp = squeeze(nanmean(VVELcli_r(1:140,:,:),1));
via = squeeze(nanmean(VVELcli_r(141:end,:,:),1));
vi = squeeze(nanmean(VVELcli_r(230:end,:,:),1));
va = squeeze(nanmean(VVELcli_r(141:229,:,:),1));
wp = squeeze(nanmean(WVELcli_r(1:140,:,:),1));
wia = squeeze(nanmean(WVELcli_r(141:end,:,:),1));
wi = squeeze(nanmean(WVELcli_r(230:end,:,:),1));
wa = squeeze(nanmean(WVELcli_r(141:229,:,:),1));
Tp = squeeze(nanmean(TEMPcli_r(1:140,:,:),1));
Tia = squeeze(nanmean(TEMPcli_r(141:end,:,:),1));
Ti = squeeze(nanmean(TEMPcli_r(230:end,:,:),1));
Ta = squeeze(nanmean(TEMPcli_r(141:229,:,:),1));

%% V cli
close all;
ftsz = 14;
ticks = 1.5;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:37),va',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([latData(35) latData(35)],[0 -700],'linewidth',1.5,'color','k')
line([latData(55) latData(55)],[0 -700],'linewidth',1.5,'color','k')
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(flipud(bl_re4));
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','cm s^-^1','Units','normalized','position',[4.5 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[-90,latData(60)]);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[-700:100:0],'FontSize',14);
ylabel('Depth (m)'); xlabel('Latitude')
print(Fig,['/Users/yysong/Desktop/figures/CESM/climo/Vcli_150E_Atlantic.png'],'-dpng','-r300');
%% W cli
close all;
ftsz = 14;
ticks = .0003;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:37),wi',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([latData(35) latData(35)],[0 -700],'linewidth',1.5,'color','k')
line([latData(55) latData(55)],[0 -700],'linewidth',1.5,'color','k')
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(flipud(bl_re4));
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','cm s^-^1','Units','normalized','position',[4.5 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[-90,latData(60)]);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[-700:100:0],'FontSize',14);
ylabel('Depth (m)'); xlabel('Latitude')
print(Fig,['/Users/yysong/Desktop/figures/CESM/climo/Wcli_150E_IndianOcean.png'],'-dpng','-r300');
%% T cli
close all;
ftsz = 14;
ticks = 8;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:37),Ti',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([latData(35) latData(35)],[0 -700],'linewidth',1.5,'color','k')
line([latData(55) latData(55)],[0 -700],'linewidth',1.5,'color','k')
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(flipud(bl_re4));
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','cm s^-^1','Units','normalized','position',[4.5 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[-90,latData(60)]);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[-700:100:0],'FontSize',14);
ylabel('Depth (m)'); xlabel('Latitude')
print(Fig,['/Users/yysong/Desktop/figures/CESM/climo/Tcli_150E_IndianOcean.png'],'-dpng','-r300');



