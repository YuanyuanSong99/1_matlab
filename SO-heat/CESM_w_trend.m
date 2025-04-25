clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
load("/Volumes/Togo4T/1_matlab/MatData/CESM/lonData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/latData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/depthData.mat");
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
clear WVEL* WISOP* WSUBM*
%% meridional mean
lats = 1:70; levs = 1:47;
var0 = cat(1,WTT(150:360,lats,levs,:),WTT(1:149,lats,levs,:));
varpac = squeeze(nanmean(var0(1:140,:,:,:),1));
varia = squeeze(nanmean(var0(141:end,:,:,:),1));
startyr = 1960;
endyr = 2020;
yrstr = '1960-2020';
var1 = permute(varpac(:,:,startyr-1959:endyr-1959),[3 1 2]);
var2 = permute(varia(:,:,startyr-1959:endyr-1959),[3 1 2]);
x = [1:size(var,1)]';
clear trd
for i = 1:size(var1,2);
    for j = 1:size(var1,3);
        par=polyfit(x,var1(:,i,j),1); % regression parameters
        trd1(i,j) = par(1); 
        par=polyfit(x,var2(:,i,j),1); % regression parameters
        trd2(i,j) = par(1); 
    end
end
h0 = trendtest(var1,0.05); % t test trend
max(trd1,[],'all')
min(trd1,[],'all')
%%
close all;
ftsz = 14;
ticks = 1*10^-6;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(levs),trd1',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
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
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:250:0],'FontSize',14);
set(gca,'YTickLabel',[-700:100:0],'FontSize',14);
ylabel('Depth (m)'); xlabel('Latitude')

Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(levs),trd2',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
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
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:250:0],'FontSize',14);
set(gca,'YTickLabel',[-700:100:0],'FontSize',14);
ylabel('Depth (m)'); xlabel('Latitude')
%% WNT
WTT1 = ncread('/Volumes/CESM-post2/WTTensmean_1920-2005.nc','WTT'); % degC/s
WTT2 = ncread('/Volumes/CESM-post2/WTTensmean_2006-2080.nc','WTT'); % degC/s
WTTmon = cat(4,WTT1(:,:,:,480+1:end),WTT2(:,:,:,1:180));
clear WTT1 WTT2
[d1 d2 d3 d4] = size(WTTmon);
WTT = squeeze(mean(reshape(WTTmon,d1,d2,d3,12,d4/12),4)); % yearly
clear WTTmon
lonData0 = ncread('/Volumes/CESM-post2/WTTensmean_2006-2080.nc','TLONG'); % degC/s
latData0 = ncread('/Volumes/CESM-post2/WTTensmean_2006-2080.nc','TLAT'); % degC/s
%% meridional mean
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
lats = 1:115; latstr = '35S-55S';
levs = 1:47;
dweit = [diff(depthData(levs));depthData(47)-2000];
var0 = ro*cp*permute(dweit,[3 2 1]).*cat(1,WTT(171:end,lats,levs,:),WTT(1:170,lats,levs,:));
varpac = squeeze(nanmean(var0(1:124,:,:,:),1));
varia = squeeze(nanmean(var0(125:end,:,:,:),1));
startyr = 1960;
endyr = 2020;
yrstr = '1960-2020';
var1 = permute(varpac(:,:,startyr-1959:endyr-1959),[3 1 2]);
var2 = permute(varia(:,:,startyr-1959:endyr-1959),[3 1 2]);
x = [1:size(var1,1)]';
clear trd
for i = 1:size(var1,2);
    for j = 1:size(var1,3);
        par=polyfit(x,var1(:,i,j),1); % regression parameters
        trd1(i,j) = par(1); 
        par=polyfit(x,var2(:,i,j),1); % regression parameters
        trd2(i,j) = par(1); 
    end
end
h1 = trendtest(var1,0.05); % t test trend
h2 = trendtest(var2,0.05); % t test trend
max(trd1,[],'all')
min(trd1,[],'all')
%%
close all;
ftsz = 20;
ticks = 1;
map = trd1*10;
h0 = h1;
Fig = figure('position',[100 100 800 400]);
xx = latData0(65,[18:94]);
contourf(xx,-depthData(levs),map(18:94,:)',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([-35 -35],[0 -2000],'linewidth',1.5,'color','k')
line([-55 -55],[0 -2000],'linewidth',1.5,'color','k')
line([xx(1)-1 xx(end)],[-700 -700],'linewidth',1.5,'color','k','linestyle','--')
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(bl_re4);
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
set(ch.Label,'String','W m^-^2 decade^-^1','Units','normalized','position',[5.5 0.5 0],'Rotation',-90,'fontsize',ftsz);
[i,j] = find(h0 == 0); % significance test
plot(latData0(65,i(1:2:end)),-depthData(j(1:2:end)),'.','MarkerSize',5,'color',[.4 .4 .4])
set(gca,'XLim',[-70,-30]);
set(gca,'XTick',[-70:10:-30],'XTickLabel',{'70^oS','60^oS', '50^oS', '40^oS', '30^oS'});
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:500:0],'FontSize',ftsz);
set(gca,'YTickLabel',[2000:-500:0],'FontSize',ftsz);
ylabel('Depth (m)'); xlabel('Latitude')
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240930_SO_warming/WTTroCp_Pac_meridionalmean_',latstr,'_trend_',yrstr,'.png'],'-dpng','-r300')

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
function [ts ts_zs] = areamean(var,lons,lats,latData);
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
function [ts] = latmean(var,lats,latData)
% average along longitude (all latitudes)
% var is lon*lat*depth*time
% ts is lon*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,2)./nansum(cos(latData(lats)'/180*pi).*var2,2),[1 3 4 2]);
end
