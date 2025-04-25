clc,clear,close all;
addpath(genpath('D:\1_matlab\help'));
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
% MOHT calculate 
inlev = 37; depthstr = '0-700m';
splon = 160; lonstr = '160E';
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:inlev)-depthData(1:inlev-1);
clear levdist
levdist(1) = 4; levdist(2:inlev) = dweit; % z distance
%% historical 
datadir1=['E:\CESM-post\LE\TEMP\yearly\']; % 指定批量数据所在的文件夹
filelist1=dir([datadir1,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='E:\CESM-post\LE\VVEL\yearly\'; % 指定批量数据所在的文件夹
filelist2=dir([datadir2,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
datadir3='E:\CESM-post\LE\VISOP\yearly\'; % 指定批量数据所在的文件夹
filelist3=dir([datadir3,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir3,filelist3(1).name]);
datadir4='E:\CESM-post\LE\VSUBM\yearly\'; % 指定批量数据所在的文件夹
filelist4=dir([datadir4,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir4,filelist4(1).name]);
k=length(filelist1)
clear OHCz* MOHT*1
lonData_r = cat(1,lonData(splon:360),lonData(1:splon-1));
%%
for s=1:40
    s 
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    VVEL = ncread(filename2,'VVEL')/100; % cm/s -> m/s
    filename3=[datadir3,filelist3(s).name];
    ncid=netcdf.open(filename3,'NC_NOWRITE');
    VISOP = ncread(filename3,'VISOP')/100; % cm/s -> m/s
    filename4=[datadir4,filelist4(s).name];
    ncid=netcdf.open(filename4,'NC_NOWRITE');
    VSUBM = ncread(filename4,'VSUBM')/100; % cm/s -> m/s
    % 0-700m MOHT
    Vvel = VVEL+VISOP+VSUBM;
    clear VVEL VISOP VSUBM
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    TemVel = Temp(:,1:90,1:inlev,1:86).*Vvel(:,1:90,1:inlev,1:86); % 1940-2005
    clear Vvel Temp
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    clear TemVel
    OHCzy = OHCz.*londist';
    clear OHCz
    OHCzy_r = cat(1,OHCzy(splon:360,:,:,:),OHCzy(1:splon-1,:,:,:));
    clear OHCzy
    MOHTp1(:,:,s) = permute(nansum(OHCzy_r(1:300-splon,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia1(:,:,s) = permute(nansum(OHCzy_r(301-splon:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTi1(:,:,s) = permute(nansum(OHCzy_r(382-splon:end,:,:,:),1),[2 3 4 1]); % io 20E-150E unit:W
    MOHTa1(:,:,s) = permute(nansum(OHCzy_r(301-splon:381-splon,:,:,:),1),[2 3 4 1]); % Atl 60W-20E unit:W
    MOHTall1(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
    
end
MOHTp1LE = nanmean(MOHTp1,3);
MOHTia1LE = nanmean(MOHTia1,3);
MOHTi1LE = nanmean(MOHTi1,3);
MOHTa1LE = nanmean(MOHTa1,3);
MOHTall1LE = nanmean(MOHTall1,3);
save("MatFile\MOHT_160E_700m\MOHTall1",'MOHTall1');
save("MatFile\MOHT_160E_700m\MOHTa1",'MOHTa1');
save("MatFile\MOHT_160E_700m\MOHTia1",'MOHTia1');
save("MatFile\MOHT_160E_700m\MOHTp1",'MOHTp1');

%%
close all
plot(mean(MOHTall1LE,2))
%% MOHT anomaly,detrend, and 8-yr filter
load("MatFile\MOHT_160E_700m\MOHTall1",'MOHTall1');
load("MatFile\MOHT_160E_700m\MOHTa1",'MOHTa1');
load("MatFile\MOHT_160E_700m\MOHTia1",'MOHTia1');
load("MatFile\MOHT_160E_700m\MOHTp1",'MOHTp1');
clear MOHTp1a MOHTia1a MOHTp1ad MOHTp1adf MOHTia1ad MOHTia1adf
for s = 1:40
MOHTp1a(:,:,s) = MOHTp1(:,:,s)-mean(MOHTp1(:,:,s),2);
MOHTia1a(:,:,s) = MOHTia1(:,:,s)-mean(MOHTia1(:,:,s),2);
dT = 1; cf = 1/8;
for i = 1:90
    MOHTp1ad(:,i,s) = detrend(MOHTp1a(i,:,s)'); % linear detrend
    MOHTp1adf(:,i,s) = lanczosfilter(MOHTp1ad(:,i,s),dT,cf,[],'low'); % 8 year filtered
    MOHTia1ad(:,i,s) = detrend(MOHTia1a(i,:,s)'); % linear detrend
    MOHTia1adf(:,i,s) = lanczosfilter(MOHTia1ad(:,i,s),dT,cf,[],'low'); % 8 year filtered
end
end
MOHTp1adf(1:4,:,:) = []; MOHTia1adf(1:4,:,:) = [];
MOHTp1adf(end-3:end,:,:) = []; MOHTia1adf(end-3:end,:,:) = [];
%%
clear parR7 h077 par87 h087 
for s = 1:40
index = TPIfz(42:79,s);
for i = 1:size(MOHTp1adf,2);
    [parR7(i,s),h077(i,s),t] = reg1_ttest(index,MOHTp1adf(:,i,s),0.1,1);
    [parR8(i,s),h088(i,s),t] = reg1_ttest(index,MOHTia1adf(:,i,s),0.1,1);
end
end
%%
close all;
map77 = mean(parR7,2)/10^15; % residual
map88 = mean(parR8,2)/10^15;
map7 = mean(par77,2)/10^15; % VVEL
map8 = mean(par88,2)/10^15;
mapE7 = -mean(par7,2)/10^15; % Ekman
mapE8 = -mean(par8,2)/10^15;
Fig = figure('position',[100 100 800 400]);
p1 = plot(latData(1:90),map77,'--','color',[0.90,0.51,0.31],'linewidth',1.5)
hold on
p2 = plot(latData(1:90),map88,'--','color',[0.17,0.73,0.69],'linewidth',1.5)
p3 = plot(latData(1:90),map7,'color',[0.90,0.51,0.31],'linewidth',1.5)
p4 = plot(latData(1:90),map8,'color',[0.17,0.73,0.69],'linewidth',1.5)
p5 = plot(latData(1:90),mapE7,'r-','linewidth',1.5)
p6 = plot(latData(1:90),mapE8,'b-','linewidth',1.5)
set(gca,'XLim',[-75,-35]);set(gca,'XTick',[-70:10:-35]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',[-.024,.032],'YTick',[-.04:.01:.04]); set(gca,'YTick',[-.04:.01:.04],'FontSize',14);
ylabel('MOHT (PW)')
legend(gca,[p4,p2,p6],'Atlantic & Indian Ocean Total Eulerian mean','Atlantic & Indian Ocean Total Residual','Atlantic & Indian Ocean Ekman',...
    'fontsize',12,'Location','north','Orientation','vertical','position',[0.28,0.63,0.18,0.3])
legend('boxoff')
ah = axes('position',get(gca,'position'),'Visible','off')
legend(ah,[p5,p1,p3],'Pacific Ekman','Pacific Total Residual','Pacific Total Eulerian mean',...
    'fontsize',12,'Location','north','Orientation','vertical','position',[0.2,0.1,0.18,0.3])
legend('boxoff')
print(Fig,['D:\figures\CESM\Yearly\LE\20240130_IPO_SO\MOHT_Eulerian_Residual_Ekman.png'],'-dpng','-r300')
%% dMOHT/dy/dlon -- VVEL
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dlonp = 300-160; dlonia = 360-dlonp;
clear dMdyp dMdyia
[d1] = length(par77);
dy = sw_dist([-90,-89],[0,0],'km')*10^3; % m
dMdyp(1,:) = nan; dMdyp(d1,:) = nan;
dMdyia(1,:) = nan; dMdyia(d1,:) = nan;
for i = 2:d1-1
    dMOHTp(i) = par77(i+1)-par77(i-1);
    dMdyp(i) = dMOHTp(i)/(2*dy);
    dMOHTia(i) = par88(i+1)-par88(i-1);
    dMdyia(i) = dMOHTia(i)/(2*dy);

end
%%
close all;
Fig = figure('position',[100 100 800 400])
map1 = dMdyp;
nump = find(map1 > 0); numn = find(map1 < 0);
map1 = dMdyp;
nump = find(map1 > 0); numn = find(map1 < 0);
h3 = bar(latData(nump)+0.5,map1(nump),'FaceColor',[.64 .08 .18],'EdgeColor',[.64 .08 .18],'barwidth',.5)
hold on
% set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:ymax1/3:ymax1]);
bar(latData(numn)+0.5,map1(numn),'FaceColor',[.64 .08 .18],'EdgeColor',[.64 .08 .18],'barwidth',.5)
% ylabel('Trend (10^4 W m^-^1 Longitude^-^1 decade^-^1)')
% plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
% plot(-50*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
set(gca,'XLim',[-75,-35]);set(gca,'XTick',[-70:10:-35]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',[-3,5]*10^7,'YTick',[-3:1:5]*10^7); set(gca,'YTick',[-3:1:5]*10^7,'FontSize',14);
xlabel('Latitude','FontSize',14),ylabel('Climo (10^6 W m^-^1 Longitude^-^1)','FontSize',14);

map1 = dMdyia;
nump = find(map1 > 0); numn = find(map1 < 0);
h4 = bar(latData(nump),map1(nump),'FaceColor',[.04 .35 .56],'EdgeColor',[.04 .35 .56],'barwidth',.5)
hold on
% set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:ymax1/3:ymax1]);
bar(latData(numn),map1(numn),'FaceColor',[.04 .35 .56],'EdgeColor',[.04 .35 .56],'barwidth',.5)
% ylabel('Trend (10^4 W m^-^1 Longitude^-^1 decade^-^1)')
% plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
% plot(-50*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
set(gca,'XLim',[-75,-35]);set(gca,'XTick',[-70:10:-35]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',[-3,5]*10^7,'YTick',[-3:1:5]*10^7); set(gca,'YTick',[-3:1:5]*10^7,'FontSize',14);
ylabel('dMOHT/dy (W m^-^1)')
legend(gca,[h3,h4],'Pacific','Atlantic-Indian Ocean','fontsize',14,'Location','north','Orientation','vertical','position',[0.2,0.82,0.18,0.008])
legend('boxoff')
print(Fig,['D:\figures\CESM\Yearly\LE\20240130_IPO_SO\dMOHTdy.png'],'-dpng','-r300')
%% dutdx+dvtdy
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
lonData_r = cat(1,lonData(splon:360),lonData(1:splon-1));
dy = sw_dist([-90,-89],[0,0],'km')*10^3; % m
% ensmean
VVELMMM = ncread('E:\CESM-post\LE\VVEL\VVELensmean_1920-2005.nc','VVEL');
VISOPMMM = ncread('E:\CESM-post\LE\VISOP\VISOPensmean_1920-2005.nc','VISOP');
VSUBMMMM = ncread('E:\CESM-post\LE\VSUBM\VSUBMensmean_1920-2005.nc','VSUBM');
VMMM = VVELMMM + VISOPMMM + VSUBMMMM;
clear VVELMMM VISOPMMM VSUBMMMM
UVELMMM = ncread('E:\CESM-post\LE\UVEL\UVELensmean_1920-2005.nc','UVEL');
UISOPMMM = ncread('E:\CESM-post\LE\UISOP\UISOPensmean_1920-2005.nc','UISOP');
USUBMMMM = ncread('E:\CESM-post\LE\USUBM\USUBMensmean_1920-2005.nc','USUBM');
UMMM = UVELMMM + UISOPMMM + USUBMMMM;
clear UVELMMM UISOPMMM USBUMMMM
TemMMM = ncread('E:\CESM-post\LE\TEMP\TEMP.ensmean.1920-2005.hist.nc','TEMP');
VTMMM = VMMM(:,1:90,1:inlev,1:86).*TemMMM(:,1:90,1:inlev,1:86);
UTMMM = UMMM(:,1:90,1:inlev,1:86).*TemMMM(:,1:90,1:inlev,1:86);
clear VMMM UMMM TemMMM
dUTM(1,:,:,:) = UTMMM(1,:,:,:)-UTMMM(end,:,:,:);
dUTM(2:360,:,:,:) = UTMMM(2:end,:,:,:)-UTMMM(1:end-1,:,:,:);
dUTMdx = dUTM./londist';
clear dUTM UTMMM
dVTM(:,1,:,:) = zeros(360,1,inlev,86);
dVTM(:,2:90,:,:) = VTMMM(:,2:90,:,:)-VTMMM(:,1:89,:,:);
dVTMdy = dVTM./dy;
clear dVTM VTMMM
UVTMall = dUTMdx+dVTMdy;

%% historical 
datadir1=['E:\CESM-post\LE\TEMP\yearly\']; % 指定批量数据所在的文件夹
filelist1=dir([datadir1,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='E:\CESM-post\LE\VVEL\yearly\'; % 指定批量数据所在的文件夹
filelist2=dir([datadir2,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
datadir3='E:\CESM-post\LE\VISOP\yearly\'; % 指定批量数据所在的文件夹
filelist3=dir([datadir3,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir3,filelist3(1).name]);
datadir4='E:\CESM-post\LE\VSUBM\yearly\'; % 指定批量数据所在的文件夹
filelist4=dir([datadir4,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir4,filelist4(1).name]);
datadir5='E:\CESM-post\LE\UVEL\yearly\'; % 指定批量数据所在的文件夹
filelist5=dir([datadir5,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir5,filelist5(1).name]);
datadir6='E:\CESM-post\LE\UISOP\yearly\'; % 指定批量数据所在的文件夹
filelist6=dir([datadir6,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir6,filelist6(1).name]);
datadir7='E:\CESM-post\LE\USUBM\yearly\'; % 指定批量数据所在的文件夹
filelist7=dir([datadir7,'*192001-200512.nc']); % 指定批量数据的类型
ncdisp([datadir7,filelist7(1).name]);
k=length(filelist1)
clear OHCz* MOHT*1
%
for s=1:40
    s 
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    VVEL = ncread(filename2,'VVEL')/100; % cm/s -> m/s
    filename3=[datadir3,filelist3(s).name];
    ncid=netcdf.open(filename3,'NC_NOWRITE');
    VISOP = ncread(filename3,'VISOP')/100; % cm/s -> m/s
    filename4=[datadir4,filelist4(s).name];
    ncid=netcdf.open(filename4,'NC_NOWRITE');
    VSUBM = ncread(filename4,'VSUBM')/100; % cm/s -> m/s
    Vvel = VVEL(:,1:90,1:inlev,1:86)+VISOP(:,1:90,1:inlev,1:86)+VSUBM(:,1:90,1:inlev,1:86); % v velocity
    clear VVEL VISOP VSUBM
    filename5=[datadir5,filelist5(s).name];
    ncid=netcdf.open(filename5,'NC_NOWRITE');
    UVEL = ncread(filename5,'UVEL')/100; % cm/s -> m/s
    filename6=[datadir6,filelist6(s).name];
    ncid=netcdf.open(filename6,'NC_NOWRITE');
    UISOP = ncread(filename6,'UISOP')/100; % cm/s -> m/s
    filename7=[datadir7,filelist7(s).name];
    ncid=netcdf.open(filename7,'NC_NOWRITE');
    USUBM = ncread(filename7,'USUBM')/100; % cm/s -> m/s
    Uvel = UVEL(:,1:90,1:inlev,1:86)+UISOP(:,1:90,1:inlev,1:86)+USUBM(:,1:90,1:inlev,1:86); % u velocity
    clear UVEL UISOP USUBM
    % temp
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    TemVel = Temp(:,1:90,1:inlev,1:86).*Vvel(:,1:90,1:inlev,1:86); % 1940-2005
    clear Vvel 
    TemUel = Temp(:,1:90,1:inlev,1:86).*Uvel(:,1:90,1:inlev,1:86); % 1940-2005
    clear Uvel Temp 
    dUT(1,:,:,:) = TemUel(1,:,:,:)-TemUel(end,:,:,:);
    dUT(2:360,:,:,:) = TemUel(2:end,:,:,:)-TemUel(1:end-1,:,:,:);
    dUTdx = dUT./londist';
    clear dUT TemUel
    dVT(:,1,:,:) = zeros(360,1,inlev,86);
    dVT(:,2:90,:,:) = TemVel(:,2:90,:,:)-TemVel(:,1:89,:,:);
    dVTdy = dVT./dy;
    clear dVT TemVel
    UVTall = dUTdx+dVTdy;
    UVTalld = (UVTall - UVTMall).*permute(levdist,[3 1 2]);
    clear UVTall 
    UVTallda = UVTalld - nanmean(UVTalld,4);
    UVT_r = cat(1,UVTallda(splon:360,:,:,:),UVTallda(1:splon-1,:,:,:));
    UVTp(:,:,s) = permute(nansum(nansum(UVT_r(1:300-splon,:,:,:),1),3),[2 4 3 1]);
    UVTia(:,:,s) = permute(nansum(nansum(UVT_r(300-splon+1:end,:,:,:),1),3),[2 4 3 1]);        
end




%%
nmem=40; nyr = 86;
UVTrp = reshape(permute(UVTp(1:55,:,:),[2 1 3]),[nyr,55*nmem]);
UVTria = reshape(permute(UVTia(1:55,:,:),[2 1 3]),[nyr,55*nmem]);
clear parp h01
dT = 1;  cf = 1/8;
parfor i = 1:size(UVTrp,2);
    ya = lanczosfilter(UVTrp(:,i),dT,cf,[],'low'); % 8 year filtered IPO index
    ya(1:4) = []; ya(end-3:end) = [];
    [parp(i),hp(i),t] = reg1_ttest(TPIfz(1:78,s),ya,0.1,1);

    yb = lanczosfilter(UVTria(:,i),dT,cf,[],'low'); % 8 year filtered IPO index
    yb(1:4) = []; yb(end-3:end) = [];
    [paria(i),hia(i),t] = reg1_ttest(TPIfz(1:78,s),yb,0.1,1);
end
parpr = reshape(parp,[55,nmem]);
hpr = reshape(hp,[55,nmem]);
pariar = reshape(paria,[55,nmem]);
hiar = reshape(hia,[55,nmem]);
%%
ftsz = 12;
close all;
map1 = mean(parpr,2);
map2 = mean(pariar,2);
ymin1 = -2.5*10^-3; ymax1 = 1.5*10^-3;
Fig = figure('position',[100 100 800 400])
h3 = plot(latData(1:55),smooth(map1,7),'Color',[.64 .08 .18],'linewidth',2)
hold on
h4 = plot(latData(1:55),smooth(map2,7),'Color',[.04 .35 .56],'linewidth',2)
plot(latData(1:55),zeros(55),'k')
% ylabel('Trend (10^4 W m^-^1 Longitude^-^1 decade^-^1)')
plot(-46*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-61*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
set(gca,'XLim',[-75,-35]);set(gca,'XTick',[-70:10:-35]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:0.5*10^-3:ymax1]); set(gca,'YTick',[ymin1:0.5*10^-3:ymax1],'FontSize',14);
ylabel('dUT/dx+dVT/dy')
legend(gca,[h3,h4],'Pacific','Atlantic-Indian Ocean','fontsize',14,'Location','north','Orientation','vertical','position',[0.2,0.82,0.18,0.008])
legend('boxoff')
print(Fig,['D:\figures\CESM\Yearly\LE\20240130_IPO_SO\dUTdx+dVTdy.png'],'-dpng','-r300')



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
function contourfSPolar(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end