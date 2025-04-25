clc,clear;
addpath(genpath('/Volumes/Togo4T/1_matlab/help/'));
load("/Volumes/Togo4T/1_matlab/MatData/CESM/lonData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/latData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
%% MOHT calculate 
inlev = 11; depthstr = '0-100m';
lats = [1:90];
f = sw_f(latData(lats));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(lats)*pi/180)/360; % distance between 2 longitudes 
levdist = gradient(depthData(1:inlev)); % z distance
%%%%%%%%%%%%%%%%%%%%% ensmean MOHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEMP1 = ncread("/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_1920-2005.nc",'TEMP');
TEMP2 = ncread("/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_2006-2022.nc",'TEMP');
TEMP = cat(4,TEMP1(:,:,1:inlev,:),TEMP2(:,:,1:inlev,1:16));
clear TEMP1 TEMP2
VVEL1 = ncread("/Volumes/CESM-post2/CESM-post/LE/VVEL/VVELensmean_1920-2005.nc",'VVEL');
VVEL2 = ncread("/Volumes/CESM-post2/CESM-post/LE/VVEL/VVELensmean_2006-2080.nc",'VVEL');
VVEL = cat(4,VVEL1,VVEL2(:,:,:,1:16));
clear VVEL1 VVEL2
VISOP1 = ncread("/Volumes/CESM-post2/CESM-post/LE/VISOP/VISOPensmean_1920-2005.nc",'VISOP');
VISOP2 = ncread("/Volumes/CESM-post2/CESM-post/LE/VISOP/VISOPensmean_2006-2080.nc",'VISOP');
VISOP = cat(4,VISOP1,VISOP2(:,:,:,1:16));
clear VISOP1 VISOP2
VSUBM1 = ncread("/Volumes/CESM-post2/CESM-post/LE/VSUBM/VSUBMensmean_1920-2005.nc",'VSUBM');
VSUBM2 = ncread("/Volumes/CESM-post2/CESM-post/LE/VSUBM/VSUBMensmean_2006-2080.nc",'VSUBM');
VSUBM = cat(4,VSUBM1,VSUBM2(:,:,:,1:16));
clear VSUBM1 VSUBM2
VRES = (VVEL+VISOP+VSUBM)/100; % cm->m
clear VVEL VISOP VSUBM
%%
TemVel = TEMP(:,lats,:,:).*VRES(:,lats,1:inlev,:);
OHCz = permute(sum((TemVel(:,:,1:inlev,:)).*londist'.*permute(levdist,[3 2 1]),3),[1 2 4 3]); % integrated along z
OHCzy = ro*cp*OHCz;
OHCzy_r = cat(1,OHCzy(151:360,:,:),OHCzy(1:150,:,:));
lonData_r = cat(1,lonData(151:360)-360,lonData(1:150));
MOHTp = permute(nansum(OHCzy_r(1:140,:,:),1),[2 3 4 1]); % Pac 150E-70W unit:W
MOHTia = permute(nansum(OHCzy_r(141:360,:,:),1),[2 3 4 1]); % Atl-IO 70W-150E unit:W
%%
startyr = 1960;
endyr = 2020;
lats = [20:60];
var1 = MOHTp(lats,startyr-1919:endyr-1919)'; % W m^-^1
var2 = MOHTia(lats,startyr-1919:endyr-1919)';
cli1 = mean(var1,1);
cli2 = mean(var2,1);
x = [1:size(var1,1)]';
clear trd1
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
%% MOHT/dlon -> W/m-2

elon = 1000*sw_dist([45,45],[1,2],'km');
map1 = 10*trd1/(140*elon); 
map2 = 10*trd2/(210*elon); 
map1c = cli1/(140*elon);
map2c = cli2/(210*elon);
% map = trd3/(90*elon); regstr = 'IO_30E-120E'
close all;
ftsz = 20;
Fig = figure('position',[700 100 800 400]); 
ymax1 = 10; ymin1 = -5; yint1 = 2.5; 
ymax2 = 10; ymin2 = -5; yint2 = 2.5; 
p1 = plot(latData(20:60),map1/10^5,'-','color',[.04 .35 .56],'LineWidth',3);
hold on
p2 = plot(latData(20:60),map2/10^5,'-','color',[.64 .08 .18],'LineWidth',3);
plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint1:ymax1]);
ylabel('Trend  (10^5 W m^-^2 decade^-^1)','FontSize',ftsz);

yyaxis right
p3 = plot(latData(20:60),map1c/10^7,'--','color',[.04 .35 .56],'LineWidth',1.5);
hold on
p4 = plot(latData(20:60),map2c/10^7,'--','color',[.64 .08 .18],'LineWidth',1.5);
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:yint2:ymax2],'ycolor','k');
ylabel('Climatology  (10^7 W m^-^2)','FontSize',ftsz,'Color','k');
set(gca,'XLim',[-70,-30],'xtick',[-70:10:-30],'fontsize',ftsz,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
set(gca,'XGrid','on','YGrid','on');
legend([p1, p3, p2, p4],['Pacific trend'],['Pacific climatology'],['Atlantic-Indian trend'],...
    ['Atlantic-Indian climatology'],'Location','northwest','numcolumns',2,'edgecolor','w')

% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240503_SO_Warming/MOHT_VRES_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
% 手动哦保存 MOHTdx_0-100m_1960-2020.png
%% dMOHT/dy m-2

dMdyp = dMOHTdy(MOHTp(lats,startyr-1919:endyr-1919));
dMdyia = dMOHTdy(MOHTia(lats,startyr-1919:endyr-1919));
cli1 = nanmean(dMdyp,2);
cli2 = nanmean(dMdyia,2);
%
var1 = dMdyp';
var2 = dMdyia';
x = [1:size(var1,1)]';
clear trd1
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
%% MOHT/dlon -> W/m-3
elon = 1000*sw_dist([45,45],[1,2],'km');

% Pacific
% map1 = trd1*10/(140*elon); 
% map2 = cli1/(140*elon);

% Atlantic-Indian Ocean
map1 = trd2*10/(210*elon); 
map2 = cli2/(210*elon);

ymin1 = -2.4; ymax1 = 2.4; yint = .8
ymin2 = -90; ymax2 = 90; yint2 = 30
latd = latData(lats);
close all;
ftsz = 20;
Fig = figure('position',[100 100 800 400])
nump = find(map1 > 0); numn = find(map1 < 0);
h1 = bar(latd(nump),map1(nump),'FaceColor',[.04 .35 .56],'EdgeColor',[.04 .35 .56],'barwidth',1)
hold on
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint:ymax1]);
h2 = bar(latd(numn),map1(numn),'FaceColor',[.64 .08 .18],'EdgeColor',[.64 .08 .18],'barwidth',1)
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
plot(-35.5*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55.5*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
set(gca,'XLim',[-70.5,-30.5],'XTick',[-70.5:10:-30.5],'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',ftsz,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',ftsz),ylabel('Trend (W m^-^3 decade^-^1)','FontSize',ftsz);

yyaxis right
numpc = find(map2 > 0); numnc = find(map2 < 0);
h3 = bar(latd(numpc),map2(numpc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
h4 = bar(latd(numnc),map2(numnc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:yint2:ymax2],'ycolor','k');
ylabel('Climatology (W m^-^3)','FontSize',ftsz);
legend([h1, h2],'Divergence trend','Convergence trend','Location','northwest','edgecolor','w')

% 手动保存 MOHTdydx_0-100m_Pac_1960-2020.png
%% VISOP
% historical 
datadir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='/Volumes/CESM-post2/CESM-post/LE/VISOP/yearly/'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1)
clear OHCz* MOHT*1
lonData_r = cat(1,lonData(180:360),lonData(1:179));
%%
splon1=151; % 150E 
splon2=291;  % 70W
for s=1:40
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VISOP')/100; % cm/s -> m/s
    % 0-700m MOHT
    TemVel = Temp(:,1:90,1:inlev,41:86).*Vvel(:,1:90,1:inlev,41:86); % 1960-2005
    clear Temp Vvel
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(splon1:360,:,:,:),OHCzy(1:splon1-1,:,:,:));
    MOHTp1(:,:,s) = permute(nansum(OHCzy_r(1:splon2-splon1,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTa1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:splon2-splon1+101,:,:,:),1),[2 3 4 1]); % Atl 70W-30E unit:W
    MOHTi1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+101:end,:,:,:),1),[2 3 4 1]); % io 30E-150E unit:W
    MOHTall1(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
MOHTp1LE = nanmean(MOHTp1,3);
MOHTia1LE = nanmean(MOHTia1,3);
MOHTi1LE = nanmean(MOHTi1,3);
MOHTa1LE = nanmean(MOHTa1,3);
MOHTall1LE = nanmean(MOHTall1,3);
clear Temp Vvel TemVel
% rcp85 2006-2022
datadir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*200601-*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='/Volumes/CESM-post2/CESM-post/LE/VISOP/yearly/'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*200601-*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1)
clear OHCz* 
%
for s=1:k
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VISOP')/100; % cm/s -> m/s
    % 0-2000m MOHT
    TemVel = Temp(:,1:90,1:inlev,1:15).*Vvel(:,1:90,1:inlev,1:15); % 2006-2020
    clear Temp Vvel
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(splon1:360,:,:,:),OHCzy(1:splon1-1,:,:,:));
    MOHTp2(:,:,s) = permute(nansum(OHCzy_r(1:splon2-splon1,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTa2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:splon2-splon1+101,:,:,:),1),[2 3 4 1]); % Atl 70W-30E unit:W
    MOHTi2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+101:end,:,:,:),1),[2 3 4 1]); % io 30E-150E unit:W
    MOHTall2(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
MOHTp2LE = nanmean(MOHTp2,3);
MOHTia2LE = nanmean(MOHTia2,3);
MOHTi2LE = nanmean(MOHTi2,3);
MOHTa2LE = nanmean(MOHTa2,3);
MOHTall2LE = nanmean(MOHTall2,3);
clear Temp Vvel TemVel
% trend and climo
MOHTp = double(cat(2,MOHTp1LE,MOHTp2LE)); % MOHT ensemble mean
MOHTia = double(cat(2,MOHTia1LE,MOHTia2LE));
MOHTi = double(cat(2,MOHTi1LE,MOHTi2LE));
MOHTa = double(cat(2,MOHTa1LE,MOHTa2LE));
MOHTall = double(cat(2,MOHTall1LE,MOHTall2LE));
MOHTcp = nanmean(MOHTp,2); % climo
MOHTcia = nanmean(MOHTia,2);
MOHTci = nanmean(MOHTi,2);
MOHTca = nanmean(MOHTa,2);
startyr = 1960;
endyr = 2020;
var1 = permute(MOHTp(:,startyr-1959:endyr-1959),[2 1]);
var2 = permute(MOHTia(:,startyr-1959:endyr-1959),[2 1]);
var3 = permute(MOHTi(:,startyr-1959:endyr-1959),[2 1]);
var4 = permute(MOHTa(:,startyr-1959:endyr-1959),[2 1]);
x = [1:size(var1,1)]';
clear trd*
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd1(i) = par1(1); 
        par2=polyfit(x,var2(:,i),1); % regression parameters
        trd2(i) = par2(1); 
        par3=polyfit(x,var3(:,i),1); % regression parameters
        trd3(i) = par3(1); 
        par4=polyfit(x,var4(:,i),1); % regression parameters
        trd4(i) = par4(1); 
end
VISOPcp=MOHTcp; VISOPtd1 = trd1; VISOPp = MOHTp;
VISOPcia=MOHTcia; VISOPtd2 = trd2; VISOPia = MOHTia;
VISOPci=MOHTci; VISOPtd3 = trd3; VISOPi = MOHTi;
VISOPca=MOHTca; VISOPtd4 = trd4; VISOPa = MOHTa;

%% VSUBM
% historical 
datadir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='/Volumes/CESM-post2/CESM-post/LE/VSUBM/yearly/'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1)
clear OHCz* MOHT*1
for s=1:40
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VSUBM')/100; % cm/s -> m/s
    % 0-700m MOHT
    TemVel = Temp(:,1:90,1:inlev,41:86).*Vvel(:,1:90,1:inlev,41:86); % 1960-2005
    clear Temp Vvel
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(splon1:360,:,:,:),OHCzy(1:splon1-1,:,:,:));
    MOHTp1(:,:,s) = permute(nansum(OHCzy_r(1:splon2-splon1,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTa1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:splon2-splon1+101,:,:,:),1),[2 3 4 1]); % Atl 70W-30E unit:W
    MOHTi1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+101:end,:,:,:),1),[2 3 4 1]); % io 30E-150E unit:W
    MOHTall1(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
MOHTp1LE = nanmean(MOHTp1,3);
MOHTia1LE = nanmean(MOHTia1,3);
MOHTi1LE = nanmean(MOHTi1,3);
MOHTa1LE = nanmean(MOHTa1,3);
MOHTall1LE = nanmean(MOHTall1,3);
clear Temp Vvel TemVel
%% rcp85 2006-2022
datadir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*200601-*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='/Volumes/CESM-post2/CESM-post/LE/VSUBM/yearly/'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*200601-*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1)
clear OHCz* 
%
for s=1:k
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VSUBM')/100; % cm/s -> m/s
    % 0-2000m MOHT
    TemVel = Temp(:,1:90,1:inlev,1:15).*Vvel(:,1:90,1:inlev,1:15); % 2006-2020
    clear Temp Vvel
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(splon1:360,:,:,:),OHCzy(1:splon1-1,:,:,:));
    MOHTp2(:,:,s) = permute(nansum(OHCzy_r(1:splon2-splon1,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTa2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:splon2-splon1+101,:,:,:),1),[2 3 4 1]); % Atl 70W-30E unit:W
    MOHTi2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+101:end,:,:,:),1),[2 3 4 1]); % io 30E-150E unit:W
    MOHTall2(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
MOHTp2LE = nanmean(MOHTp2,3);
MOHTia2LE = nanmean(MOHTia2,3);
MOHTi2LE = nanmean(MOHTi2,3);
MOHTa2LE = nanmean(MOHTa2,3);
MOHTall2LE = nanmean(MOHTall2,3);
clear Temp Vvel TemVel
% trend and climo
MOHTp = double(cat(2,MOHTp1LE,MOHTp2LE)); % MOHT ensemble mean
MOHTia = double(cat(2,MOHTia1LE,MOHTia2LE));
MOHTi = double(cat(2,MOHTi1LE,MOHTi2LE));
MOHTa = double(cat(2,MOHTa1LE,MOHTa2LE));
MOHTall = double(cat(2,MOHTall1LE,MOHTall2LE));
MOHTcp = nanmean(MOHTp,2); % climo
MOHTcia = nanmean(MOHTia,2);
MOHTci = nanmean(MOHTi,2);
MOHTca = nanmean(MOHTa,2);
startyr = 1960;
endyr = 2020;
var1 = permute(MOHTp(:,startyr-1959:endyr-1959),[2 1]);
var2 = permute(MOHTia(:,startyr-1959:endyr-1959),[2 1]);
var3 = permute(MOHTi(:,startyr-1959:endyr-1959),[2 1]);
var4 = permute(MOHTa(:,startyr-1959:endyr-1959),[2 1]);
x = [1:size(var1,1)]';
clear trd1
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd1(i) = par1(1); 
        par2=polyfit(x,var2(:,i),1); % regression parameters
        trd2(i) = par2(1); 
        par3=polyfit(x,var3(:,i),1); % regression parameters
        trd3(i) = par3(1); 
        par4=polyfit(x,var4(:,i),1); % regression parameters
        trd4(i) = par4(1); 
end
VSUBMcp=MOHTcp; VSUBMtd1 = trd1; VSUBMp = MOHTp;
VSUBMcia=MOHTcia; VSUBMtd2 = trd2; VSUBMia = MOHTia;
VSUBMci=MOHTci; VSUBMtd3 = trd3; VSUBMi = MOHTi;
VSUBMca=MOHTca; VSUBMtd4 = trd4; VSUBMa = MOHTa;

%% VVEL
% historical 
datadir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='/Volumes/CESM-post2/CESM-post/LE/VVEL/yearly/'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1);
clear OHCz* MOHT*1
for s=1:40
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VVEL')/100; % cm/s -> m/s
    % 0-700m MOHT
    TemVel = Temp(:,1:90,1:inlev,41:86).*Vvel(:,1:90,1:inlev,41:86); % 1960-2005
    clear Temp Vvel
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(splon1:360,:,:,:),OHCzy(1:splon1-1,:,:,:));
    MOHTp1(:,:,s) = permute(nansum(OHCzy_r(1:splon2-splon1,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTa1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:splon2-splon1+101,:,:,:),1),[2 3 4 1]); % Atl 70W-30E unit:W
    MOHTi1(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+101:end,:,:,:),1),[2 3 4 1]); % io 30E-150E unit:W
    MOHTall1(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
MOHTp1LE = nanmean(MOHTp1,3);
MOHTia1LE = nanmean(MOHTia1,3);
MOHTi1LE = nanmean(MOHTi1,3);
MOHTa1LE = nanmean(MOHTa1,3);
MOHTall1LE = nanmean(MOHTall1,3);
clear Temp Vvel TemVel
%% rcp85 2006-2022
datadir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*200601-*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='/Volumes/CESM-post2/CESM-post/LE/VVEL/yearly/'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*200601-*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1);
clear OHCz* 
%
for s=1:k
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VVEL')/100; % cm/s -> m/s
    % 0-2000m MOHT
    TemVel = Temp(:,1:90,1:inlev,1:15).*Vvel(:,1:90,1:inlev,1:15); % 2006-2020
    clear Temp Vvel
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(splon1:360,:,:,:),OHCzy(1:splon1-1,:,:,:));
    MOHTp2(:,:,s) = permute(nansum(OHCzy_r(1:splon2-splon1,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTa2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+1:splon2-splon1+101,:,:,:),1),[2 3 4 1]); % Atl 70W-30E unit:W
    MOHTi2(:,:,s) = permute(nansum(OHCzy_r(splon2-splon1+101:end,:,:,:),1),[2 3 4 1]); % io 30E-150E unit:W
    MOHTall2(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
MOHTp2LE = nanmean(MOHTp2,3);
MOHTia2LE = nanmean(MOHTia2,3);
MOHTi2LE = nanmean(MOHTi2,3);
MOHTa2LE = nanmean(MOHTa2,3);
MOHTall2LE = nanmean(MOHTall2,3);
clear Temp Vvel TemVel
% trend and climo
MOHTp = double(cat(2,MOHTp1LE,MOHTp2LE)); % MOHT ensemble mean
MOHTia = double(cat(2,MOHTia1LE,MOHTia2LE));
MOHTi = double(cat(2,MOHTi1LE,MOHTi2LE));
MOHTa = double(cat(2,MOHTa1LE,MOHTa2LE));
MOHTall = double(cat(2,MOHTall1LE,MOHTall2LE));
MOHTcp = nanmean(MOHTp,2); % climo
MOHTcia = nanmean(MOHTia,2);
MOHTci = nanmean(MOHTi,2);
MOHTca = nanmean(MOHTa,2);
startyr = 1960;
endyr = 2020;
var1 = permute(MOHTp(:,startyr-1959:endyr-1959),[2 1]);
var2 = permute(MOHTia(:,startyr-1959:endyr-1959),[2 1]);
var3 = permute(MOHTi(:,startyr-1959:endyr-1959),[2 1]);
var4 = permute(MOHTa(:,startyr-1959:endyr-1959),[2 1]);
x = [1:size(var1,1)]';
clear trd1
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd1(i) = par1(1); 
        par2=polyfit(x,var2(:,i),1); % regression parameters
        trd2(i) = par2(1); 
        par3=polyfit(x,var3(:,i),1); % regression parameters
        trd3(i) = par3(1); 
        par4=polyfit(x,var4(:,i),1); % regression parameters
        trd4(i) = par4(1); 
end
%
close all;
plot(trd1,'r')
hold on
plot(trd2,'b')
VVELcp=MOHTcp; VVELtd1 = trd1; VVELp = MOHTp;
VVELcia=MOHTcia; VVELtd2 = trd2; VVELia = MOHTia;
VVELci=MOHTci; VVELtd3 = trd3; VVELi = MOHTi;
VVELca=MOHTca; VVELtd4 = trd4; VVELa = MOHTa;
%% save
VISOPc = cat(2,VISOPcp,VISOPcia,VISOPci,VISOPca);
VISOPtrd = cat(2,VISOPtd1',VISOPtd2',VISOPtd3',VISOPtd4');
VISOPt = cat(3,VISOPp,VISOPia,VISOPi,VISOPa);
VSUBMc = cat(2,VSUBMcp,VSUBMcia,VSUBMci,VSUBMca);
VSUBMtrd = cat(2,VSUBMtd1',VSUBMtd2',VSUBMtd3',VSUBMtd4');
VSUBMt = cat(3,VSUBMp,VSUBMia,VSUBMi,VSUBMa);
VVELc = cat(2,VVELcp,VVELcia,VVELci,VVELca);
VVELtrd = cat(2,VVELtd1',VVELtd2',VVELtd3',VVELtd4');
VVELt = cat(3,VVELp,VVELia,VVELi,VVELa);
%
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VVELc.mat'],'VVELc');
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VVELtrd.mat'],'VVELtrd');
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VVELt.mat'],'VVELt');
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VISOPc.mat'],'VISOPc');
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VISOPtrd.mat'],'VISOPtrd');
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VISOPt.mat'],'VISOPt');
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VSUBMc.mat'],'VSUBMc');
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VSUBMtrd.mat'],'VSUBMtrd');
save(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VSUBMt.mat'],'VSUBMt');


%% load
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VVELc.mat'],'VVELc');
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VVELtrd.mat'],'VVELtrd');
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VVELt.mat'],'VVELt');
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VISOPc.mat'],'VISOPc');
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VISOPtrd.mat'],'VISOPtrd');
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VISOPt.mat'],'VISOPt');
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VSUBMc.mat'],'VSUBMc');
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VSUBMtrd.mat'],'VSUBMtrd');
load(['/Volumes/Togo4T/1_matlab/MatData/CESM/MOHT_150E_70W_700m/VSUBMt.mat'],'VSUBMt');
%%
VRESc = VVELc + VISOPc + VSUBMc;
VREStrd = VVELtrd + VISOPtrd + VSUBMtrd;
VRESt = VVELt + VISOPt + VSUBMt;
%% MOHT 

% map1 = VRESc(:,1); map2 = VREStrd(:,1); clor = [.04 .35 .56];
% regstr1 = 'Pacific'; regstr2 = 'Pac'; 
% 
map1 = VRESc(:,2); map2 = VREStrd(:,2); clor = [.64 .08 .18];
regstr1 = 'Atlantic-Indian Ocean'; regstr2 = 'AtlIO';
% % % % % % % 
% map1 = VRESc(:,3); map2 = VREStrd(:,3); clor = [.64 .08 .18];
% regstr1 = 'Indian Ocean'; regstr2 = 'IO';
% % % % % % % 
% map1 = VRESc(:,4); map2 = VREStrd(:,4); clor = [.64 .08 .18];
% regstr1 = 'Atlantic'; regstr2 = 'Atl';

close all;
Fig = figure('position',[700 100 800 400]);
ymin = -2; ymax = 2.5; yint = 0.5;
p2 = plot(latData(10:60),map2(10:60)*10/10^13,'-','color',clor,'LineWidth',3);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);
hold on
set(gca,'XLim',[-80,-30],'xtick',[-80:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('Trend  (10^-^2 PW decade^-^1)','FontSize',14);
plot(-35*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
p1 = plot(latData(10:60),map1(10:60)/10^15,'-','color','k','LineWidth',3);
hold on
% p4 = plot(climap1,'-','color','k','LineWidth',3);
plot(latData(10:60),zeros(1,51),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12),ylabel('Climo (PW)','FontSize',14);
legend([p1, p2],['Climo (',regstr1,')'],['Trend (',regstr1,')'],'Location','northwest')
legend('boxoff')

yyaxis left
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240503_SO_Warming/MOHT_VRES_',regstr2,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% MOHT/dlon -> W/m-1
% lonstr = '150E_70W'
% dlon = splon2-splon1;
% elon = 1000*sw_dist([45,45],[1,2],'km');
% dx = dlon*elon;
% map1 = VRESc(:,1)/dx; map2 = VREStrd(:,1)/dx; clor = [.04 .35 .56];
% regstr1 = 'Pacific'; regstr2 = 'Pac'; 
% % 
% dlon = 360-(300-splon);
% elon = 1000*sw_dist([45,45],[1,2],'km');
% dx = dlon*elon;
% map1 = VRESc(:,2)/dx; map2 = VREStrd(:,2)/dx; clor = [.64 .08 .18];
% regstr1 = 'Atlantic-Indian Ocean'; regstr2 = 'AtlIO';
% % % % % % % % % % % 
% dlon = 120;
% elon = 1000*sw_dist([45,45],[1,2],'km');
% dx = dlon*elon;
% map1 = VRESc(:,3)/dx; map2 = VREStrd(:,3)/dx; clor = [.64 .08 .18];
% regstr1 = 'Indian Ocean'; regstr2 = 'IO';
% % % % % % % % % %  
dlon = 100;
elon = 1000*sw_dist([45,45],[1,2],'km');
dx = dlon*elon
map1 = VRESc(:,4)/dx; map2 = VREStrd(:,4)/dx; clor = [.64 .08 .18];
regstr1 = 'Atlantic'; regstr2 = 'Atl';

close all;
Fig = figure('position',[700 100 800 400]);
% ymax = 1; ymin = -1; yint = 0.5; % 100m
% ymax = 2.5; ymin = -2.5; yint = 0.5; % 300m
ymax = 320; ymin = -320; yint = 80; % 700m
p2 = plot(latData(20:60),map2(20:60)/10^3,'-','color',clor,'LineWidth',3);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);
hold on
set(gca,'XLim',[-70,-30],'xtick',[-70:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('Trend  (kW m^-^1 yr^-^1)','FontSize',14);
plot(-35*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
p1 = plot(latData(20:60),map1(20:60)/10^6,'-','color','k','LineWidth',3);
hold on
% p4 = plot(climap1,'-','color','k','LineWidth',3);
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12),ylabel('Climo (10^3 kW m^-^1)','FontSize',14);
legend([p1, p2],['Climo (',regstr1,')'],['Trend (',regstr1,')'],'Location','northwest')
legend('boxoff')

yyaxis left
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);

print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240503_SO_Warming/MOHTdlon_VRES_',regstr2,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% dMOHT/dy m-1
% dlon = splon2-splon1; regstr = 'Pac';
% elon = 1000*sw_dist([45,45],[1,2],'km');
% dx = dlon*elon;
% dMdy = dMOHTdy(VRESt(:,:,1));
% dMdycli = dMOHTdy(VRESc(:,1));
% ymin1 = -.9; ymax1 = .9;  ymin2 = -900; ymax2 = 900;
% ymin3 = -3; ymax3 = 3;  ymin4 = -.9; ymax4 = .9;
% 
% dlon = 220; regstr = 'AtlIO';
% dx = dlon*elon;
% dMdy = dMOHTdy(VRESt(:,:,2));
% dMdycli = dMOHTdy(VRESc(:,2));
% ymin1 = -.9; ymax1 = .9;  ymin2 = -900; ymax2 = 900;
% ymin3 = -3; ymax3 = 3;  ymin4 = -.9; ymax4 = .9;
% % % % % % % % 
% dlon = 120; regstr = 'IO';
% dx = dlon*elon;
% dMdy = dMOHTdy(VRESt(:,:,3));
% dMdycli = dMOHTdy(VRESc(:,3));
% ymin1 = -.9; ymax1 = .9;  ymin2 = -900; ymax2 = 900;
% ymin3 = -3; ymax3 = 3;  ymin4 = -3; ymax4 = 3;
% % % % % % %
dlon = 100; regstr = 'Atl';
dMdy = dMOHTdy(VRESt(:,:,4));
dMdycli = dMOHTdy(VRESc(:,4));
ymin1 = -.9; ymax1 = .9;  ymin2 = -900; ymax2 = 900;
ymin3 = -3; ymax3 = 3;  ymin4 = -3; ymax4 = 3;

startyr = 1960; endyr = 2020;
var1 = permute(dMdy(:,startyr-1959:endyr-1959),[2 1]);
x = [1:size(var1,1)]';
clear trd
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd(i) = par1(1); 
end
% dMOHT/dy/dx
close all;
Fig = figure('position',[100 100 800 400])
map1 = trd(20:60)/dx;
latd = latData(20:60);
nump = find(map1 > 0); numn = find(map1 < 0);
h3 = bar(latd(nump),map1(nump),'FaceColor',[.64 .08 .18],'EdgeColor',[.64 .08 .18],'barwidth',1)
hold on
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:ymax1/3:ymax1]);
bar(latd(numn),map1(numn),'FaceColor',[.04 .35 .56],'EdgeColor',[.04 .35 .56],'barwidth',1)
ylabel('Trend (W m^-^2 yr^-^1)')
plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
map2 = dMdycli(20:60)/dx;
numpc = find(map2 > 0); numnc = find(map2 < 0);
h1 = bar(latd(numpc),map2(numpc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
h2 = bar(latd(numnc),map2(numnc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
plot(latData(20:60),zeros(1,41),'k','LineWidth',1)
set(gca,'XLim',[-70,-30],'XTick',[-70:10:-30],'XTickLabel',[-70:10:30],'FontSize',14,'XGrid','on','YGrid','on');
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:ymax2/3:ymax2],'YColor',[.5 .5 .5]);
xlabel('Latitude','FontSize',14),ylabel('Climo (W m^-^2)','FontSize',14);
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240503_SO_Warming/MOHTdydx_VRES_',regstr,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% dMOHT/dy
% close all;
Fig = figure('position',[100 100 800 400])
map1 = trd(1:60)*10/10^7;
nump = find(map1 > 0); numn = find(map1 < 0);
h3 = bar(latData(nump),map1(nump),'FaceColor',[.64 .08 .18],'EdgeColor',[.64 .08 .18],'barwidth',1)
hold on
set(gca,'YLim',[ymin3,ymax3],'YTick',[ymin3:ymax3/3:ymax3]);
bar(latData(numn),map1(numn),'FaceColor',[.04 .35 .56],'EdgeColor',[.04 .35 .56],'barwidth',1)
ylabel('Trend (10^-^2 GW m^-^1 decade^-^1)')
plot(-35*ones(11),[ymin3:(ymax3-ymin3)/10:ymax3],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-50*ones(11),[ymin3:(ymax3-ymin3)/10:ymax3],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
map2 = dMdycli(1:60)/10^9;
numpc = find(map2 > 0); numnc = find(map2 < 0);
h1 = bar(latData(numpc),map2(numpc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
h2 = bar(latData(numnc),map2(numnc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
plot(latData(1:60),zeros(1,60),'k','LineWidth',1)
set(gca,'XLim',[-80,-30],'XTick',[-80:10:-30],'XTickLabel',[-80:10:30],'FontSize',14,'XGrid','on','YGrid','on');
set(gca,'YLim',[ymin4,ymax4],'YTick',[ymin4:ymax4/3:ymax4],'YColor',[.5 .5 .5]);
xlabel('Latitude','FontSize',14),ylabel('Climo (GW m^-^1)','FontSize',14);
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240120_SO_Warming/MOHTdy_VRES_',regstr,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')



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
    load('/Volumes/Togo/1_matlab/help/colorbar_mat/bl_re4.mat');
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
