
clc,clear;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
load("/Volumes/Togo4T/1_matlab/MatData/LICOM/lonData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/LICOM/latData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/LICOM/depthData.mat");

for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
end

[ctrl_p ctrl_ia ctrl_zm] = LICOM_OHCzm('/Volumes/Togo4T/data/LICOMpost/Tem_2000m.era5.grid.nc',dx);
[heat_p heat_ia heat_zm] = LICOM_OHCzm('/Volumes/Togo4T/data/LICOMpost/Tem_2000m.era5.HEAT.grid.nc',dx);
[wind_p wind_ia wind_zm] = LICOM_OHCzm('/Volumes/Togo4T/data/LICOMpost/Tem_2000m.era5.WND.grid.nc',dx);
%% average of OHC, is equal to integration of 0-700m at each longitude and latitude
% str = 'Pacific_4'
% trd1 = ctrl_p./(140*dx(1:90)*1000);
% trd2 = heat_p./(140*dx(1:90)*1000);
% trd3 = wind_p./(140*dx(1:90)*1000);
% trd4 = (heat_p+wind_p)./(140*dx(1:90)*1000);

str = 'Atl-IO_3'
trd1 = ctrl_ia./(210*dx(1:90)*1000);
trd2 = heat_ia./(210*dx(1:90)*1000);
trd3 = wind_ia./(210*dx(1:90)*1000);
% trd4 = (heat_ia+wind_ia)./(210*dx(1:90)*1000);

% str = 'ZonalM_3'
% trd1 = ctrl_zm./(350*dx(1:90)*1000);
% trd2 = heat_zm./(350*dx(1:90)*1000);
% trd3 = wind_zm./(350*dx(1:90)*1000);
% trd4 = (heat_zm+wind_zm)./(350*dx(1:90)*1000);

close all;
ftsz = 20;
Fig = figure('position',[700 100 800 400]);
ymin1 = -.3; ymax1 = .9; yint1 = .3;
p1 = plot(latData(20:60),trd1(20:60),'-','color','k','LineWidth',3);
hold on
p2 = plot(latData(20:60),trd2(20:60),'-','color',[.64 .08 .18],'LineWidth',3);
p3 = plot(latData(20:60),trd3(20:60),'-','color',[.04 .35 .56],'LineWidth',3);
% p4 = plot(latData(20:60),trd4(20:60),'--','color','k','LineWidth',3);
plot(latData(20:60),zeros(length(latData(20:60))),'-k','LineWidth',1.5);
% plot(latData(find(trd1 == max(trd1)))*ones(1,21),[ymin1:(ymax1-ymin1)/20:ymax1],'-','color','k','LineWidth',1);
% plot(latData(find(trd2 == max(trd2)))*ones(1,21),[ymin1:(ymax1-ymin1)/20:ymax1],'-','color',[.64 .08 .18],'LineWidth',1);
% plot(latData(find(trd3 == max(trd3(1:50))))*ones(1,21),[ymin1:(ymax1-ymin1)/20:ymax1],'-','color',[.04 .35 .56],'LineWidth',1);
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint1:ymax1],'YColor','k');
set(gca,'XLim',[-70.5,-30.5],'xtick',[-70.5:10:-30.5],'fontsize',ftsz,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('OHC Trend  (W/yr)','FontSize',ftsz);
plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
xlabel('Latitude','FontSize',ftsz),ylabel('Trend (W m^-^2)','FontSize',ftsz);
legend([p1, p2, p3],['CTRL'],['HEAT'],['WND'],'Location','northwest','numcolumns',1,'edgecolor','w')
% legend([p1, p2, p3, p4],['CTRL'],['HEAT'],['WND'],['Sum of HEAT & WND'],'Location','northwest','numcolumns',1,'edgecolor','w')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/LICOM/Yearly/20241021_SO_warming/OHCzm_perDx_',str,'.png'],'-dpng','-r300')
%% OHC integrated 
% str = 'Pacific_4'
% trd1 = ctrl_p/10^7;
% trd2 = heat_p/10^7;
% trd3 = wind_p/10^7;
% trd4 = (heat_p+wind_p)/10^7;

% str = 'Atl-IO_4'
% trd1 = ctrl_ia/10^7;
% trd2 = heat_ia/10^7;
% trd3 = wind_ia/10^7;
% trd4 = (heat_ia+wind_ia)/10^7;

str = 'ZonalM_4'
trd1 = ctrl_zm/10^7;
trd2 = heat_zm/10^7;
trd3 = wind_zm/10^7;
trd4 = (heat_zm+wind_zm)/10^7;

close all;
ftsz = 20;
Fig = figure('position',[700 100 800 400]);
ymin1 = -.5; ymax1 = 2; yint1 = .5;
p1 = plot(latData(20:60),trd1(20:60),'-','color','k','LineWidth',3);
hold on
p2 = plot(latData(20:60),trd2(20:60),'-','color',[.64 .08 .18],'LineWidth',3);
p3 = plot(latData(20:60),trd3(20:60),'-','color',[.04 .35 .56],'LineWidth',3);
p4 = plot(latData(20:60),trd4(20:60),'--','color','k','LineWidth',3);
plot(latData(20:60),zeros(length(latData(20:60))),'-k','LineWidth',1.5);
% plot(latData(find(trd1 == max(trd1)))*ones(1,21),[ymin1:(ymax1-ymin1)/20:ymax1],'-','color','k','LineWidth',1);
% plot(latData(find(trd2 == max(trd2)))*ones(1,21),[ymin1:(ymax1-ymin1)/20:ymax1],'-','color',[.64 .08 .18],'LineWidth',1);
% plot(latData(find(trd3 == max(trd3(1:50))))*ones(1,21),[ymin1:(ymax1-ymin1)/20:ymax1],'-','color',[.04 .35 .56],'LineWidth',1);
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint1:ymax1],'YColor','k');
set(gca,'XLim',[-70.5,-30.5],'xtick',[-70.5:10:-30.5],'fontsize',ftsz,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('OHC Trend  (W/yr)','FontSize',ftsz);
plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-55*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
xlabel('Latitude','FontSize',ftsz),ylabel('Trend (10^7 W m^-^1)','FontSize',ftsz);
% legend([p1, p2, p3],['CTRL'],['HEAT'],['WND'],'Location','northwest','numcolumns',1,'edgecolor','w')
legend([p1, p2, p3, p4],['CTRL'],['HEAT'],['WND'],['Sum of HEAT & WND'],'Location','northwest','numcolumns',1,'edgecolor','w')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/LICOM/Yearly/20241021_SO_warming/OHCzm_int_',str,'.png'],'-dpng','-r300')





function [map_p map_ia map_zm] = LICOM_OHCzm(filename,dx)
%    filename = '/Volumes/Togo4T/data/LICOMpost/Tem_2000m.era5.WND.grid.nc'
    addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
    load("/Volumes/Togo4T/1_matlab/MatData/LICOM/lonData.mat");
    load("/Volumes/Togo4T/1_matlab/MatData/LICOM/latData.mat");
    load("/Volumes/Togo4T/1_matlab/MatData/LICOM/depthData.mat");
    nlon = length(lonData); nlat = length(latData);
    nlev = length(depthData);
    ncid=netcdf.open(filename,'NOWRITE');
    ncdisp(filename);
    lonData = ncread(filename,'lon');
    latData = ncread(filename,'lat');
    depthData = ncread(filename,'lev');
    Tem1 = ncread(filename,'ts');
    % zonal mean OHC
    cp = 4096 % J (kg oC)
    ro = 1025 % kg m-3
    %------------------------- 0-700m -----------------------------------------
    depthstr = '0-700m';
    depthData_r = -depthData;
    dweit = depthData_r(2:21)-depthData_r(1:20); % depth weight
    Tem700 = Tem1(:,:,21,:)+(Tem1(:,:,22,:)-Tem1(:,:,21,:))*(700-depthData_r(21))/(depthData_r(22)-depthData_r(21));
    splon = 150; lonstr = '150E';
%    dy = 111*1000;
    OHC_0 = cp*ro*dx*1000.*cat(3,Tem1(:,:,1,:)*5,Tem1(:,:,2:21,:).*(permute(dweit,[3 2 1])),Tem700*(700-depthData_r(21))); % J/m2
    %
    OHC_0r = reshape(OHC_0(:,1:90,:,:),360*90*22,82)';
    x = [1:size(OHC_0r,1)]';
    clear OHC_trd_0
    parfor i = 1:size(OHC_0r,2)
        par = polyfit(x,OHC_0r(:,i),1); % regression parameters trend /year
        OHC_trd_0(i) = par(1);
    end
    OHC_trd = reshape(OHC_trd_0,360,90,22);
    OHC_trdr = cat(1,OHC_trd(splon:360,:,:),OHC_trd(1:splon-1,:,:));
    OHC_zm_p = squeeze(nansum(nansum(OHC_trdr(1:290-splon,:,:),1),3)); % Pac unit:W 150E-70w
    OHC_zm_ia = squeeze(nansum(nansum(OHC_trdr(291-splon:end,:,:),1),3)); % IO+Atl 
    OHC_zm = squeeze(nansum(nansum(OHC_trdr,1),3)); % zonal mean

    map_zm = OHC_zm/(60*60*24*365);
    map_p = OHC_zm_p/(60*60*24*365);
    map_ia = OHC_zm_ia/(60*60*24*365);
end