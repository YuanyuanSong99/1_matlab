load('MatFile/Tempadlf');
Tempadlf(:,:,:,1:4) = [];
Tempadlf(:,:,:,end-3:end) = []; % 1944-2016
% Phase 1:  1944-1976; Phase 2: 1977-1997; Phase 3: 1998-2014
Temcom1 = mean(Tempadlf(:,:,:,1:33),4); % Phase 1:  1944-1976
Temcom2 = mean(Tempadlf(:,:,:,34:54),4); % Phase 2: 1977-1997
Temcom3 = mean(Tempadlf(:,:,:,55:71),4); % Phase 3: 1998-2014
Temcom = cat(4,Temcom1,Temcom2,Temcom3);
%% longitude mean (zonal mean)
% Indian Ocean
phases = {'1944-1976','1977-1997','1998-2014'};
lonmI1 = permute(nanmean(Temcom(40:100,83:116,:,:),1),[2 3 4 1]); % 40E-100E, 8S-25N
lonmI2 = permute(nanmean(Temcom(20:146,20:82,:,:),1),[2 3 4 1]);
lonmI = cat(1,lonmI2,lonmI1);
lonmI1m = permute(nanmean(Temp(40:100,83:116,:,:),1),[2 3 4 1]); % climotology
lonmI2m = permute(nanmean(Temp(20:146,20:82,:,:),1),[2 3 4 1]);
lonmIm = cat(1,lonmI2m,lonmI1m);
%
close all;
for k = 1:3
    Fig = figure('position',[100 100 800 400]);
    contourf(latData(20:116,1),-depthData,lonmI(:,:,k)',[-2:0.02:2],'linestyle','none')
    hold on
    caxis([-0.1,0.1])
%     [ac ah] = contour(latData(20:116,1),-depthData,(nanmean(lonmIm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    ch = colorbar;
    set(ch,'Ticks',[-1:0.02:0.1]);
    set(gca,'XLim',[-70,20]);
    set(gca,'XTick',[-60:20:20]);
    set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN'},'FontSize',14);
    set(gca,'YLim',[-2000,0]);
    set(gca,'YTick',[-2000:200:0],'FontSize',14);
    set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
    ylabel('Depth (m)')
    title([phases{k},' Com.zonal mean temperature in Indian Ocean'])
    print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.phase',num2str(k),'.Indian.zonalmean.tem.png'],'-dpng','-r600')
end
%% Pacific 
lonmP1 = permute(nanmean(Temcom(146:293,20:82,:,:),1),[2 3 4 1]); % 146E-67W, 70S-8S
lonmP2 = permute(nanmean(Temcom(100:293,83:101,:,:),1),[2 3 4 1]); % 100E-67W, 8S-10N
lonmP3 = permute(nanmean(Temcom(100:270,102:109,:,:),1),[2 3 4 1]); % 100E-90W, 10N-18N
lonmP4 = permute(nanmean(Temcom(100:255,110:156,:,:),1),[2 3 4 1]); % 100E-105W, 18N-65N
lonmP = cat(1,lonmP1,lonmP2,lonmP3,lonmP4);
max(lonmP,[],'all')
min(lonmP,[],'all')
lonmP1m = permute(nanmean(Temp(146:293,20:82,:,:),1),[2 3 4 1]); % climotology
lonmP2m = permute(nanmean(Temp(100:293,83:101,:,:),1),[2 3 4 1]); 
lonmP3m = permute(nanmean(Temp(100:270,102:109,:,:),1),[2 3 4 1]);
lonmP4m = permute(nanmean(Temp(100:255,110:156,:,:),1),[2 3 4 1]); 
lonmPm = cat(1,lonmP1m,lonmP2m,lonmP3m,lonmP4m);
%
close all;
for k = 1:3
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:156,1),-depthData,lonmP(:,:,k)',[-0.5:0.01:0.5],'linestyle','none')
caxis([-0.1,0.1])
hold on
[ac ah] = contour(latData(20:156,1),-depthData,(nanmean(lonmPm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.1:0.02:0.1]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:200:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
ylabel('Depth (m)')
title([phases{k},' Com.zonal mean temperature in Pacific'])
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.phase',num2str(k),'.Pacific.zonalmean.tem.png'],'-dpng','-r600')
end
%% Atlantic
TemAtl = cat(1,Temcom(255:360,:,:,:),Temcom(1:20,:,:,:));
lonAtl = [lonData(255:360)-360;lonData(1:20)];
lonmA1 = permute(nanmean(TemAtl(39:end,20:101,:,:),1),[2 3 4 1]); % 67W-20E, -70S-10N
lonmA2 = permute(nanmean(TemAtl(16:end,102:109,:,:),1),[2 3 4 1]); % 90W-20E, 11N-18N 
lonmA3 = permute(nanmean(TemAtl(:,110:161,:,:),1),[2 3 4 1]); % 105W-20E, 19N-70N
lonmA = cat(1,lonmA1,lonmA2,lonmA3);
max(lonmA,[],'all')
min(lonmA,[],'all')
TemAtlm = cat(1,Temp(255:360,:,:,:),Temp(1:20,:,:,:)); % climotology
lonmA1m = permute(nanmean(TemAtlm(39:end,20:101,:,:),1),[2 3 4 1]); 
lonmA2m = permute(nanmean(TemAtlm(16:end,102:109,:,:),1),[2 3 4 1]); 
lonmA3m = permute(nanmean(TemAtlm(:,110:161,:,:),1),[2 3 4 1]); 
lonmAm = cat(1,lonmA1m,lonmA2m,lonmA3m);
%
close all;
for k = 1:3
Fig = figure('position',[100 100 800 400]);
contourf(latData(20:161,1),-depthData,lonmA(:,:,k)',[-1:0.01:1],'linestyle','none')
caxis([-0.1,0.1])
hold on
[ac ah] = contour(latData(20:161,1),-depthData,(nanmean(lonmAm,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.1:0.02:0.1]);
set(gca,'XLim',[-70,70]);
set(gca,'XTick',[-60:20:70]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:200:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
ylabel('Depth (m)')
title([phases{k},' Com.zonal mean temperature in Atlantic'])
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.phase',num2str(k),'.Atlantic.zonalmean.tem.png'],'-dpng','-r600')
end
%% meridional mean 
for k = 1:3
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[50:111],'40S-20N',k,phases,lonData,depthData);
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[20:40],'70S-50S',k,phases,lonData,depthData);
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[40:60],'50S-30S',k,phases,lonData,depthData);
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[60:116],'30S-25N',k,phases,lonData,depthData);
    meridionalmean(Temcom(:,:,:,k),Temp,latData,[116:141],'25N-50N',k,phases,lonData,depthData);
end
%% different depth
for k = 1:3
    diffdepth(Temcom(:,:,:,k),1,k,phases,lonData,latData,'0000',1.2,0.4)
    diffdepth(Temcom(:,:,:,k),7,k,phases,lonData,latData,'0050',1.2,0.4)
    diffdepth(Temcom(:,:,:,k),12,k,phases,lonData,latData,'0100',1.2,0.4)
    diffdepth(Temcom(:,:,:,k),17,k,phases,lonData,latData,'0200',1.2,0.3)
    diffdepth(Temcom(:,:,:,k),23,k,phases,lonData,latData,'0500',1.2,0.2)
    diffdepth(Temcom(:,:,:,k),29,k,phases,lonData,latData,'0800',1.2,0.1)
    diffdepth(Temcom(:,:,:,k),32,k,phases,lonData,latData,'1000',1.2,0.1)
    diffdepth(Temcom(:,:,:,k),37,k,phases,lonData,latData,'1500',1.2,0.05)
    diffdepth(Temcom(:,:,:,k),41,k,phases,lonData,latData,'2000',1.2,0.05)
end
%% 100-400m mean
Temcom_m = nanmean(Temcom(:,:,12:21,:),3);
for k = 1:3
    Fig = diffdepth(Temcom_m(:,:,k),1,k,phases,lonData,latData,'100-400',1.2,0.4)
    m_line([220 340 340 220 220 ],[-50 -50 30 30 -50 ],'linewidth',1.5,'color','k')
    print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.phase',num2str(k),'.diffdepth.100-400m.tem.png'],'-dpng','-r300')
end
%% zonal mean 140W-20W
lonm1 = permute(nanmean(Temcom(220:340,:,1:23,:),1),[2 3 4 1]);
lonm1m = permute(nanmean(Temp(220:340,:,1:23,:),1),[2 3 4 1]); % climotology
close all;
for k = 1:3
Fig = figure('position',[100 100 800 400]);
contourf(latData,-depthData(1:23),lonm1(:,:,k)',[-0.5:0.01:0.5],'linestyle','none')
caxis([-0.2,0.2])
hold on
line([-50 -50 30 30 -50],[-100 -400 -400 -100 -100],'linewidth',1.5,'color','k')
[ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.2:0.04:0.2]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([phases{k},' Com.zonal mean temperature 140W-20W'])
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.phase',num2str(k),'.140W_20W.zonalmean.tem.png'],'-dpng','-r600')
end
%% meridional mean 50S-30N
close all;
for k = 1:3
    names = '50S-30N';
    latm = latmean(Temcom(:,:,1:23,k),[40:121],latData);
    size(latm)
    latmc = latmean(Temp,[40:121],latData); % climotology
%     close all;
    Fig = figure('position',[100 100 800 400]);
    [c h] = contourf(lonData,-depthData(1:23),latm',[-0.5:0.01:0.5],'linestyle','none');
    caxis([-0.2,0.2])
    hold on
    line([220 220 340 340 220],[-100 -400 -400 -100 -100],'linewidth',1.5,'color','k')
    [ac ah] = contour(lonData,-depthData(1:23),(nanmean(latmc(:,1:23,:),3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    ch = colorbar;
    set(ch,'Ticks',[-0.2:0.04:0.2]);
    set(gca,'XLim',[0,360]);
    set(gca,'XTick',[0:60:360]);
    set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
    set(gca,'YLim',[-500,0]);
    set(gca,'YTick',[-500:100:0],'FontSize',14);
    set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
    ylabel('Depth (m)')
    title([phases{k},' Com.meridional mean temperature (',names,')'])
    print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.phase',num2str(k),'.meridionalmean',names,'.tem.png'],'-dpng','-r300')
end
%% AMO, Atlantic Nino
% AMO index
Ts1d = permute(Tempa(:,:,1,:),[4 1 2 3]);
for i = 1:360
    for j = 1:180
        Tsdetrend(:,i,j) = detrend(Ts1d(:,i,j)); % detrend
    end
end
%
Tsdetrendr = permute(Tsdetrend,[2 3 1]);
[AMO_zs AMO] = areamean(Tsdetrendr,290:360,91:161,latData); % 70W-0,0-70N
dT = 1; % interval
cf = 1/8;
AMOf = lanczosfilter(AMO,dT,cf,[],'low'); % 8 year filtered AMO index
%%
close all;
Fig = figure('position',[2700 100 800 400])
plot(AMO,'--','color',[0,46,166]/255,'LineWidth',1.5)
hold on
plot(AMOf,'color',[0,46,166]/255,'LineWidth',1.5)
plot(zeros(1,81),'k','LineWidth',1)
plot([5:37],ones(33)*mean(AMOf(5:37)),'color',[252,198,48]/255,'LineWidth',4)
plot([38:58],ones(21)*mean(AMOf(38:58)),'color',[252,198,48]/255,'LineWidth',4)
plot([59:75],ones(17)*mean(AMOf(59:75)),'color',[252,198,48]/255,'LineWidth',4)
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[-0.6,0.6]);
set(gca,'YTick',[-0.6:0.2:0.6]);
xlabel('YEAR','FontSize',12),ylabel('AMO index','FontSize',14);
legend('unfiltered','filtered','Location','north','Orientation','horizontal')
legend('boxoff');
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\AMOindex.png'],'-dpng','-r300')
%% Atlantic nino
[Anino_zs Anino] = areamean(Tsdetrendr,340:360,87:94,latData); % 20W-0,3S-3N
Aninof = lanczosfilter(Anino,dT,cf,[],'low'); % 8 year filtered AMO index
%
close all;
Fig = figure('position',[2700 100 800 400])
plot(Anino,'--','color',[0,46,166]/255,'LineWidth',1.5)
hold on
plot(Aninof,'color',[0,46,166]/255,'LineWidth',1.5)
plot(zeros(1,81),'k','LineWidth',1)
plot([5:37],ones(33)*mean(Aninof(5:37)),'color',[252,198,48]/255,'LineWidth',4)
plot([38:58],ones(21)*mean(Aninof(38:58)),'color',[252,198,48]/255,'LineWidth',4)
plot([59:75],ones(17)*mean(Aninof(59:75)),'color',[252,198,48]/255,'LineWidth',4)
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[-0.6,0.6]);
set(gca,'YTick',[-0.6:0.2:0.6]);
xlabel('YEAR','FontSize',12),ylabel('Atlantic Nino index','FontSize',14);
legend('unfiltered','filtered','Location','north','Orientation','horizontal')
legend('boxoff');
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\ANINOindex.png'],'-dpng','-r300')
%% 10m wind
filename = 'E:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\10mSI-MO\uwnd.10m.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
uwnd = ncread(filename,'uwnd');
uwnd = permute(mean(reshape(uwnd(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
ua = double(uwnd - mean(uwnd(:,:,42:71),3)); % remove climatology from 1981 to 2010
ua1 = permute(ua,[3 1 2]);
for i = 1:360
    for j = 1:180
        uad(:,i,j) = detrend(ua1(:,i,j)); % detrend
        uadf(:,i,j) = lanczosfilter(uad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
%
uadr = permute(uad,[2 3 1]);
uadfr = permute(uadf,[2 3 1]);
%%
filename = 'E:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\10mSI-MO\vwnd.10m.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
vwnd = ncread(filename,'vwnd');
vwnd = permute(mean(reshape(vwnd(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
va = double(vwnd - mean(vwnd(:,:,42:71),3)); % remove climatology from 1981 to 2010
va1 = permute(va,[3 1 2]);
for i = 1:360
    for j = 1:180
        vad(:,i,j) = detrend(va1(:,i,j)); % detrend
        vadf(:,i,j) = lanczosfilter(vad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
%
vadr = permute(vad,[2 3 1]);
vadfr = permute(vadf,[2 3 1]);
%% correlation
clear par11 h01
for i = 1:size(uadfr,1);
    for j = 1:size(uadfr,2);
        [par11(i,j),h01(i,j),t] = corr_eff(TPI_filtered(1:71),uadfr(i,j,5:end-1),0.05); % filtered
        [par12(i,j),h02(i,j),t] = corr_eff(TPI_filtered(1:71),vadfr(i,j,5:end-1),0.05); % filtered
    end
end
%%
close all;
Fig = figure('position',[100 100 800 400]);
m_proj('equidistant Cylindrical','lat',[-80 80],'long',[0 360]);
hold on
d = 8;
[x,y] = meshgrid(lonData,latData);
[hp,ht] = m_vec(2,x(1:d:end,1:d:end),y(1:d:end,1:d:end),par11(1:d:end,1:d:end)',par12(1:d:end,1:d:end)',...
    'headlength',8,'headwidth',10,'EdgeColor','b','linewidth',0.8);
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',9,'box','on','linestyle','none','fontsize',12);
m_arrow_scale2(hp,par11,par12,1,'label','1','position',[0.855 0.87 0.05 0.05])
print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.corr.uvwnd.png'],'-dpng','-r300')




function [ts] = lonmean(var,lats,latData)
% average along latitude (all longitudes)
% var is lon*lat*depth*time
% ts is lat*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,1)./nansum(cos(latData(lats)'/180*pi).*var2,1),[2 3 4 1]);
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
function contourVARra(var,val,ctr,lon1,lon2,lat1,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[lat1 lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
end
function meridionalmean(Tempadlf,Temp,latData,lats,names,k,phases,lonData,depthData)
    latm = latmean(Tempadlf,lats,latData);
    size(latm)
    latmc = latmean(Temp,lats,latData); % climotology
    max(latm,[],'all')
    min(latm,[],'all')
    %
%     close all;
    Fig = figure('position',[100 100 800 400]);
    contourf(lonData,-depthData,latm',[-0.5:0.01:0.5],'linestyle','none')
    caxis([-0.1,0.1])
    hold on
    [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    ch = colorbar;
    set(ch,'Ticks',[-1:0.02:1]);
    set(gca,'XLim',[0,360]);
    set(gca,'XTick',[0:60:360]);
    set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
    set(gca,'YLim',[-2000,0]);
    set(gca,'YTick',[-2000:200:0],'FontSize',14);
    set(gca,'YTickLabel',[2000:-200:0],'FontSize',14);
    ylabel('Depth (m)')
    title([phases{k},' Com.meridional mean temperature (',names,')'])
%     print(Fig,['E:\figures\IAP\Yearly\220724_8yrs_LinearDetrend\IPO.phase',num2str(k),'.meridionalmean',names,'.tem.png'],'-dpng','-r300')
end
function [Fig] = diffdepth(Temcom,nlevel,k,phases,lonData,latData,names,val,ticks)
    Fig = figure('position',[100 100 600 300]);
    contourVARra(Temcom(:,:,nlevel),[-val:ticks/10:val],1,0,360,-90,90,lonData,latData);
    caxis([-ticks,ticks]);
    ch = colorbar;
    set(ch,'Ticks',[-ticks:ticks/5:ticks]);
    % testdots(h05,'k',lonData,latData)
    title([phases{k},' Com.T',' ',names,'m']);
end
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end





    