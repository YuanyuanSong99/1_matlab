% era5 forcing Data 
% clc,clear,close all;
addpath G:\1_matlab\help;
addpath G:\1_matlab\help\seawater\;
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
addpath G:\1_matlab\help;
filename1 = 'G:\data\LICOM\Tem_2000m.jra55.grid.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
lonData = ncread(filename1,'lon');
latData = ncread(filename1,'lat');
depthData = ncread(filename1,'lev');
% save('MatFile/lonData.mat','lonData');
% save('MatFile/latData.mat','latData');
% save('MatFile/depthData.mat','depthData');
Tem1 = ncread(filename1,'ts');
% filename2 = 'H:\CESM-post\LE\Temperature\rcp85_sub_ensmean.nc';
% ncid=netcdf.open(filename2,'NOWRITE');
% ncdisp(filename2);
% Tem2 = ncread(filename2,'TEMP');
% TemMMMall = cat(4,Tem1,Tem2(:,1:90,:,:));
% ------------------------- 0-700m ----------------------------------------
depthstr = '0-700';
dweit = -(depthData(2:21)-depthData(1:20));
lats = 1:60;
T700 = nanmean(Tem1(:,:,21:22,:),3);
Temsub = permute(nansum(cat(3,Tem1(:,:,1,:)*5,Tem1(:,:,2:21,:).*permute(dweit,[3 2 1 4]),T700*79),3)/700,[1 2 4 3]);

%% trend 2d
startyr = 1959;
endyr = 2020;
var = permute(Temsub(:,:,startyr-1958:endyr-1958),[3 1 2]);
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
%%
close all;
ftsz = 12; ticks = 0.1;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd*10,[-ticks*5:0.01:ticks*5],ticks,lonData,latData,12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
% m_line([0:1:360],-50,'linewidth',2,'color','k'); 
m_line(-60,[-90:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-90:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K decade^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots(h0,[.4 .4 .4],lonData,latData);
% print(Fig,['G:\figures\LICOM\Yearly\20231101\jra55.Tem',depthstr,'m_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% Time series
lats = [40:55]; % latitude
[spacz_long spac_long] = areamean(Tsub_r,1:120,lats,latData); 
[siaz_long sia_long] = areamean(Tsub_r,121:360,lats,latData); 
DI_rawLE = sia_long-spac_long;


%%
ncid=netcdf.open('H:\CESM-post\LE\Temperature\his_sst_ensmean.nc','NOWRITE');
ncdisp('H:\CESM-post\LE\Temperature\his_sst_ensmean.nc');
tsMMM = ncread('H:\CESM-post\LE\Temperature\his_sst_ensmean.nc','TEMP'); 
sstMMM = permute(tsMMM(:,:,1,:),[1 2 4 3]);
datadir='H:\CESM-post\LE\Temperature\his_sst_yr1x1\'; %指定批量数据所在的文件夹
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
%% MOHT calculate (0-700m) 
f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:end)-depthData(1:end-1);
levdist(1) = 4; levdist(2:37) = dweit; % z distance
datadir1='H:\CESM-post\LE\Temperature\his_sub_yr1x1\'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='H:\CESM-post\LE\VVEL\his_vvel_yr1x1\'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1);
clear OHCz*
%
for s=1:40
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VVEL')/100; % cm/s -> m/s
    % 0-700m MOHT
    TemVel = Temp.*Vvel(:,1:90,:,:);
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,1:37,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(180:300,:,:,:),OHCzy(301:360,:,:,:),OHCzy(1:179,:,:,:));
    MOHT71(:,:,s) = permute(nansum(OHCzy_r(1:121,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHT72(:,:,s) = permute(nansum(OHCzy_r(122:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHT7all(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
%%
close all;
Fig = figure('position',[100 100 600 300])
map1 = mean(mean(MOHT71,2),3)/10^15;
map2 = mean(mean(MOHT72,2),3)/10^15;
% map(85:96,:) = nan;
plot(latData(1:90),map1,'r','linewidth',1.5)
hold on
plot(latData(1:90),map2,'b','linewidth',1.5)
set(gca,'XLim',[-90,90],'YLim',[-2,2],'XGrid','on','YGrid','on');
set(gca,'XTick',[-80:20:80],'YTick',[-2:0.5:2]);
ylabel('MOHT (PW)');xlabel('Latitude')
% print(Fig,['G:\figures\CESM\Yearly\LE\MOHT_Tv_climotology.png'],'-dpng','-r300')  
% 跟实际差距较大
%% trend
MOHT1em = nanmean(MOHT71,3);
MOHT2em = nanmean(MOHT72,3);
startyr = 1940;
endyr = 2005;
var1 = permute(MOHT1em(:,startyr-1919:endyr-1919),[2 1]);
var2 = permute(MOHT2em(:,startyr-1919:endyr-1919),[2 1]);
x = [1:size(var1,1)]';
clear trd
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd1(i) = par1(1); 
        par2=polyfit(x,var2(:,i),1); % regression parameters
        trd2(i) = par2(1); 
end
h0 = trendtest(var1,0.05); % t test trend
max(trd,[],'all')
min(trd,[],'all')
%%
close all;
plot(trd1,'r')
hold on
plot(trd2,'b')
set(gca,'XLim',[0,90],'XGrid','on','YGrid','on');
set(gca,'XTicklabel',[-90:10:0]);



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
    load('G:/1_matlab/help/colorbar_mat/bl_re4.mat');
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