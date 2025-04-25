clc,clear,close all;
addpath(genpath('D:\1_matlab\help'));
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
% LE ensmean (EX)
addpath D:\1_matlab\help;
% MOHT calculate 
inlev = 60; depthstr = '0-5000m'; depthData = depthData(1:inlev);
splon = 150; lonstr = '150E';

f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:inlev)-depthData(1:inlev-1);
clear levdist
levdist(1) = 4; levdist(2:inlev) = dweit; % z distance
% historical 
datadir1='E:\CESM-post\LE\TEMP\yearly\'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*.B20TRC5CNBDRD*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='E:\CESM-post\LE\VVEL\historical\192001-200512\'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1);
clear OHCz* MOHT*1
lonData_r = cat(1,lonData(180:360),lonData(1:179));
for s=1
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VVEL')/100; % cm/s -> m/s
    % 0-700m MOHT
    TemVel = Temp(:,:,1:inlev,1:86).*Vvel(:,:,1:inlev,1:86); % 1960-2005
    OHCz = permute(ro*cp*sum((TemVel(:,:,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(splon:360,:,:),OHCzy(1:splon-1,:,:));
    MOHTp1(:,:,s) = permute(nansum(OHCzy_r(1:300-splon,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia1(:,:,s) = permute(nansum(OHCzy_r(301-splon:end,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTi1(:,:,s) = permute(nansum(OHCzy_r(382-splon:end,:,:),1),[2 3 4 1]); % io 20E-150E unit:W
    MOHTa1(:,:,s) = permute(nansum(OHCzy_r(301-splon:381-splon,:,:),1),[2 3 4 1]); % Atl 60W-20E unit:W
    MOHTall1(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
MOHTp1LE = nanmean(MOHTp1,3);
MOHTia1LE = nanmean(MOHTia1,3);
MOHTi1LE = nanmean(MOHTi1,3);
MOHTa1LE = nanmean(MOHTa1,3);
MOHTall1LE = nanmean(MOHTall1,3);
%%
close all
plot(latData,-nanmean(MOHTall1LE,2));
set(gca,'ylim',[-2,2]*10^15)
%% rcp85 2006-2022
datadir1='E:\CESM-post\LE\TEMP\rcp85\200601-208012\'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='E:\CESM-post\LE\VVEL\rcp85\200601-202212\'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1);
clear OHCz* 
%
for s=1:39
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VVEL')/100; % cm/s -> m/s
    % 0-2000m MOHT
    TemVel = Temp(:,1:90,1:inlev,1:15).*Vvel(:,1:90,1:inlev,1:15); % 2006-2020
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(splon:360,:,:,:),OHCzy(1:splon-1,:,:,:));
    MOHTp2(:,:,s) = permute(nansum(OHCzy_r(1:300-splon,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia2(:,:,s) = permute(nansum(OHCzy_r(301-splon:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTi2(:,:,s) = permute(nansum(OHCzy_r(382-splon:end,:,:,:),1),[2 3 4 1]); % io 20E-150E unit:W
    MOHTa2(:,:,s) = permute(nansum(OHCzy_r(301-splon:381-splon,:,:,:),1),[2 3 4 1]); % Atl 60W-20E unit:W
    MOHTall2(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
MOHTp2LE = nanmean(MOHTp2,3);
MOHTia2LE = nanmean(MOHTia2,3);
MOHTi2LE = nanmean(MOHTi2,3);
MOHTa2LE = nanmean(MOHTa2,3);
MOHTall2LE = nanmean(MOHTall2,3);
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
h0 = trendtest(var1,0.05); % t test trend
max(trd1,[],'all')
min(trd1,[],'all')
%
close all;
plot(trd1,'r')
hold on
plot(trd2,'b')
%% MOHT 

map1 = MOHTcp; map2 = trd1; clor = [.04 .35 .56];
regstr1 = 'Pacific'; regstr2 = 'Pac'; 
% 
map1 = MOHTcia; map2 = trd2; clor = [.64 .08 .18];
regstr1 = 'Atlantic-Indian Ocean'; regstr2 = 'AtlIO';
% % 
map1 = MOHTci; map2 = trd3; clor = [.64 .08 .18];
regstr1 = 'Indian Ocean'; regstr2 = 'IO';
% % 
map1 = MOHTca; map2 = trd4; clor = [.64 .08 .18];
regstr1 = 'Atlantic'; regstr2 = 'Atl';

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
plot(-50*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')

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
print(Fig,['D:\figures\CESM\Yearly\LE\0_LEM\20240120_SO_Warming\MOHT_',regstr2,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% MOHT/dlon 
dlon = 300-splon;
map1 = MOHTcp/dlon; map2 = trd1/dlon; clor = [.04 .35 .56];
regstr1 = 'Pacific'; regstr2 = 'Pac'; 
% % 
dlon = 360-(300-splon);
map1 = MOHTcia/dlon; map2 = trd2/dlon; clor = [.64 .08 .18];
regstr1 = 'Atlantic-Indian Ocean'; regstr2 = 'AtlIO';
% % % % % 
dlon = splon-20;
map1 = MOHTci/dlon; map2 = trd3/dlon; clor = [.64 .08 .18];
regstr1 = 'Indian Ocean'; regstr2 = 'IO';
% % % % %  
dlon = 360-300;
map1 = MOHTca/dlon; map2 = trd4/dlon; clor = [.64 .08 .18];
regstr1 = 'Atlantic'; regstr2 = 'Atl';

close all;
Fig = figure('position',[700 100 800 400]);
% ymax = 1.5; ymin = -1.5; yint = 0.5; % 100m
% ymax = 2.5; ymin = -2.5; yint = 0.5; % 300m
ymax = 3; ymin = -3; yint = 1; % 700m
p2 = plot(latData(10:60),map2(10:60)*10/10^11,'-','color',clor,'LineWidth',3);
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);
hold on
set(gca,'XLim',[-80,-30],'xtick',[-80:10:-30],'fontsize',12,'XGrid','on','YGrid','on');
set(gca,'XTicklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS'});
ylabel('Trend  (10^-^4 PW decade^-^1 Longitude^-^1)','FontSize',14);
plot(-35*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-50*ones(11),[ymin:(ymax-ymin)/10:ymax],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
p1 = plot(latData(10:60),map1(10:60)/10^13,'-','color','k','LineWidth',3);
hold on
% p4 = plot(climap1,'-','color','k','LineWidth',3);
plot(latData(10:60),zeros(1,51),'k','LineWidth',1)
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('Latitude','FontSize',12),ylabel('Climo (10^-^2 PW Longitude^-^1)','FontSize',14);
legend([p1, p2],['Climo (',regstr1,')'],['Trend (',regstr1,')'],'Location','northwest')
legend('boxoff')

yyaxis left
set(gca,'YLim',[ymin,ymax],'YTick',[ymin:yint:ymax],'YColor',clor);

print(Fig,['D:\figures\CESM\Yearly\LE\0_LEM\20240120_SO_Warming\MOHTdlon_',regstr2,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% dMOHT/dy
dlon = 300-splon; regstr = 'Pac';
dMdy = dMOHTdy(MOHTp);
dMdycli = dMOHTdy(MOHTcp);
ymin1 = -39; ymax1 = 39;  ymin2 = -39; ymax2 = 39;
ymin3 = -9; ymax3 = 9;  ymin4 = -3; ymax4 = 3;
% 
dlon = 360-(300-splon); regstr = 'AtlIO';
dMdy = dMOHTdy(MOHTia);
dMdycli = dMOHTdy(MOHTcia);
ymin1 = -39; ymax1 = 39;  ymin2 = -39; ymax2 = 39;
ymin3 = -9; ymax3 = 9;  ymin4 = -3; ymax4 = 3;
% % % % % 
dlon = splon-20; regstr = 'IO';
dMdy = dMOHTdy(MOHTi);
dMdycli = dMOHTdy(MOHTci);
ymin1 = -126; ymax1 = 126;  ymin2 = -138; ymax2 = 138;
ymin3 = -9; ymax3 = 9;  ymin4 = -12; ymax4 = 12;
% % % % %
dlon = 360-300; regstr = 'Atl';
dMdy = dMOHTdy(MOHTa);
dMdycli = dMOHTdy(MOHTca);
ymin1 = -126; ymax1 = 126;  ymin2 = -138; ymax2 = 138;
ymin3 = -9; ymax3 = 9;  ymin4 = -12; ymax4 = 12;

startyr = 1960; endyr = 2020;
var1 = permute(dMdy(:,startyr-1959:endyr-1959),[2 1]);
x = [1:size(var1,1)]';
clear trd
for i = 1:size(var1,2);
        par1=polyfit(x,var1(:,i),1); % regression parameters
        trd(i) = par1(1); 
end
% dMOHT/dy/dlon
close all;
Fig = figure('position',[100 100 800 400])
map1 = trd(1:60)*10/10^4/dlon;
nump = find(map1 > 0); numn = find(map1 < 0);
h3 = bar(latData(nump),map1(nump),'FaceColor',[.64 .08 .18],'EdgeColor',[.64 .08 .18],'barwidth',1)
hold on
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:ymax1/3:ymax1]);
bar(latData(numn),map1(numn),'FaceColor',[.04 .35 .56],'EdgeColor',[.04 .35 .56],'barwidth',1)
ylabel('Trend (10^4 W m^-^1 Longitude^-^1 decade^-^1)')
plot(-35*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')
plot(-50*ones(11),[ymin1:(ymax1-ymin1)/10:ymax1],'linewidth',1.5,'Color','k',LineStyle='--')

yyaxis right
map2 = dMdycli(1:60)/10^6/dlon;
numpc = find(map2 > 0); numnc = find(map2 < 0);
h1 = bar(latData(numpc),map2(numpc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
h2 = bar(latData(numnc),map2(numnc),'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'facealpha',.3,'LineWidth',1);
plot(latData(1:60),zeros(1,60),'k','LineWidth',1)
set(gca,'XLim',[-80,-30],'XTick',[-80:10:-30],'XTickLabel',[-80:10:30],'FontSize',14,'XGrid','on','YGrid','on');
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:ymax2/3:ymax2],'YColor',[.5 .5 .5]);
xlabel('Latitude','FontSize',14),ylabel('Climo (10^6 W m^-^1 Longitude^-^1)','FontSize',14);
print(Fig,['D:\figures\CESM\Yearly\LE\0_LEM\20240120_SO_Warming\MOHTdydlon_',regstr,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
% dMOHT/dy
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
print(Fig,['D:\figures\CESM\Yearly\LE\0_LEM\20240120_SO_Warming\MOHTdy_',regstr,'_',lonstr,'_',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')


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
    load('D:/1_matlab/help/colorbar_mat/bl_re4.mat');
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
