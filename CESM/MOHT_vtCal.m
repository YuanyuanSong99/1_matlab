clc,clear;
addpath(genpath('/Volumes/CESM-post2/1_matlab/help'));
load("MatFile/lonData.mat");
load("MatFile/latData.mat");
load("MatFile/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
%  VT trend
VNT1 = read_nlev_1nc("/Volumes/CESM-post2/CESM-post/LE/VNT/VNTensmean_1920-2005.nc",'VNT');
VNT2 = read_nlev_1nc("/Volumes/CESM-post2/CESM-post/LE/VNT/VNTensmean_2006-2080.nc",'VNT');
VNT = cat(4,VNT1,VNT2(:,:,:,1:16));
size(VNT1)
size(VNT2)
clear VNT1 VNT2
depthstr = '0-700';
dweit = depthData(2:37)-depthData(1:36);
londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
lats = 1:60;
vntMMMsub = permute(nansum(londist'.*cat(3,VNT(:,:,1,:)*5,VNT(:,:,2:37,:).*permute(dweit,[3 2 1 4])),3),[1 2 4 3]);
%% VNT trend 2d
startyr = 1960;
endyr = 2020;
var = permute(vntMMMsub(:,lats,startyr-1919:endyr-1919),[3 1 2]);
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
ftsz = 14; ticks = 1;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd,[-ticks*5:ticks/10:ticks*5],ticks,lonData,latData(lats),12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-55,'linewidth',2,'color','k'); 
m_line(-70,[-55:1:-35],'linewidth',2,'color','k'); 
m_line(150,[-55:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'m/s · ^oC',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots(h0,[.4 .4 .4],lonData,latData(lats));
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240503_SO_warming/VNT',depthstr,'m_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%% MOHT calculate ensemble mean
inlev = 60; depthstr = '0-700m'; depth = depthData(1:inlev);
f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
dweit = depth(2:inlev)-depth(1:inlev-1);
clear levdist
levdist(1) = 4; levdist(2:inlev) = dweit; % z distance
OHCz = permute(sum((VNT(:,:,1:inlev,:)).*londist'.*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z
OHCzy = ro*cp*OHCz;
OHCzy_r = cat(1,OHCzy(181:360,:,:),OHCzy(1:180,:,:));
lonData_r = cat(1,lonData(181:360)-360,lonData(1:180));
MOHTp = permute(nansum(OHCzy_r(1:121,:,:),1),[2 3 4 1]); % Pac  unit:W
MOHTi = permute(nansum(OHCzy_r(211:301,:,:),1),[2 3 4 1]); % io 30E-120E unit:W
MOHTa = permute(nansum(OHCzy_r(121:191,:,:),1),[2 3 4 1]); % Atl 60W-10E unit:W
MOHTall = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
%%
startyr = 1960;
endyr = 2020;
lats = [20:60];
var1 = MOHTp(lats,startyr-1959:endyr-1959)';
var2 = MOHTi(lats,startyr-1959:endyr-1959)';
var3 = MOHTa(lats,startyr-1959:endyr-1959)';
x = [1:size(var1,1)]';
clear trd1 trd2 trd3
for i = 1:size(var1,2);
    par1=polyfit(x,var1(:,i),1); % regression parameters
    trd1(i) = par1(1);
    par2=polyfit(x,var2(:,i),1); % regression parameters
    trd2(i) = par2(1);
    par3=polyfit(x,var3(:,i),1); % regression parameters
    trd3(i) = par3(1);
end
%
close all;
figure(4)
plot(latData(lats),trd1,'r')
hold on
plot(latData(lats),trd2,'b')
plot(latData(lats),trd3,'k')
legend('pac','IO','atl')
%% MOHT calculate each member
inlev = 37; depthstr = '0-700m'; depthData = depthData(1:inlev);
f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:inlev)-depthData(1:inlev-1);
clear levdist
levdist(1) = 4; levdist(2:inlev) = dweit; % z distance
% historical 
datadir1='/Volumes/CESM-post2/CESM-post/LE/VNT/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*.B20TRC5CNBDRD*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
k=length(filelist1);
clear OHCz* MOHT*1
lonData_r = cat(1,lonData(180:360),lonData(1:179));
for s=2:11
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    TemVel = ncread(filename1,'VNT');
    OHCz = permute(sum((TemVel(:,:,:,:)).*londist'.*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = ro*cp*OHCz;
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
MOHTall1LE = nanmean(MOHTall1(:,:,2:10),3);
%
% close all
plot(latData,nanmean(MOHTall1LE,2),'k');
set(gca,'ylim',[-8,8]*10^10)
hold on
plot(latData,nanmean(MOHTp1LE,2),'r');
plot(latData,nanmean(MOHTi1LE,2),'g');
plot(latData,nanmean(MOHTa1LE,2),'b');
%% uninterpolated
clc,clear,close all;
addpath(genpath('/Volumes/CESM-post2/1_matlab/help'));
load("MatFile/depthData.mat");
nlev = length(depthData);
% MOHT calculate 
inlev = 11; depthstr = '0-2000m'; depthData = depthData(1:inlev);
splon = 170; lonstr = '150E';
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
dweit = depthData(2:inlev)-depthData(1:inlev-1);
clear levdist
levdist(1) = 4; levdist(2:inlev) = dweit; % z distance
% historical 
datadir1='/Volumes/CESM-post2/CESM-post/LE/VNTyr/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*.B20TRC5CNBDRD*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
TLONG = ncread([datadir1,filelist1(1).name],'TLONG');
ULAT = ncread([datadir1,filelist1(1).name],'ULAT');
londist = 2*pi*6.371393*10^6*cos(ULAT*pi/180)/360; % distance between 2 longitudes 
k=length(filelist1);
clear OHCz* MOHT*1
lonData_r = cat(1,TLONG(splon:end,:),TLONG(1:splon-1,:));
d300 = abs(lonData_r(:,1)-300); 
num300 = find(d300 == min(d300)); % 60W number
d20 = abs(lonData_r(:,1)-20); 
num20 = find(d20 == min(d20)); % 20E number
%
clear MOHT* midVT* OHC*
for s=1
    s 
    filename1=[datadir1,filelist1(s+1).name]; % start from the second
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    TemVel = ncread(filename1,'VNT').*londist; % the direct output VNT is V*T at each longitude
    midVT(:,:,1,:) = TemVel(:,:,1,:);
    midVT(:,:,2:inlev,:) = (TemVel(:,:,2:inlev,:) + TemVel(:,:,1:inlev-1,:))/2;
    OHCz = ro*cp*permute(sum(midVT.*londist.*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    MOHTall1(:,:,s) = permute(nansum(OHCz,1),[2 3 4 1]); % integrated along x
    OHCzy = OHCz;
    OHCzy_r = cat(1,OHCzy(splon:end,:,:),OHCzy(1:splon-1,:,:));
    MOHTp1(:,:,s) = permute(nansum(OHCzy_r(1:num300,:,:),1),[2 3 4 1]); % Pac 150E unit:W
    MOHTia1(:,:,s) = permute(nansum(OHCzy_r(num300+1:end,:,:),1),[2 3 4 1]); % IO+Atl 60W unit:W
    MOHTi1(:,:,s) = permute(nansum(OHCzy_r(num20:end,:,:),1),[2 3 4 1]); % io 20E-150E unit:W
    MOHTa1(:,:,s) = permute(nansum(OHCzy_r(num300:num20,:,:),1),[2 3 4 1]); % Atl 60W-20E unit:W
end

MOHTp1LE = nanmean(MOHTp1,3);
MOHTia1LE = nanmean(MOHTia1,3);
MOHTi1LE = nanmean(MOHTi1,3);
MOHTa1LE = nanmean(MOHTa1,3);
MOHTall1LE = nanmean(MOHTall1(:,:,:),3);
%
close all
plot(ULAT(1,:),nanmean(MOHTall1LE,2),'k');
% set(gca,'ylim',[-8,8]*10^10)
hold on
plot(ULAT(1,:),nanmean(MOHTp1LE,2),'r');
plot(ULAT(1,:),nanmean(MOHTi1LE,2),'g');
plot(ULAT(1,:),nanmean(MOHTa1LE,2),'b');






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
