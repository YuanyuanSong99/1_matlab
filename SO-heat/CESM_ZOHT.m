clc,clear,close all;
addpath(genpath('/Volumes/Togo/1_matlab/help'));
load("/Volumes/Togo4T/1_matlab/MatData/CESM/lonData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/latData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
% LE ensmean (EX)
% MOHT calculate 
inlev = 37; depthstr = '0-700m';
splon = 150; lonstr = '150E';
% UVEL
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:inlev)-depthData(1:inlev-1);
latdist = sw_dist([-90,-89],[0,0],'km')*10^3; % m
clear levdist
levdist(1) = 4; levdist(2:inlev) = dweit; % z distance
%% historical 
datadir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='/Volumes/CESM-post2/CESM-post/LE/UVEL/yearly/'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
datadir3='/Volumes/CESM-post2/CESM-post/LE/UISOP/yearly/'; %指定批量数据所在的文件夹
filelist3=dir([datadir3,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir3,filelist3(1).name]);
datadir4='/Volumes/CESM-post2/CESM-post/LE/USUBM/yearly/'; %指定批量数据所在的文件夹
filelist4=dir([datadir4,'*192001-200512.nc']); %指定批量数据的类型
ncdisp([datadir4,filelist4(1).name]);
% rcp85
datadir5='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist5=dir([datadir5,'*200601*.nc']); %指定批量数据的类型
datadir6='/Volumes/CESM-post2/CESM-post/LE/UVEL/yearly/'; %指定批量数据所在的文件夹
filelist6=dir([datadir6,'*200601*.nc']); %指定批量数据的类型
datadir7='/Volumes/CESM-post2/CESM-post/LE/UISOP/yearly/'; %指定批量数据所在的文件夹
filelist7=dir([datadir7,'*200601*.nc']); %指定批量数据的类型
datadir8='/Volumes/CESM-post2/CESM-post/LE/USUBM/yearly/'; %指定批量数据所在的文件夹
filelist8=dir([datadir8,'*200601*.nc']); %指定批量数据的类型k=length(filelist1)
clear OHCz* MOHT*1
lonData_r = cat(1,lonData(160:360),lonData(1:159)+360);
%%
lats = [1:90];
intlev = 37; % 700m
clear ZOHT
for s=1:40
    s 
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    UVEL1 = ncread(filename2,'UVEL')/100; % cm/s -> m/s
    filename3=[datadir3,filelist3(s).name];
    ncid=netcdf.open(filename3,'NC_NOWRITE');
    UISOP1 = ncread(filename3,'UISOP')/100; % cm/s -> m/s
    filename4=[datadir4,filelist4(s).name];
    ncid=netcdf.open(filename4,'NC_NOWRITE');
    USUBM1 = ncread(filename4,'USUBM')/100; % cm/s -> m/s
    Uvel1 = UVEL1(:,lats,1:intlev,:)+UISOP1(:,lats,1:intlev,:)+USUBM1(:,lats,1:intlev,:);
    clear UVEL1 UISOP1 USUBM1
    filename6=[datadir6,filelist6(s).name];
    ncid=netcdf.open(filename6,'NC_NOWRITE');
    UVEL2 = ncread(filename6,'UVEL')/100; % cm/s -> m/s
    filename7=[datadir7,filelist7(s).name];
    ncid=netcdf.open(filename7,'NC_NOWRITE');
    UISOP2 = ncread(filename7,'UISOP')/100; % cm/s -> m/s
    filename8=[datadir8,filelist8(s).name];
    ncid=netcdf.open(filename8,'NC_NOWRITE');
    USUBM2 = ncread(filename8,'USUBM')/100; % cm/s -> m/s
    Uvel2 = UVEL2(:,lats,1:intlev,1:16)+UISOP2(:,lats,1:intlev,1:16)+USUBM2(:,lats,1:intlev,1:16); % 2005-2020
    clear UVEL2 UISOP2 USUBM2
    Uvel = cat(4,Uvel1,Uvel2);
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp1 = ncread(filename1,'TEMP');
    filename5=[datadir5,filelist5(s).name];
    ncid=netcdf.open(filename5,'NC_NOWRITE');
    Temp2 = ncread(filename5,'TEMP');  
    Temp = cat(4,Temp1(:,lats,1:intlev,:),Temp2(:,lats,1:intlev,1:16));
    clear Temp1 Temp2
    % 0-700m MOHT
    TemUel = Temp.*Uvel; % 1920-2020
    OHCzy = latdist*permute(ro*cp*sum((TemUel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    ZOHT(:,:,s) = permute(nansum(OHCzy(:,30:55,:),2),[1 3 2]); % 60S-35S
end
clear Uvel Temp TemUel
%% ensmean
Temp1 = ncread('/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_1920-2005.nc','TEMP');
UVEL1 = ncread('/Volumes/CESM-post2/CESM-post/LE/UVEL/UVELensmean_1920-2005.nc','UVEL')/100; % cm/s -> m/s
UISOP1 = ncread('/Volumes/CESM-post2/CESM-post/LE/UISOP/UISOPensmean_1920-2005.nc','UISOP')/100; % cm/s -> m/s
USUBM1 = ncread('/Volumes/CESM-post2/CESM-post/LE/USUBM/USUBMensmean_1920-2005.nc','USUBM')/100; % cm/s -> m/s
Uvel1 = UVEL1 + UISOP1 + USUBM1;
clear UVEL1 UISOP1 USUBM1
%
Temp2 = ncread('/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_2006-2022.nc','TEMP');
UISOP2 = ncread('/Volumes/CESM-post2/CESM-post/LE/UISOP/UISOPensmean_2006-2080.nc','UISOP')/100; % cm/s -> m/s
USUBM2 = ncread('/Volumes/CESM-post2/CESM-post/LE/USUBM/USUBMensmean_2006-2080.nc','USUBM')/100; % cm/s -> m/s
UVEL2 = ncread('/Volumes/CESM-post2/CESM-post/LE/UVEL/UVELensmean_2006-2080.nc','UVEL')/100; % cm/s -> m/s
Uvel2 = UVEL2 + UISOP2 + USUBM2;
clear UVEL2 UISOP2 USUBM2
%%
intlev = 37; % 700m
Temp = cat(4,Temp1(:,:,1:intlev,:),Temp2(:,:,1:intlev,1:16));
Uvel = cat(4,Uvel1(:,:,1:intlev,:),Uvel2(:,:,1:intlev,1:16));
TemUel = Temp(:,1:90,1:inlev,:).*Uvel(:,1:90,1:inlev,:); % 1920-2020
OHCzy = latdist*permute(ro*cp*sum((TemUel(:,1:90,:,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z
%%
% OHCzy_r = cat(1,OHCzy(151:360,:,:),OHCzy(1:150,:,:));
% OHCrr = reshape(OHCzy_r,[360 3 30 86]);
% ZOHTe3 = squeeze(sum(OHCrr,2));
% [trd h0] = lntrend3d(ZOHTe3,.05);
% trd1 = sum()
%
lonData_r = cat(1,lonData(151:360,:,:),lonData(1:150,:,:));
ZOHTem = permute(nansum(OHCzy(:,1:55,:),2),[1 3 2]); % all

% trend
startyr = 1960;
endyr = 2020;
var = ZOHTem(:,startyr-1919:endyr-1919)';
x = [1:size(var,1)]';
clear trd
for i = 1:size(var,2);
    
        par=polyfit(x,var(:,i),1); % regression parameters
        trd(i) = par(1); 
    
end
% max(trd,[],'all')
% min(trd,[],'all')
%
close all;
Fig = figure('position',[100 100 800 400]);
p1 = plot(lonData,trd/10^12,'-','color','b','linewidth',2)
hold on
plot(lonData,zeros(360),'-','color','k','linewidth',1)
plot(150*ones(11),[0:4/10:4],'-','color','k','linewidth',2)
plot(290*ones(11),[0:4/10:4],'-','color','k','linewidth',2)
% text(40,3.8,'Pacific','fontsize',14);
% text(200,3.8,'Atlantic-Indian Ocean','fontsize',14);
set(gca,'XLim',[0,360]);set(gca,'XTick',[0:60:360]);
set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^oE','120^oW','60^oW'},'FontSize',14);
set(gca,'YLim',[0,4],'YTick',[0:1:4]);
ylabel('ZOHT (10^1^2 W yr^-^1)');xlabel('Longitude')
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240615_SO_warming/ZOHTensmean_trend.png'],'-dpng','-r300')
%% zoht anomaly time series 150E, 70W
ln1 = mean(ZOHTem(151,41:101),1); %1960-2020
ln2 = mean(ZOHTem(290,41:101),1); %1960-2020
p1 = ln1-mean(ln1(1:11)); % baseline 1960-1970
p2 = ln2-mean(ln2(1:11)); % baseline 1960-1970
bsln = '1960-1970';
yb1 = polyfit([1:61],ln1/10^15,1);
yb2 = polyfit([1:61],ln2/10^15,1);

close all;
ftsz = 20;
Fig = figure('position',[100 100 800 400]);
plot(ln1/10^15,'-','color',[.04 .35 .56],'linewidth',3)
hold on
plot(ln2/10^15,'-','color',[.64 .08 .18],'linewidth',3)

set(gca,'XLim',[1,61]);set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',ftsz);
set(gca,'YLim',[.8,1.6],'YTick',[.8:.2:1.6]);
ylabel('ZOHT (PW)');xlabel('Year')
legend(['150^oE'],['70^oW'],'location','north','fontsize',ftsz,'orientation','horizontal')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240809_SO_warming/ZOHTensmean_150E_70W.png'],'-dpng','-r300')
%%
close all;
ftsz = 20;
Fig = figure('position',[100 100 800 400]);
plot(p1/10^15,'-','color',[.04 .35 .56],'linewidth',3)
hold on
plot(p2/10^15,'-','color',[.64 .08 .18],'linewidth',3)
plot(zeros(61),'-','color','k','linewidth',1)
set(gca,'XLim',[1,61]);set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',ftsz);
set(gca,'YLim',[-.02,.12],'YTick',[-.02:.02:.12]);
ylabel('ZOHT anomaly (PW)');xlabel('Year')
text(40,-.01,['Baseline ',bsln],'fontsize',ftsz,'Color','k')
legend(['150^oE'],['70^oW'],'location','north','fontsize',ftsz,'orientation','horizontal')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240809_SO_warming/ZOHTbslnensmean_150E_70W.png'],'-dpng','-r300')
%% UNT trend 2d
startyr = 1960;
endyr = 2005;
lats = 1:90;
TU = permute(nansum(cat(3,TemUel(:,:,1,:)*5,TemUel(:,:,2:37,:).*permute(dweit,[3 2 1 4])),3),[1 2 4 3]);

var = permute(TU(:,lats,startyr-1919:endyr-1919),[3 1 2]);
x = [1:size(var,1)]';
clear trd
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
% h0 = trendtest(var,0.05); % t test trend
max(trd,[],'all')
min(trd,[],'all')
%% Tem plot
close all;
ftsz = 14; ticks = 1;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd,[-ticks*5:ticks/10:ticks*5],ticks,lonData,latData(lats),12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-40,'linewidth',2,'color','k'); 
m_line([0:1:360],-45,'linewidth',2,'color','k'); 
m_line([0:1:360],-50,'linewidth',2,'color','k'); 
m_line([0:1:360],-55,'linewidth',2,'color','k'); 
m_line(-70,[-55:1:-35],'linewidth',2,'color','k'); 
m_line(150,[-55:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'m/s · ^oC',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots(h0,[.4 .4 .4],lonData,latData(lats));
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240503_SO_warming/UNT',depthstr,'m_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%%
ZOHTda = ZOHT-nanmean(ZOHT,2)-ZOHTem;
ZOHTdar = permute(ZOHTda,[2 1 3]);
clear ZOHTadf
for s = 1:40
dT = 1; cf = 1/8;
for i = 1:360
    ZOHTadf(:,i,s) = lanczosfilter(ZOHTdar(:,i,s),dT,cf,[],'low'); % 8 year filtered
end
end
ZOHTadf(1:4,:,:) = []; 
ZOHTadf(end-3:end,:,:) = []; 
clear parZ hZ par87 h087 
for s = 1:40
index = TPIfz(1:78,s);
for i = 1:size(ZOHTadf,2);
    [parZ(i,s),hZ(i,s),t] = reg1_ttest(index,ZOHTadf(:,i,s),0.1,1);
end
end
%%
map = mean(parZ,2)/10^15;
mapr = cat(1,map(160:360),map(1:159));
close all;
Fig = figure('position',[100 100 800 400]);
p1 = plot(mapr,'-','color','b','linewidth',1.5)
hold on
plot(zeros(360),'-','color','k','linewidth',1)
plot(140*ones(7),[-0.024:0.01:0.04],'--','color','k','linewidth',1.5)
text(40,0.025,'Pacific','fontsize',14);
text(200,0.025,'Atlantic-Indian Ocean','fontsize',14);
set(gca,'XLim',[0,360]);set(gca,'XTick',[20:60:360]);
set(gca,'XTickLabel',{'180^oE','120^oW','60^oW','0^o','60^oE','120^oE'},'FontSize',14);
set(gca,'YLim',[-.024,.032],'YTick',[-.04:.01:.04]); set(gca,'YTick',[-.04:.01:.04],'FontSize',14);
ylabel('ZOHT (PW)');xlabel('Longitude')
% legend(gca,[p4,p2,p6],'Atlantic & Indian Ocean Total Eulerian mean','Atlantic & Indian Ocean Total Residual','Atlantic & Indian Ocean Ekman',...
%     'fontsize',12,'Location','north','Orientation','vertical','position',[0.28,0.63,0.18,0.3])
% legend('boxoff')
% ah = axes('position',get(gca,'position'),'Visible','off')
% legend(ah,[p5,p1,p3],'Pacific Ekman','Pacific Total Residual','Pacific Total Eulerian mean',...
%     'fontsize',12,'Location','north','Orientation','vertical','position',[0.2,0.1,0.18,0.3])
% legend('boxoff')
print(Fig,['/Volumes/Togo/figures/CESM/Yearly/LE/20240130_IPO_SO/ZOHT_Eulerian.png'],'-dpng','-r300')
%% dZOHT/dx
dz(1) = 0;
dz(2:360) = (map(2:end)-map(1:end-1))*10^15;
dzdx = dz/(londist(40));
dzdxr = cat(2,dzdx(160:360),dzdx(1:159));
close all;
Fig = figure('position',[100 100 800 400]);
nump = find(dzdxr > 0); numn = find(dzdxr < 0);
xx = [1:360];
h3 = bar(xx(nump),dzdxr(nump),'FaceColor',[.64 .08 .18],'EdgeColor',[.64 .08 .18],'barwidth',1)
hold on
bar(xx(numn),dzdxr(numn),'FaceColor',[.04 .35 .56],'EdgeColor',[.04 .35 .56],'barwidth',1)
plot(zeros(360))
plot(140*ones(15),[-8:1:6]*10^7,'-','color','k','linewidth',2)
text(40,4*10^7,'Pacific','fontsize',14);
text(200,4*10^7,'Atlantic-Indian Ocean','fontsize',14);
set(gca,'XLim',[0,360]);set(gca,'XTick',[20:60:360]);
set(gca,'XTickLabel',{'180^oE','120^oW','60^oW','0^o','60^oE','120^oE'},'FontSize',14);
% set(gca,'YLim',[-.0024,.0032],'YTick',[-.04:.01:.04]); set(gca,'YTick',[-.04:.01:.04],'FontSize',14);
ylabel('dZOHT/dx (W m^-^1)');xlabel('Longitude')
print(Fig,['/Volumes/Togo/figures/CESM/Yearly/LE/20240130_IPO_SO/dZOHTdx.png'],'-dpng','-r300')







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





