% IAP
% Trend 0-700m Tem150
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
datadir='/Volumes/Togo4T/data/IAP/temperature/Yearly/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
loniap = h5read(filetemp,'/lon');
latiap = h5read(filetemp,'/lat');
depthiap = h5read(filetemp,'/level');
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
% OHC
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthiap(2:27)-depthiap(1:26); % depth weight
OHCiap = cp*ro*permute(nansum(cat(3,Temp(:,:,1,:),Temp(:,:,2:27,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]); % J/m2
%------------------------- 700-2000m -----------------------------------------
dweit = depthiap(28:41)-depthiap(27:40); % depth weight
OHCiap_2000 = cp*ro*permute(nansum(cat(3,Temp(:,:,27,:),Temp(:,:,28:41,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]); % J/m2


%% EN4.2.1
clearvars -except cp ro OHC* lon* lat* depth*
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
datadir='/Volumes/Togo4T/data/EN4/annual/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'EN4*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonen4 = h5read(filetemp,'/lon');
laten4 = h5read(filetemp,'/lat');
depthen4 = h5read(filetemp,'/level');
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
T1en4 = Temp(:,:,1,1);
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthen4(2:24)-depthen4(1:23); % depth weight
Temp700 = Temp(:,:,24,:)+(Temp(:,:,25,:)-Temp(:,:,24,:))*(700-depthen4(24))/(depthen4(25)-depthen4(24));
OHCen4 = cp*ro*permute(nansum(cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:24,:).*(permute(dweit,[3 2 1])),Temp700*(700-depthen4(24))),3),[1 2 4 3]); % J/m2
% ------------------------- 700-2000m -----------------------------------------
dweit = depthen4(26:31)-depthen4(25:30); % depth weight
Temp2000 = Temp(:,:,30,:)+(Temp(:,:,31,:)-Temp(:,:,30,:))*(700-depthen4(30))/(depthen4(31)-depthen4(30));
OHCen4_2000 = cp*ro*permute(nansum(cat(3,Temp700*(depthen4(25)-700),...
    Temp(:,:,25:30,:).*(permute(dweit,[3 2 1])),Temp2000*(2000-depthen4(30))),3),[1 2 4 3]); % J/m2

%% Ishii
clearvars -except cp ro OHC* lon* lat* depth* T1en4
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
datadir='/Volumes/Togo4T/data/Ishii/v2017_annual/temp/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'Ishii*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonish = h5read(filetemp,'/lon');
latish = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
Temps = permute(flipud(Temp),[3 2 1 4]);
depthish = flipud(depthData);
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthish(2:16)-depthish(1:15); % depth weight
OHCish = cp*ro*permute(nansum(cat(3,Temps(:,:,1,:),Temps(:,:,2:16,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]); % J/m2
% ------------------------- 700-2000m -----------------------------------------
dweit = depthish(17:26)-depthish(16:25); % depth weight
OHCish_2000 = cp*ro*permute(nansum(cat(3,Temps(:,:,16,:),Temps(:,:,17:26,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]); % J/m2
clearvars -except cp ro OHC* lon* lat* depth* T1en4
%% combine 700m
startyr = 1960; endyr = 2020;
OHC1 = OHCiap(:,8:end,startyr-1939:endyr-1939);
OHC2 = OHCen4(:,:,startyr-1949:endyr-1949); % lonen4 laten4
OHC3 = OHCish(:,8:end,startyr-1954:endyr-1954);
OHCall = (OHC1+OHC2+OHC3)/3;
var = permute(OHCall,[3 1 2]);
x = [1:size(var,1)]';
clear trd 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.01); % t test trend
h0(find(isnan(T1en4(:,:,1,1)) == 1)) = nan;
max(trd,[],'all')
min(trd,[],'all')
%% SO plot
close all;
mapr = cat(1,map,map(1,:));
lonr = [lonen4;lonen4(1)];
ftsz = 20; ticks = 1.5;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(mapr,[-ticks*5:ticks/10:ticks*5],ticks,lonr,laten4,ftsz)
m_coast('patch','w');
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2',5.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots0(h0(1:2:end,1:2:end),[.4 .4 .4],lonen4(1:2:end,1:2:end),laten4(1:2:end,1:2:end));
%print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240705/OHC',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.3.png'],'-dpng','-r300')
%% combine 700-2000m
depthstr = '700-2000m'
startyr = 1960; endyr = 2020;
OHC1 = OHCiap_2000(:,8:end,startyr-1939:endyr-1939);
OHC2 = OHCen4_2000(:,:,startyr-1949:endyr-1949); % lonen4 laten4
OHC3 = OHCish_2000(:,8:end,startyr-1954:endyr-1954);
OHCall = (OHC1+OHC2+OHC3)/3;
var = permute(OHCall,[3 1 2]);
x = [1:size(var,1)]';
clear trd 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.01); % t test trend
h0(find(isnan(T1en4(:,:,1,1)) == 1)) = nan;
max(trd,[],'all')
min(trd,[],'all')
%% SO plot
close all;
map = trd/(60*60*24*365); % J/m2 /yr -> W/m2
mapr = cat(1,map,map(1,:));
lonr = [lonen4;lonen4(1)];
ftsz = 20; ticks = 1.5;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(mapr,[-ticks*5:ticks/10:ticks*5],ticks,lonr,laten4,ftsz)
m_coast('patch','w');
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2',5.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots0(h0(1:2:end,1:2:end),[.4 .4 .4],lonen4(1:2:end,1:2:end),laten4(1:2:end,1:2:end));
%print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240705/OHC',depthstr,'_trend_',num2str(startyr),'_',num2str(endyr),'.3.png'],'-dpng','-r300')
%% meridional mean
% IAP
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
datadir='/Volumes/Togo4T/data/IAP/temperature/Yearly/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
loniap = h5read(filetemp,'/lon');
latiap = h5read(filetemp,'/lat');
depthiap = h5read(filetemp,'/level');
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Tempiap(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
% EN4.2.1
clearvars -except Temp* lon* lat* depth*
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
datadir='/Volumes/Togo4T/data/EN4/annual/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'EN4*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonen4 = h5read(filetemp,'/lon');
laten4 = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Teme0(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
Teme1 = permute(Teme0,[3 1 2 4]);
Teme2 = interp1(depthData,Teme1,depthiap);
Tempen4 = permute(Teme2,[2 3 1 4]);
% Ishii
clearvars -except Temp* lon* lat* depth*
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
datadir='/Volumes/Togo4T/data/Ishii/v2017_annual/temp/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'Ishii*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonish = h5read(filetemp,'/lon');
latish = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temi0(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
Temi1 = interp1(depthData,Temi0,depthiap); 
Tempish = permute(Temi1,[3 2 1 4]);
% clearvars -except Temp* lon* lat* depth*

%% combine
startyr = 1960; endyr = 2020;
Temp1 = Tempiap(:,8:end,:,startyr-1939:endyr-1939);
Temp2 = Tempen4(:,:,:,startyr-1949:endyr-1949);
Temp3 = Tempish(:,8:end,:,startyr-1954:endyr-1954);
Tempall = (Temp1+Temp2+Temp3)/3;
% meridional mean
lats = 29:49; % 55S-35S
latstr = '35S_55S';
[Tempa_mm] = latmean(Tempall,lats,laten4);
% ts is lon*depth*time
startyr = 1960;
endyr = 2020;
yrstr = '1960-2020';
var = permute(Tempa_mm(:,:,startyr-1959:endyr-1959),[3 1 2]);
x = [1:size(var,1)]';
clear trd
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.01); % t test trend
max(trd,[],'all')
min(trd,[],'all')
%% 150 longitude
close all;
ftsz = 20;
ticks = 0.1;
map = trd*10;
map_r = cat(1,map(150:360,:),map(1:150,:));
h0_r = cat(1,h0(150:360,:),h0(1:150,:));
lonData_r = [lonen4(150:360);lonen4(1:150)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthiap,map_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([290 290],[0 -2000],'linewidth',2,'color','k')
line([150 509],[-700 -700],'linewidth',1.5,'color','k','linestyle','--')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = [];
colormap(flipud(bl_re4));
ch = colorbar('location','eastoutside');
[i,j] = find(h0_r == 0); % significance test
plot(lonData_r(i(1:3:end)),-depthiap(j(1:3:end)),'.','MarkerSize',5,'color',[.4 .4 .4])
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','K decade^-^1','Units','normalized','position',[6 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[150,360+150]);
set(gca,'XTick',[150:60:360+150]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',ftsz);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:500:0],'FontSize',ftsz);
set(gca,'YTickLabel',[2000:-500:0],'FontSize',ftsz);
ylabel('Depth (m)')
print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240705/meridionalmean_',latstr,'_trend_',yrstr,'_150E.png'],'-dpng','-r300')
%% zonal mean 

datadir='/Volumes/Togo4T/data/Ishii/v2017_annual/temp/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'Ishii*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonish = h5read(filetemp,'/lon');
latish = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
k=length(filelist);
clear Temp
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
Temp=permute(Temp,[3 2 1 4]);
Temp_r = cat(1,Temp(150:360,:,:,:),Temp(1:149,:,:,:));
Temp_Pac = squeeze(nanmean(Temp_r(1:151,:,:,:),1));
Temp_AIO = squeeze(nanmean(Temp_r(152:360,:,:,:),1));

var1 = permute(Temp_Pac,[3 1 2]);
var2 = permute(Temp_AIO,[3 1 2]);
x = [1:size(var1,1)]';
clear trd1 trd2
for i = 1:size(var1,2);
    for j = 1:size(var1,3);
        par1=polyfit(x,var1(:,i,j),1); % regression parameters
        trd1(i,j) = par1(1); 

        par2=polyfit(x,var2(:,i,j),1); % regression parameters
        trd2(i,j) = par2(1); 
    end
end
%
close all
figure(1)
contourf(latData(1:60),-depthData(1:27),trd1(1:60,1:27)',[-0.02:0.002:0.02])
    caxis([-0.01,0.01]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    colorbar

figure(2)
contourf(latData(1:60),-depthData(1:27),trd2(1:60,1:27)',[-0.02:0.002:0.02])
    caxis([-0.01,0.01]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    colorbar


function [ts] = latmean(var,lats,latData)
% average along longitude (all latitudes)
% var is lon*lat*depth*time
% ts is lon*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,2)./nansum(cos(latData(lats)'/180*pi).*var2,2),[1 3 4 2]);
end
function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','none','xaxisloc','top');
end
function contourfSPolar(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',80,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
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
        elseif abs(T) < tinv(1-alp/2,d1-2);
            h0(i,j) = 0;
        end
    end
end
end
function testdots1(h,clor,lonData,latData)
% used to dot significant area
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',2,'color',clor); % F-test dots
end
function testdots0(h,clor,lonData,latData)
% used to dot insignificant area
    dots = mod(find(h == 0),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 0)/length(lonData))),'.','markersize',5,'color',clor); % F-test dots
end
