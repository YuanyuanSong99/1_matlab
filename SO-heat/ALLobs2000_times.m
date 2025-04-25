% 700 - 2000 m
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
% IAP data
datadir='/Volumes/Togo4T/data/IAP/temperature/Yearly/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% anomaly 
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
Tempa = double(Temp - nanmean(Temp(:,:,:,42:71),4)); % remove climatology from 1981-2010 
sst_iap = permute(Tempa(:,:,1,:),[1 2 4 3]);
%% tem
dweit = depthData(27:41)-depthData(26:40); % depth weight
lats = [35:55]; % 55S-35S
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
Temp = double(Temp);
Tempz = Temp(:,:,27:41,:).*(permute(dweit,[3 2 1]));
VVz = VV(:,:,27:41,:).*(permute(dweit,[3 2 1]));
VVz(find(isnan(Tempz(:,:,:,30)) == 1)) = nan;
Tempzy = Tempz*111*1000; % 111 km / latitude
VVzy = VVz*111*1000; % 111 km / latitude
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(150:360,:,:,:),Tempzyx(1:149,:,:,:)); 
VVzyx(find(isnan(Tempzyx(:,:,1,1)) == 1)) = nan;
Vzyxsub_r_0 = cat(1,VVzyx(150:360,:,:),VVzyx(1:149,:,:)); 
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,:,:),1),2),3));
spacV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,lats,:),1),2),3)); % Square
spac0V = spac_0./spacV;
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,:,:),1),2),3));
siaV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,lats,:),1),2),3)); % Square
sia0V = sia_0./siaV; % /m2
sia_IAP = sia_0;
spac_IAP = spac_0; 
spacV_IAP = spacV; 
siaV_IAP = siaV; 

%% EN4.2.1 data
clearvars -except sia_IAP spac_IAP siaV_IAP spacV_IAP
datadir='/Volumes/Togo4T/data/EN4/annual/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'EN4*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% anomaly 
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
dweit = depthData(26:31)-depthData(25:30); % depth weight
lats = [29:49]; % 55S-35S
Temp = double(Temp);
Temp700 = Temp(:,:,24,:)+(Temp(:,:,25,:)-Temp(:,:,24,:))*(700-depthData(24))/(depthData(25)-depthData(24));
Temp2000 = Temp(:,:,30,:)+(Temp(:,:,31,:)-Temp(:,:,30,:))*(700-depthData(30))/(depthData(31)-depthData(30));
Tempz = cat(3,Temp700*(700-depthData(24)),Temp(:,:,25:30,:).*(permute(dweit,[3 2 1])),Temp2000*(2000-depthData(30)));
Tempzy = Tempz*111*1000; % 111 km / latitude
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(151:360,:,:,:),Tempzyx(1:150,:,:,:)); 
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,:,:),1),2),3));
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,:,:),1),2),3));
sia_EN4 = sia_0;
spac_EN4 = spac_0;
% spac0S_EN4 = spac0S;
% sia0S_EN4 = sia0S;
%% Ishii Data
clearvars -except sia_IAP spac_IAP siaV_IAP spacV_IAP sia_EN4 spac_EN4
datadir='/Volumes/Togo4T/data/Ishii/v2017_annual/temp/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'Ishii*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% anomaly 
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end

Temps = permute(flipud(Temp),[3 2 1 4]);
depthData = flipud(depthData);
%% Tem
lats = [35:55]; % 55S-35S
dweit = depthData(16:26)-depthData(15:25); % depth weight
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
Temps = double(Temps);
Tempz = Temps(:,:,16:26,:).*(permute(dweit,[3 2 1]));
VVz = VV(:,:,16:26).*(permute(dweit,[3 2 1]));
Tempzy = Tempz*111*1000; % 111 km / latitude
VVzy = VVz*111*1000;
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(151:360,:,:,:),Tempzyx(1:150,:,:,:)); 
VVzyx(find(isnan(Tempzyx(:,:,:,1)) == 1)) = nan;
Vzyxsub_r_0 = cat(1,VVzyx(151:360,:,:),VVzyx(1:150,:,:)); 
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,:,:),1),2),3));
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,:,:),1),2),3));
sia_Ishii = sia_0;
spac_Ishii = spac_0;
%% CESM1 LE
clearvars -except sia_IAP spac_IAP siaV_IAP spacV_IAP sia_EN4 spac_EN4 ...
     sia_Ishii spac_Ishii 
load("/Volumes/Togo4T/1_matlab/MatData/CESM/lonData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/latData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
% LE ensmean (EX)
addpath /Volumes/Togo4T/1_matlab/help;
filename1 = '/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_1920-2005.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
TemMMM1 = ncread(filename1,'TEMP');
depthData = ncread(filename1,'z_t')/100; % cm -> m
filename2 = '/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_2006-2100.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
TemMMM2 = ncread(filename2,'TEMP');
TemMMMall = cat(4,TemMMM1,TemMMM2(:,:,1:47,:));
Temp = TemMMMall;
clear TemMMM1 TemMMM2
%% original OHC time series
depthstr = '0-700m';
dweit = depthData(37:46)-depthData(36:45);
lats = [35:55]; % 55S-35S
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
Tempz = Temp(:,:,37:46,:).*(permute(dweit,[3 2 1]));
Tempzy = Tempz*111*1000; % 111 km / latitude
clear dx Tempzyx 
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(152:360,:,:,:),Tempzyx(1:151,:,:,:)); 
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,:,:),1),2),3));
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,:,:),1),2),3));
sia_LE = sia_0;
spac_LE = spac_0;
%% CESM each member
datadir='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir,'*1920*']); %指定批量数据的类型
filelist2=dir([datadir,'*200601-210012*']); %指定批量数据的类型
load("/Volumes/Togo4T/1_matlab/MatData/CESM/lonData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/latData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/depthData.mat");
for j = 1:length(latData);
    dx(1,j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
end
%
clear CESM*
for s = 1:40
    s
    filename1 = [datadir,filelist1(s).name];
    filename2 = [datadir,filelist2(s).name];
    [CESMpac(:,s) CESMia(:,s)] = caleach(filename1,filename2,depthData,dx);
end

%%
stdpac = std(CESMpac',1);
mpac = mean(CESMpac,2)';
stdia = std(CESMia',1);
mia = mean(CESMia,2)';

close all
figure(1)
hold on
for s = 1:40
    plot(CESMpac(:,s),'b')
    plot(CESMia(:,s),'r')
end
plot(mpac,'k','LineWidth',2)
plot(mpac+stdpac,'k','LineWidth',2)
plot(mia,'k','LineWidth',2)

mspac1 = mpac+stdpac;
mspac2 = mpac-stdpac;

%% plot Pac 1960-2020 ZJ or GJ/m2
close all;
depthstr = '700-2000m'
bsln = '1960-1970'; % baseline
% ln1 = (spac_IAP-mean(spac_IAP(51:81)))/10^21;
% ln2 = (spac_EN4-mean(spac_EN4(41:71)))/10^21;
% ln3 = (spac_Ishii-mean(spac_Ishii(36:66)))/10^21;
% ln4 = (spac_LE-mean(spac_LE(71:101)))/10^21;
ln1 = (spac_IAP(21:81)-mean(spac_IAP(21:31)))/spacV_IAP; % 1960-2020
ln2 = (spac_EN4(11:71)-mean(spac_EN4(11:21)))/spacV_IAP;
ln3 = (spac_Ishii(6:66)-mean(spac_Ishii(6:16)))/spacV_IAP;
ln4 = (spac_LE(41:end-1)-mean(spac_LE(41:51)))/spacV_IAP;
ln5 = (spac_LE(41:end-1)-mean(spac_LE(41:51))+stdpac(1:end-1)')/spacV_IAP;
ln6 = (spac_LE(41:end-1)-mean(spac_LE(41:51))-stdpac(1:end-1)')/spacV_IAP;
%
close all
ftsz = 20;
Fig = figure('position',[700 50 700 400]);
% 2100
% x1 = [1:1.7:103];
% x2 = [104:1:182];
% x = [x1,x2];
% h = fill([x,fliplr(x)],[ln5',fliplr(ln6')],[.9 .9 .9])
% 2020
x1 = [1:61];
x = x1;
h = fill([x,fliplr(x)],[ln5(1:61)',fliplr(ln6(1:61)')],[.9 .9 .9])

hold on
p1 = plot(x1,ln1,'-','color',[0.78,0.58,0.11],'LineWidth',3);
p2 = plot(x1,ln2,'-','color',[.04 .35 .56],'LineWidth',3);
p3 = plot(x1,ln3,'-','color',[.64 .08 .18],'LineWidth',3);
p4 = plot(x,ln4(1:length(x)),'-','color',[.6 .6 .6],'LineWidth',3);
plot(x,zeros(length(x)),'k-','LineWidth',1.5)
ylabel('Temperature  (K)','FontSize',ftsz);
set(gca,'XGrid','on','YGrid','on');
xlabel('Year','FontSize',ftsz);
% 2100
% plot(76.5*ones(11),[-0.4:0.64:6],'k--','LineWidth',1.5)
% set(gca,'XLim',[1,182]);
% xt1 = [1:1.7*20:103]; xt2 = [134,160,182];
% set(gca,'XTick',[xt1,xt2]);
% set(gca,'XTickLabel',[1960,1980,2000,2020,2050,2075,2100],'FontSize',ftsz);
% set(gca,'YLim',[-.4,7],'YTick',[0:1:7],'YColor','k');
% text(120,0.3,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')
% text(20,5.3,'historical','fontsize',ftsz-2,'Color','k')
% text(110,5.3,'rcp85','fontsize',ftsz-2,'Color','k')
% 2020
set(gca,'YLim',[-.1,.5],'YTick',[-.1:.1:.5],'YColor','k');
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',ftsz);
text(40,-0.2,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')

legend([p1, p2, p3, p4],'IAP','EN4.2.1','Ishii','CESM1-LE','Location','north','Orientation','horizontal')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240812/TimeSeries2020_OHC_Pacific_PW_',depthstr,'_35S_90S_&CESM_baseline',bsln,'.png'],'-dpng','-r300')
%% plot Atl-IO 1960-2020
close all;
bsln = '1960-1970'; % baseline
% ln1 = (sia_IAP-mean(sia_IAP(51:81)))/10^21;
% ln2 = (sia_EN4-mean(sia_EN4(41:71)))/10^21;
% ln3 = (sia_Ishii-mean(sia_Ishii(36:66)))/10^21;
% ln4 = (sia_LE-mean(sia_LE(71:101)))/10^21;
ln1 = (sia_IAP(21:81)-mean(sia_IAP(21:31)))/siaV_IAP; % 1960-2020
ln2 = (sia_EN4(11:71)-mean(sia_EN4(11:21)))/siaV_IAP;
ln3 = (sia_Ishii(6:66)-mean(sia_Ishii(6:16)))/siaV_IAP;
ln4 = (sia_LE(41:end-1)-mean(sia_LE(41:51)))/siaV_IAP;
ln5 = (sia_LE(41:end-1)-mean(sia_LE(41:51))+stdia(1:end-1)')/siaV_IAP;
ln6 = (sia_LE(41:end-1)-mean(sia_LE(41:51))-stdia(1:end-1)')/siaV_IAP;
%
close all
ftsz = 20;
Fig = figure('position',[700 50 700 400]);
% 2100
% x1 = [1:1.7:103];
% x2 = [104:1:182];
% x = [x1,x2];
% h = fill([x,fliplr(x)],[ln5',fliplr(ln6')],[.9 .9 .9])
% 2020
x1 = [1:61];
x = x1;
h = fill([x,fliplr(x)],[ln5(1:61)',fliplr(ln6(1:61)')],[.9 .9 .9])

hold on
p1 = plot(x1,ln1,'-','color',[0.78,0.58,0.11],'LineWidth',3);
p2 = plot(x1,ln2,'-','color',[.04 .35 .56],'LineWidth',3);
p3 = plot(x1,ln3,'-','color',[.64 .08 .18],'LineWidth',3);
p4 = plot(x,ln4(1:length(x)),'-','color',[.6 .6 .6],'LineWidth',3);
plot(x,zeros(length(x)),'k-','LineWidth',1.5)
ylabel('Temperature  (K)','FontSize',ftsz);
set(gca,'XGrid','on','YGrid','on');
xlabel('Year','FontSize',ftsz);
% 2100
% plot(76.5*ones(11),[-0.4:0.64:6],'k--','LineWidth',1.5)
% set(gca,'XLim',[1,182]);
% xt1 = [1:1.7*20:103]; xt2 = [134,160,182];
% set(gca,'XTick',[xt1,xt2]);
% set(gca,'XTickLabel',[1960,1980,2000,2020,2050,2075,2100],'FontSize',ftsz);
% set(gca,'YLim',[-.4,7],'YTick',[0:1:7],'YColor','k');
% text(120,0.3,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')
% text(20,5.3,'historical','fontsize',ftsz-2,'Color','k')
% text(110,5.3,'rcp85','fontsize',ftsz-2,'Color','k')
% 2020
set(gca,'YLim',[-.1,.5],'YTick',[-.1:.1:.5],'YColor','k');
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',ftsz);
text(40,-0.2,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')

legend([p1, p2, p3, p4],'IAP','EN4.2.1','Ishii','CESM1-LE','Location','north','Orientation','horizontal')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240812/TimeSeries2020_OHC_Atl-IO_PW_',depthstr,'_35S_90S_&CESM_baseline',bsln,'.png'],'-dpng','-r300')
%% plot PAC & Atl-IO
close all;
Fig = figure('position',[700 50 550 800]);
p1 = plot(spac_IAP(21:81)/10^15,':','color','k','LineWidth',2);
hold on
p2 = plot(spac_EN4(11:71)/10^15,':','color',[.04 .35 .56],'LineWidth',2);
p3 = plot(spac_Ishii(6:66)/10^15,':','color',[.64 .08 .18],'LineWidth',2);
p4 = plot(spac_LE(41:101)/10^15,':','color',[0.78,0.58,0.11],'LineWidth',2);

p5 = plot(sia_IAP(21:81)/10^15,'-','color','k','LineWidth',2);
p6 = plot(sia_EN4(11:71)/10^15,'-','color',[.04 .35 .56],'LineWidth',2);
p7 = plot(sia_Ishii(6:66)/10^15,'-','color',[.64 .08 .18],'LineWidth',2);
p8 = plot(sia_LE(41:101)/10^15,'-','color',[0.78,0.58,0.11],'LineWidth',2);
set(gca,'YLim',[110,200],'YTick',[110:10:200],'YColor','k');
ylabel('0-700m OHC  (PW)','FontSize',14);
% plot([1:61],polyval(yb,[1:61]),'--','Color',[.04 .35 .56],'linewidth',2)
% text(50,218,[num2str(roundn(yb(1),-4)*10),' PW per decade'],'fontsize',14,'Color',[.04 .35 .56])
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12,'FontSize',14);
text(28,148,'Pacific','FontSize',12);
l1 = legend([p1, p2, p3, p4],'IAP','EN4.2.1','Ishii','CESM-LE');
set(l1,'Position',[0.155,0.4,0.722,0.031],'Orientation','horizontal');
legend('boxoff')
text(20,198,'Atlantic-Indian Ocean','FontSize',12);
ah = axes('position',get(gca,'position'),'visible','off');
l2 = legend(ah,[p5, p6, p7, p8],'IAP','EN4.2.1','Ishii','CESM1-LE');
set(l2,'Position',[0.155,0.85,0.722,0.031],'Orientation','horizontal','fontsize',12);
legend('boxoff')
% prin t(Fig,['/Volumes/Togo4T/figures/SOtrend_all_Data/20240101/TimeSeries_OHC_Pac_Atl-IO_PW_',depthstr,'m_35S_55S_&CESM.png'],'-dpng','-r300')


function [spac_0 sia_0] = caleach(filename1,filename2,depthData,dx)
% LE ensmean (EX)
addpath /Volumes/Togo4T/1_matlab/help;
TemMMM1 = ncread(filename1,'TEMP');
TemMMM2 = ncread(filename2,'TEMP');
TemMMMall = cat(4,TemMMM1(:,:,36:46,41:end),TemMMM2(:,:,36:46,:));
Temp = TemMMMall;
clear TemMMM1 TemMMM2
% original OHC time series
depthstr = '0-700m';
dweit = depthData(37:46)-depthData(36:45);
lats = [35:55]; % 55S-35S
Tempz = Temp(:,:,2:end,:).*(permute(dweit,[3 2 1]));
Tempzy = Tempz*111*1000; % 111 km / latitude
Tempzyx = Tempzy.*dx*1000; % m
Tzyxsub_r_0 = cat(1,Tempzyx(152:360,:,:,:),Tempzyx(1:151,:,:,:)); 
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,:,:),1),2),3));
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,:,:),1),2),3));
end