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
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthData(2:27)-depthData(1:26); % depth weight

lats = [35:55]; % 55S-35S
levs = [1:27]; % 700 m
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
VV(find(isnan(Temp(:,:,:,1)) == 1)) = nan;
Temp = double(Temp);
Tempz = cat(3,Temp(:,:,1,:),Temp(:,:,2:27,:).*(permute(dweit,[3 2 1])));
VVz = cat(3,VV(:,:,1),VV(:,:,2:27).*(permute(dweit,[3 2 1])));
Tempzy = Tempz*111*1000; % 111 km / latitude
VVzy = VVz*111*1000; % 111km / latitude
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(150:360,:,:,:),Tempzyx(1:149,:,:,:)); 
VVzyx(find(isnan(Tempzyx(:,:,:,1)) == 1)) = nan;
Vzyxsub_r_0 = cat(1,VVzyx(150:360,:,:),VVzyx(1:149,:,:)); 
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),2),3));
spacV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,lats,levs),1),2),3)); % Volume
spac0V = spac_0./spacV;
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),2),3));
siaV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,lats,levs),1),2),3));
sia0V = sia_0./siaV; % /m3
sdiff = spac0V-sia0V;
max(spac0V,[],'all')
sia_IAP = sia0V;
spac_IAP = spac0V;
spacV_IAP = spacV;
siaV_IAP = siaV;
%------------------------ zonal mean -------------------------------------
Tzm_pac = squeeze(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),3));
Szm_pac = squeeze(nansum(nansum(Vzyxsub_r_0(1:151,lats,levs,:),1),3));
zm_pac = Tzm_pac./Szm_pac';
Tzm_aio = squeeze(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),3));
Szm_aio = squeeze(nansum(nansum(Vzyxsub_r_0(152:360,lats,levs,:),1),3));
zm_aio = Tzm_aio./Szm_aio';
x = [1:size(zm_pac(:,1960-1939:2020-1939),2)]';
clear trd 
for i = 1:size(zm_pac,1);
    par1=polyfit(x,zm_pac(i,1960-1939:2020-1939),1); % regression parameters
    trd1_IAP(i) = par1(1); 
    par2=polyfit(x,zm_aio(i,1960-1939:2020-1939),1); % regression parameters
    trd2_IAP(i) = par2(1); 
end
%h0 = trendtest(var,0.01); % t test trend
close all
figure(1)
plot(latData(lats),trd1_IAP,'b')
hold on
plot(latData(lats),trd2_IAP,'r')

yyaxis right
bar(latData(lats),(trd2_IAP-trd1_IAP)./trd1_IAP)
%------------------ 2020 - mean(1960-1970) -------------------------------
zmd_pac_IAP = zm_pac(:,2020-1939)-mean(zm_pac(:,1960-1939:1970-1939),2);
zmd_aio_IAP = zm_aio(:,2020-1939)-mean(zm_aio(:,1960-1939:1970-1939),2);

figure(2)
plot(latData(lats),zmd_pac_IAP,'b')
hold on
plot(latData(lats),zmd_aio_IAP,'r')

yyaxis right
bar(latData(lats),(zmd_aio_IAP-zmd_pac_IAP)./zmd_pac_IAP)
%% -------------- zonal mean MAX contrast (3 data)---------------------------
latsM = [46:47];
spac_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,latsM,levs,:),1),2),3));
spacVmax = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,latsM,levs),1),2),3)); % Volume
spac0Vmax = spac_0max./spacVmax;
sia_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,latsM,levs,:),1),2),3));
siaVmax = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,latsM,levs),1),2),3));
sia0Vmax = sia_0max./siaVmax; % /m3

sia_IAPmax = sia0Vmax;
spac_IAPmax = spac0Vmax;
spacV_IAPmax = spacVmax;
siaV_IAPmax = siaVmax;

figure(3)
plot(sia_IAPmax,'r')
hold on
plot(spac_IAPmax,'b')

%% EN4.2.1 data
clearvars -except sia_IAP spac_IAP siaS_IAP spacS_IAP spacV_IAP siaV_IAP...
    zmd_pac_IAP zmd_aio_IAP trd1_IAP trd2_IAP sia_IAPmax spac_IAPmax ...
    spacV_IAPmax siaV_IAPmax

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

%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthData(2:24)-depthData(1:23); % depth weight

lats = [29:49]; % 55S-35S
levs = [1:25]; % 700 m
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
VV(find(isnan(Temp(:,:,:,1)) == 1)) = nan;
Temp = double(Temp);
Temp700 = Temp(:,:,24,:)+(Temp(:,:,25,:)-Temp(:,:,24,:))*(700-depthData(24))/(depthData(25)-depthData(24));
Tempz = cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:24,:).*(permute(dweit,[3 2 1])),Temp700*(700-depthData(24)));
VV700 = VV(:,:,24,:)+(VV(:,:,25,:)-VV(:,:,24,:))*(700-depthData(24))/(depthData(25)-depthData(24));
VVz = cat(3,VV(:,:,1)*5,VV(:,:,2:24).*(permute(dweit,[3 2 1])),VV700*(700-depthData(24)));
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
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),2),3));
spacV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,lats,levs),1),2),3)); % Volume
spac0V = spac_0./spacV;
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),2),3));
siaV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,lats,levs),1),2),3));
sia0V = sia_0./siaV;
sdiff = spac0V-sia0V;
max(spac0V,[],'all')
sia_EN4 = sia0V;
spac_EN4 = spac0V;

%------------------------ zonal mean ------------------------------------
Tzm_pac = squeeze(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),3));
Szm_pac = squeeze(nansum(nansum(Vzyxsub_r_0(1:151,lats,levs,:),1),3));
zm_pac = Tzm_pac./Szm_pac';
Tzm_aio = squeeze(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),3));
Szm_aio = squeeze(nansum(nansum(Vzyxsub_r_0(152:360,lats,levs,:),1),3));
zm_aio = Tzm_aio./Szm_aio';
x = [1:size(zm_pac(:,1960-1949:2020-1949),2)]';
clear trd 
for i = 1:size(zm_pac,1);
    par1=polyfit(x,zm_pac(i,1960-1949:2020-1949),1); % regression parameters
    trd1_EN4(i) = par1(1); 
    par2=polyfit(x,zm_aio(i,1960-1949:2020-1949),1); % regression parameters
    trd2_EN4(i) = par2(1); 
end
%h0 = trendtest(var,0.01); % t test trend
close all
plot(latData(lats),trd1_EN4,'b')
hold on
plot(latData(lats),trd2_EN4,'r')

yyaxis right
bar(latData(lats),(trd2_EN4-trd1_EN4)./trd1_EN4)
%------------------------ 2020 - mean(1960-1970)--------------------------
zmd_pac_EN4 = zm_pac(:,2020-1949)-mean(zm_pac(:,1960-1949:1970-1949),2);
zmd_aio_EN4 = zm_aio(:,2020-1949)-mean(zm_aio(:,1960-1949:1970-1949),2);
figure(2)
plot(latData(lats),zmd_pac_EN4,'b')
hold on
plot(latData(lats),zmd_aio_EN4,'r')

yyaxis right
bar(latData(lats),(zmd_aio_EN4-zmd_pac_EN4)./zmd_pac_EN4)
%% -------------- zonal mean MAX contrast (3 data)---------------------------
latsM = [40:41];
spac_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,latsM,levs,:),1),2),3));
spacVmax = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,latsM,levs),1),2),3)); % Volume
spac0Vmax = spac_0max./spacVmax;
sia_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,latsM,levs,:),1),2),3));
siaVmax = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,latsM,levs),1),2),3));
sia0Vmax = sia_0max./siaVmax; % /m3

sia_EN4max = sia0Vmax;
spac_EN4max = spac0Vmax;
spacV_EN4max = spacVmax;
siaV_EN4max = siaVmax;

figure(3)
plot(sia_EN4max,'r')
hold on
plot(spac_EN4max,'b')
%% Ishii Data
clearvars -except sia_IAP spac_IAP siaS_IAP spacS_IAP sia_EN4 spac_EN4...
    spacV_IAP siaV_IAP zmd_pac_IAP zmd_aio_IAP zmd_pac_EN4 zmd_aio_EN4...
    trd1_IAP trd2_IAP trd1_EN4 trd2_EN4 sia_IAPmax spac_IAPmax ...
    spacV_IAPmax siaV_IAPmax sia_EN4max spac_EN4max
datadir='/Volumes/Togo4T/data/Ishii/v2017_annual/temp/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'Ishii*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% anomaly 
k=length(filelist); clear Temp
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end

Temps = permute(flipud(Temp),[3 2 1 4]);
depthData = flipud(depthData);
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
dweit = depthData(2:16)-depthData(1:15); % depth weight

lats = [35:55]; % 55S-35S
levs = [1:16]; % 700 m
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
VV(find(isnan(Temp(:,:,:,50)) == 1)) = nan;
Temps = double(Temps);
Tempz = cat(3,Temps(:,:,1,:),Temps(:,:,2:16,:).*(permute(dweit,[3 2 1])));
VVz = cat(3,VV(:,:,1),VV(:,:,2:16).*(permute(dweit,[3 2 1])));
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
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),2),3));
spacV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,lats,levs),1),2),3)); % Volume
spac0V = spac_0./spacV;
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),2),3));
siaV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,lats,levs),1),2),3));
sia0V = sia_0./siaV;
sia_Ishii = sia_0/siaV_IAP; 
spac_Ishii = spac_0/spacV_IAP;
%--------------------- zonal mean ----------------------------------------
Tzm_pac = squeeze(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),3));
Szm_pac = squeeze(nansum(nansum(Vzyxsub_r_0(1:151,lats,levs,:),1),3));
zm_pac = Tzm_pac./Szm_pac';
Tzm_aio = squeeze(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),3));
Szm_aio = squeeze(nansum(nansum(Vzyxsub_r_0(152:360,lats,levs,:),1),3));
zm_aio = Tzm_aio./Szm_aio';
x = [1:size(zm_pac(:,1960-1954:2020-1954),2)]';
clear trd 
for i = 1:size(zm_pac,1);
    par1=polyfit(x,zm_pac(i,1960-1954:2020-1954),1); % regression parameters
    trd1_Ishii(i) = par1(1); 
    par2=polyfit(x,zm_aio(i,1960-1954:2020-1954),1); % regression parameters
    trd2_Ishii(i) = par2(1); 
end
%h0 = trendtest(var,0.01); % t test trend
close all
plot(latData(lats),trd1_Ishii,'b')
hold on
plot(latData(lats),trd2_Ishii,'r')

yyaxis right
bar(latData(lats),(trd2_Ishii-trd1_Ishii)./trd1_Ishii)
%-------------- 2020 - mean(1960-1970) ----------------------------------
zmd_pac_Ishii = zm_pac(:,2020-1954)-mean(zm_pac(:,1960-1954:1970-1954),2);
zmd_aio_Ishii = zm_aio(:,2020-1954)-mean(zm_aio(:,1960-1954:1970-1954),2);
figure(2)
plot(latData(lats),zmd_pac_Ishii,'b')
hold on
plot(latData(lats),zmd_aio_Ishii,'r')

yyaxis right
bar(latData(lats),(zmd_aio_Ishii-zmd_pac_Ishii)./zmd_pac_Ishii)
%% -------------- zonal mean MAX contrast (3 data)---------------------------
latsM = [46:47];
spac_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,latsM,levs,:),1),2),3));
spacVmax = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,latsM,levs),1),2),3)); % Volume
spac0Vmax = spac_0max./spacVmax;
sia_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,latsM,levs,:),1),2),3));
siaVmax = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,latsM,levs),1),2),3));
sia0Vmax = sia_0max./siaVmax; % /m3

sia_Ishiimax = sia0Vmax;
spac_Ishiimax = spac0Vmax;
spacV_Ishiimax = spacVmax;
siaV_Ishiimax = siaVmax;

figure(3)
plot(sia_Ishiimax,'r')
hold on
plot(spac_Ishiimax,'b')
%%
close all
plot(spac_IAP,'b-');
hold on
plot(spac_EN4,'b--');
plot(spac_Ishii,'b.');
plot(sia_IAP,'r-');
plot(sia_EN4,'r--');
plot(sia_Ishii,'r.');
%% 2020 - mean(1960-1970)
close all
y1=(zmd_aio_IAP-zmd_pac_IAP)./zmd_pac_IAP;
y2=(zmd_aio_EN4-zmd_pac_EN4)./zmd_pac_EN4;
y3=(zmd_aio_Ishii-zmd_pac_Ishii)./zmd_pac_Ishii;
plot(latData(lats),y1,'r')
hold on
plot(latData(lats),y2,'b')
plot(latData(lats),y3,'g')
plot(latData(lats),(y1+y2+y3)/3,'k','linewidth',2)
%% trd
figure(2)
y1=(trd2_IAP-trd1_IAP)./trd1_IAP;
y2=(trd2_EN4-trd1_EN4)./trd1_EN4;
y3=(trd2_Ishii-trd1_Ishii)./trd1_Ishii;
yall = [y1;y2;y3];
plot(latData(lats),y1,'r')
hold on
plot(latData(lats),y2,'b')
plot(latData(lats),y3,'g')
plot(latData(lats),(y1+y2+y3)/3,'k','LineWidth',2)
plot(latData(lats),std(yall),'m')
%% CESM1 LE
clearvars -except sia_IAP spac_IAP siaS_IAP spacS_IAP sia_EN4 spac_EN4...
    spacV_IAP siaV_IAP zmd_pac_IAP zmd_aio_IAP zmd_pac_EN4 zmd_aio_EN4...
    trd1_IAP trd2_IAP trd1_EN4 trd2_EN4 sia_IAPmax spac_IAPmax ...
    spacV_IAPmax siaV_IAPmax sia_EN4max spac_EN4max sia_Ishiimax ...
    spac_Ishiimax zmd_pac_Ishii zmd_aio_Ishii trd1_Ishii trd2_Ishii ...
    spac_Ishii sia_Ishii
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

depthstr = '0-700m';
dweit = depthData(2:37)-depthData(1:36);
lats = [35:55]; % 55S-35S
levs = [1:37]; % 700 m
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
VV(find(isnan(Temp(:,:,:,50)) == 1)) = nan;
Tempz = cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:37,:).*(permute(dweit,[3 2 1])));
VVz = cat(3,VV(:,:,1)*5,VV(:,:,2:37).*(permute(dweit,[3 2 1])));
Tempzy = Tempz*111*1000; % 111 km / latitude
VVzy = VVz*111*1000;
clear dx Tempzyx VVzyx
for j = 1:length(latData);
    dx(j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 每个纬度上，经度之间距离不一样
    Tempzyx(:,j,:,:) = Tempzy(:,j,:,:)*dx(j)*1000; % m
    VVzyx(:,j,:) = VVzy(:,j,:)*dx(j)*1000; % m
end
Tzyxsub_r_0 = cat(1,Tempzyx(152:360,:,:,:),Tempzyx(1:151,:,:,:)); 
Vzyxsub_r_0 = cat(1,VVzyx(152:360,:,:),VVzyx(1:151,:,:)); 
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),2),3));
spacV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,lats,levs),1),2),3)); % Volume
spac0V = spac_0./spacV;
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),2),3));
siaV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,lats,levs),1),2),3));
sia0V = sia_0./siaV;
sia_LE = sia0V;
spac_LE = spac0V;
siaV_LE = siaV;
spacV_LE = spacV;
%% -------------- zonal mean MAX contrast (3 data)---------------------------
latsM = [46:47];
spac_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,latsM,levs,:),1),2),3));
spacVmax = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:151,latsM,levs),1),2),3)); % Volume
spac0Vmax = spac_0max./spacVmax;
sia_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,latsM,levs,:),1),2),3));
siaVmax = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(152:360,latsM,levs),1),2),3));
sia0Vmax = sia_0max./siaVmax; % /m3

sia_LEmax = sia0Vmax;
spac_LEmax = spac0Vmax;
spacV_LEmax = spacVmax;
siaV_LEmax = siaVmax;

figure(3)
plot(sia_LEmax,'r')
hold on
plot(spac_LEmax,'b')
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
    [CESMpac(:,s) CESMia(:,s) CESMpacMAX(:,s) CESMiaMAX(:,s)] = caleach(filename1,filename2,depthData,dx);
end
CESMpac = CESMpac/spacV_LE;
CESMia = CESMia/siaV_LE;
CESMpacMAX = CESMpacMAX/spacV_LEmax;
CESMiaMAX = CESMiaMAX/siaV_LEmax;
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

% zonal mean max
stdpacMAX = std(CESMpacMAX',1);
mpacMAX = mean(CESMpacMAX,2)';
stdiaMAX = std(CESMiaMAX',1);
miaMAX = mean(CESMiaMAX,2)';
%% slope 1960-2020
x1=[1:61]; clear yb_pacLE
for s = 1:40
    ln1 = CESMpacMAX(1960-1959:2020-1959,s);
    yb = polyfit(x1,ln1,1);
    yb_pacLE(s) = yb(1);
    ln2 = CESMiaMAX(1960-1959:2020-1959,s);
    yb = polyfit(x1,ln2,1);
    yb_iaLE(s) = yb(1);
end
MEANpac = mean(yb_pacLE)
STDpac = std(yb_pacLE)
MEANia = mean(yb_iaLE)
STDia = std(yb_iaLE)
% contrast (a bit tricky)
ctrst_LE = (MEANia-MEANpac)/MEANpac   % selected
ctst_mean = mean((yb_iaLE-yb_pacLE)./yb_pacLE)
ctst_std = std((yb_iaLE-yb_pacLE)./yb_pacLE)  % selected
%% plot 2100
close all;
bsln = '1960-1970'; % baseline
ftsz = 20;
Fig = figure('position',[700 50 600 600]);
ln1 = (spac_LE(41:end-1)-mean(spac_LE(41:51))+stdpac(1:end-1)');
ln2 = (spac_LE(41:end-1)-mean(spac_LE(41:51))-stdpac(1:end-1)');
ln3 = (spac_LE(41:end-1)-mean(spac_LE(41:51)));
x = [1:140];
h1 = fill([x,fliplr(x)],[ln1',fliplr(ln2')],[0.65,0.65,0.65])
hold on
p1 = plot(x,ln3(1:length(x)),'-','color',[.04 .35 .56],'LineWidth',2);
ln4 = (sia_LE(41:end-1)-mean(sia_LE(41:51)));
ln5 = (sia_LE(41:end-1)-mean(sia_LE(41:51))+stdia(1:end-1)');
ln6 = (sia_LE(41:end-1)-mean(sia_LE(41:51))-stdia(1:end-1)');
p2 = plot(x,ln4(1:length(x)),'-','color',[.64 .08 .18],'LineWidth',2);
h2 = fill([x,fliplr(x)],[ln5',fliplr(ln6')],[0.65,0.65,0.65])
plot(x,zeros(length(x)),'k-','LineWidth',1.5)
set(h1,'facecolor',[0.26,0.69,0.97],'facealpha',0.5,'linewidth',1)
set(h2,'FaceColor',[0.96,0.32,0.32],'facealpha',0.5,'linewidth',1)
xlabel('Year','FontSize',ftsz);
ylabel('Temperature (K)')
% set(gca,'YLim',[-1,8],'YTick',[-1:1:8],'YColor','k');
set(gca,'XLim',[1,140]);
set(gca,'XTick',[1:20:140,140]);
set(gca,'XTickLabel',[1960:20:2100],'FontSize',ftsz);
text(90,-0.25,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')

legend([p1, p2],'Pacific','Atlantic-Indian','Location','north','Orientation','horizontal')
legend('boxoff')
%print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240930/TimeSeriesCESM_2100.png'],'-dpng','-r300')
%% plot Pac 1960-2020 ZJ or GJ/m2
close all;
bsln = '1960-1970'; % baseline
% ln1 = (spac_IAP-mean(spac_IAP(51:81)))/10^21;
% ln2 = (spac_EN4-mean(spac_EN4(41:71)))/10^21;
% ln3 = (spac_Ishii-mean(spac_Ishii(36:66)))/10^21;
% ln4 = (spac_LE-mean(spac_LE(71:101)))/10^21;
ln1 = (spac_IAP(21:81)-mean(spac_IAP(21:31))); % 1960-2020
ln2 = (spac_EN4(11:71)-mean(spac_EN4(11:21)));
ln3 = (spac_Ishii(6:66)-mean(spac_Ishii(6:16)));
ln4 = (spac_LE(41:end-1)-mean(spac_LE(41:51)));
ln5 = (spac_LE(41:end-1)-mean(spac_LE(41:51))+stdpac(1:end-1)');
ln6 = (spac_LE(41:end-1)-mean(spac_LE(41:51))-stdpac(1:end-1)');
% slope
x1=[1:61]
yb = polyfit(x1,(ln1+ln2+ln3)/3,1)
yb1 = polyfit(x1,ln1,1)
yb2 = polyfit(x1,ln2,1)
yb3 = polyfit(x1,ln3,1)
STD=std([roundn(yb1(1),-3);roundn(yb2(1),-3);roundn(yb3(1),-3)])

%%
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
text(40,-0.05,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')

legend([p1, p2, p3, p4],'IAP','EN4.2.1','Ishii','CESM1-LE','Location','north','Orientation','horizontal')
legend('boxoff')
% print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240812/TimeSeries2020_Tem_Pacific_PW_',depthstr,'_35S_90S_&CESM_baseline',bsln,'.png'],'-dpng','-r300')
%% plot Atl-IO 1960-2020
close all;
bsln = '1960-1970'; % baseline
% ln1 = (sia_IAP-mean(sia_IAP(51:81)))/10^21;
% ln2 = (sia_EN4-mean(sia_EN4(41:71)))/10^21;
% ln3 = (sia_Ishii-mean(sia_Ishii(36:66)))/10^21;
% ln4 = (sia_LE-mean(sia_LE(71:101)))/10^21;
ln1 = (sia_IAP(21:81)-mean(sia_IAP(21:31))); % 1960-2020
ln2 = (sia_EN4(11:71)-mean(sia_EN4(11:21)));
ln3 = (sia_Ishii(6:66)-mean(sia_Ishii(6:16)));
ln4 = (sia_LE(41:end-1)-mean(sia_LE(41:51)));
ln5 = (sia_LE(41:end-1)-mean(sia_LE(41:51))+stdia(1:end-1)');
ln6 = (sia_LE(41:end-1)-mean(sia_LE(41:51))-stdia(1:end-1)');
% slope
x1=[1:61]
yb = polyfit(x1,(ln1+ln2+ln3)/3,1)
yb1 = polyfit(x1,ln1,1)
yb2 = polyfit(x1,ln2,1)
yb3 = polyfit(x1,ln3,1)
STD=std([roundn(yb1(1),-3);roundn(yb2(1),-3);roundn(yb3(1),-3)])

%%
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
text(40,-0.05,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')

legend([p1, p2, p3, p4],'IAP','EN4.2.1','Ishii','CESM1-LE','Location','north','Orientation','horizontal')
legend('boxoff')
% print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240812/TimeSeries2020_Tem_Atl-IO_PW_',depthstr,'_35S_90S_&CESM_baseline',bsln,'.png'],'-dpng','-r300')


%% plot Zonal Maximum Pac 1960-2020 ZJ or GJ/m2
close all;
bsln = '1960-1970'; % baseline
ln1 = (spac_IAPmax(21:81)-mean(spac_IAPmax(21:31))); % 1960-2020
ln2 = (spac_EN4max(11:71)-mean(spac_EN4max(11:21)));
ln3 = (spac_Ishiimax(6:66)-mean(spac_Ishiimax(6:16)));
% ln4 = (spac_LEmax(41:end-1)-mean(spac_LEmax(41:51)));
% ln5 = (spac_LEmax(41:end-1)-mean(spac_LEmax(41:51))+stdpacMAX(1:end-1)');
% ln6 = (spac_LEmax(41:end-1)-mean(spac_LEmax(41:51))-stdpacMAX(1:end-1)');
%
close all
ftsz = 20;
Fig = figure('position',[700 50 700 400]);

% 2020
x1 = [1:61];
x = x1;
%h = fill([x,fliplr(x)],[ln5(1:61)',fliplr(ln6(1:61)')],[.9 .9 .9])

hold on
p1 = plot(x1,ln1,'-','color',[0.78,0.58,0.11],'LineWidth',3);
p2 = plot(x1,ln2,'-','color',[.04 .35 .56],'LineWidth',3);
p3 = plot(x1,ln3,'-','color',[.64 .08 .18],'LineWidth',3);
%p4 = plot(x,ln4(1:length(x)),'-','color',[.6 .6 .6],'LineWidth',3);
yb = polyfit(x1,(ln1+ln2+ln3)/3,1);
plot(x1,polyval(yb,x1),'k-','linewidth',4)
yb1 = polyfit(x1,ln1,1)
yb2 = polyfit(x1,ln2,1)
yb3 = polyfit(x1,ln3,1)
STD=std([roundn(yb1(1),-3);roundn(yb2(1),-3);roundn(yb3(1),-3)])
text(35,-0.08,['Slope=',num2str(roundn(yb(1),-3)*10),'\pm',num2str(roundn(STD,-3)*10),' K/decade'],'fontsize',ftsz-2,'Color','k')
plot(x,zeros(length(x)),'k-','LineWidth',1.5)
ylabel('Temperature  (K)','FontSize',ftsz);
set(gca,'XGrid','on','YGrid','on');
xlabel('Year','FontSize',ftsz);

% 2020
set(gca,'YLim',[-.2,.8],'YTick',[-.2:.2:.8],'YColor','k');
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',ftsz);
text(35,-0.15,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')

legend([p1, p2, p3],'IAP','EN4.2.1','Ishii','Location','north','Orientation','horizontal')
legend('boxoff')
%print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20241122/TimeSeries2020_Tem_Pacific_PW_',depthstr,'_43.5S_44.5S_baseline',bsln,'.png'],'-dpng','-r300')

%% plot Zonal Maximum Atl-IO 1960-2020
close all;
bsln = '1960-1970'; % baseline
ln1 = (sia_IAPmax(21:81)-mean(sia_IAPmax(21:31))); % 1960-2020
ln2 = (sia_EN4max(11:71)-mean(sia_EN4max(11:21)));
ln3 = (sia_Ishiimax(6:66)-mean(sia_Ishiimax(6:16)));
% ln4 = (sia_LEmax(41:end-1)-mean(sia_LEmax(41:51)));
% ln5 = (sia_LEmax(41:end-1)-mean(sia_LEmax(41:51))+stdiaMAX(1:end-1)');
% ln6 = (sia_LEmax(41:end-1)-mean(sia_LEmax(41:51))-stdiaMAX(1:end-1)');
%
close all
ftsz = 20;
Fig = figure('position',[700 50 700 400]);

% 2020
x1 = [1:61];
x = x1;
%h = fill([x,fliplr(x)],[ln5(1:61)',fliplr(ln6(1:61)')],[.9 .9 .9])

hold on
p1 = plot(x1,ln1,'-','color',[0.78,0.58,0.11],'LineWidth',3);
p2 = plot(x1,ln2,'-','color',[.04 .35 .56],'LineWidth',3);
p3 = plot(x1,ln3,'-','color',[.64 .08 .18],'LineWidth',3);
%p4 = plot(x,ln4(1:length(x)),'-','color',[.6 .6 .6],'LineWidth',3);
yb = polyfit(x1,(ln1+ln2+ln3)/3,1);
plot(x1,polyval(yb,x1),'k-','linewidth',4)
yb1 = polyfit(x1,ln1,1)
yb2 = polyfit(x1,ln2,1)
yb3 = polyfit(x1,ln3,1)
STD=std([roundn(yb1(1),-3);roundn(yb2(1),-3);roundn(yb3(1),-3)])
text(35,-0.08,['Slope=',num2str(roundn(yb(1),-3)*10),'\pm',num2str(roundn(STD,-3)*10),' K/decade'],'fontsize',ftsz-2,'Color','k')
plot(x,zeros(length(x)),'k-','LineWidth',1.5)
ylabel('Temperature  (K)','FontSize',ftsz);
set(gca,'XGrid','on','YGrid','on');
xlabel('Year','FontSize',ftsz);

% 2020
set(gca,'YLim',[-.2,.8],'YTick',[-.2:.2:.8],'YColor','k');
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',ftsz);
text(35,-0.15,['Baseline ',bsln],'fontsize',ftsz-2,'Color','k')

legend([p1, p2, p3],'IAP','EN4.2.1','Ishii','Location','north','Orientation','horizontal')
legend('boxoff')
%print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20241122/TimeSeries2020_Tem_Atl-IO_PW_',depthstr,'_43.5S_44.5S_baseline',bsln,'.png'],'-dpng','-r300')


function [spac_0 sia_0 spac_0max sia_0max] = caleach(filename1,filename2,depthData,dx)
% LE ensmean (EX)
addpath /Volumes/Togo4T/1_matlab/help;
TemMMM1 = ncread(filename1,'TEMP');
TemMMM2 = ncread(filename2,'TEMP');
TemMMMall = cat(4,TemMMM1(:,:,1:37,41:end),TemMMM2(:,:,1:37,:));
Temp = TemMMMall;
clear TemMMM1 TemMMM2
% original OHC time series
depthstr = '0-700m';
dweit = depthData(2:37)-depthData(1:36);
lats = [35:55]; % 55S-35S
levs = [1:37]; % 700 m
Tempz = cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:37,:).*(permute(dweit,[3 2 1])));
Tempzy = Tempz*111*1000; % 111 km / latitude
Tempzyx = Tempzy.*dx*1000; % m
Tzyxsub_r_0 = cat(1,Tempzyx(152:360,:,:,:),Tempzyx(1:151,:,:,:)); 
spac_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,lats,levs,:),1),2),3));
sia_0 = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,lats,levs,:),1),2),3));
% zonal mean max
latsM = [46:47];
spac_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:151,latsM,levs,:),1),2),3));
sia_0max = squeeze(nansum(nansum(nansum(Tzyxsub_r_0(152:360,latsM,levs,:),1),2),3));
end