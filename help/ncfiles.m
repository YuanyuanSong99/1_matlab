%% see 1 nc file information
% clc,clear,close all;
filename = 'E:\data\IAPocean\Temperature\CZ16_1_2000m_Temp_year_1940_month_01.nc';
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
level = ncread(filename,'depth_std');
% lonData = ncread(filename,'lon');
% tp = ncread('G:\data\hadley\1950-2020\HadISST_ice.nc','time');
% sst = ncread(filename,'sst');
%% CESM vertical levels
levelData = ncread(filename,'lev');
p0 = ncread(filename,'P0');
hyam = ncread(filename,'hyam');
hybm = ncread(filename,'hybm');
ps = ncread(filename,'PS');
psp = mean(ps(97:129,128:166,335:end),3); % EAT area
for h = 1:length(hyam);
    p(:,:,h) = hyam(h)*p0+hybm(h)*psp;
end
mean(mean(p,1),2)
% z500 = ncread('E:\sea ice\data\ERA-Interim\sst_daily_4time_2.5\sst.daily.1979.nc','z');
%% contourf
close all;
% lons 
lats = latData(1:180);
var = sst(:,1:180);
val = [min(var,[],'all'):(max(var,[],'all')-min(var,[],'all'))/20:max(var,[],'all')];
m_proj('stereographic','lon',150,'lat',90,'radius',90,'rotangle',0);
hold on
m_contourf(lonData,lats,var',val,'linestyle','none');
% caxis([-ctr,ctr]);
load('G:/1_matlab/help/colorbar_mat/bl_re2.mat');
colormap(double(bl_re2)/255);
m_coast('linewidth',0.8,'color','k');
m_grid('xtick',13,'ytick',0,'yticklabel',[],'box','on','linestyle','--');

%% read lots of nc files
clc,clear,close all;
datadir='D:\sea ice\data\NCEP\Geopotential_height\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
a=filelist(1).name; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
b=filelist(2).name; %查看你要读取的文件的编号。filelist（2）.name在window下为第二个标号的数据
k=length(filelist);
z500_djf=cell(1,k-1);
for s=1:k-1
    filename1=[datadir,filelist(s).name];
    filename2=[datadir,filelist(s+1).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    u1 = ncread(filename1,'hgt'); %读入变量
    u2 = ncread(filename2,'hgt');
    timeData1 = ncread(filename1,'time');
    timeData2 = ncread(filename2,'time');
    lonData = ncread(filename1,'lon');
    latData = ncread(filename1,'lat');
    netcdf.close(ncid);  %关闭文件
    t_num1=length(timeData1);
    t_num2=length(timeData2);
    if t_num1 == 366
        u_dec = u1(:,1:37,6,336:366);
    else
        u_dec = u1(:,1:37,6,335:365);
    end
    if t_num2 == 366
        u_jf = u2(:,1:37,6,1:60);
    else
        u_jf = u2(:,1:37,6,1:59);
    end
    z500_djf{1,s}=cat(4,u_dec,u_jf); %40个DJF
    z500_djf{1,s}=permute(z500_djf{1,s},[1 2 4 3]);
%     l=size(u_djf{1,s},3);
%     for i=1:l/4
%         z_djf_daily{1,s}(:,:,i)=mean(u_djf{1,s}(:,:,4*(i-1)+1:4*i),3); %40个DJF，daily
%     end
end
save('D:\1_matlab\NCEP\z500_djf_daily.mat','z500_djf');
% save('D:\1_matlab\NCEP\lon.mat','lonData');
% latData=latData(1:37);
% save('D:\1_matlab\NCEP\lat.mat','latData');