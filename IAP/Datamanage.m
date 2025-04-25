clc,clear,close all;
addpath G:\1_matlab\help;
datadir='G:\data\IAP\temperature\monthly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'CZ16*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE');
lonData = ncread(filetemp,'lon');
latData = ncread(filetemp,'lat');
depthData = ncread(filetemp,'depth_std');
% save('MatFile/lonData.mat','lonData');
% save('MatFile/latData.mat','latData');
% save('MatFile/depthData.mat','depthData');
%%  
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp(:,:,:,s) = ncread(filename1,'temp'); %读入变量
    netcdf.close(ncid);  %关闭文件    
end
%%
Temp1 = reshape(Temp,)