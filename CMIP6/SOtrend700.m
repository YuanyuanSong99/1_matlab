clc,clear,close all;
addpath(genpath('D:\1_matlab\help'));
filename1 = 'F:\CMIP6\historical\ptem\orggrid\thetao_year_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412_horgrid.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
Tem = ncread(filename1,'thetao');
%%
lat = ncread(filename1,'lat');
lat(1:60)
lon = ncread(filename1,'lon');
lev = ncread(filename1,'lev');
%
levdp = find(lev-700 > 0);
lev800 = lev(levdp(1))  % > 700m
lev600 = lev(levdp(1)-1)  % < 700m
Tem700 = Tem(:,1:60,1:levdp(1),111:end); % 90S-30S; 1-800m; 1960-2014;
clear Tem
Tem700a = Tem700 - nanmean(Tem700(:,:,:,22:51),4); % remove climatology from 1981-2010
depthstr = '0-700';
dweit = lev(2:lev800)-lev(1:lev600); % depth weight
Tsubraw = permute(nansum(cat(3,Temp(:,lats,1,:),Temp(:,lats,2:27,:).*(permute(dweit,[3 2 1]))),3)/700,[1 2 4 3]); 
