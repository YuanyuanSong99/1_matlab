% LE ensmean (EX)
addpath G:\1_matlab\help;
filename1 = 'H:\CESM-post\LE\Temperature\his_ensmean.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
TemMMM1 = ncread(filename1,'TEMP');
TemMMM1(:,:,:,87) = [];
filename2 = 'H:\CESM-post\LE\Temperature\rcp85_sub_ensmean.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
TemMMM2 = ncread(filename2,'TEMP');
TemMMMall = cat(4,TemMMM1,TemMMM2(:,1:90,:,:));
TemMMMa = TemMMMall-mean(TemMMMall,4);
lonData = ncread(filename2,'lon');
latData = ncread(filename2,'lat');
depthData = ncread(filename1,'z_t')/100; % unit: cm -> m
dweit = depthData(2:end)-depthData(1:end-1);
TemMMMsub = permute(nansum(TemMMMa(:,:,1,:)*5+TemMMMa(:,:,2:37,:).*permute(dweit,[3 2 1 4]),3)/707,[1 2 4 3]);
lats = [29:44]; % latitude
Tsub_r = cat(1,TemMMMsub(160:300,:,:),TemMMMsub(301:360,:,:),TemMMMsub(1:159,:,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
[spacz_long spac_long] = areamean(Tsub_r,1:141,lats,latData); 
[siaz_long sia_long] = areamean(Tsub_r,142:360,lats,latData); 
DI_rawLE = spac_long-sia_long;
dT = 1;  cf = 1/8;
DIf_LE = lanczosfilter(DI_rawLE,dT,cf,[],'low'); % 8 year filtered IPO index

% TPI-LE
filename1 = 'H:\CESM-post\LE\Temperature\his_sst_ensmean.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
sstMMM1 = ncread(filename1,'TEMP');
sstMMM1(:,:,:,87) = [];
filename2 = 'H:\CESM-post\LE\Temperature\rcp85_sst_ensmean.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
sstMMM2 = ncread(filename2,'TEMP');
sstMMMall = cat(4,sstMMM1,sstMMM2);
sstMMMa = permute(sstMMMall-mean(sstMMMall,4),[1 2 4 3]);
% TPI index ---------------------------------------------------------------
[ts1_zs ts1] = areamean(sstMMMa,141:216,116:136,latData); % 25N-45N,140E-145W
[ts2_zs ts2] = areamean(sstMMMa,171:271,80:101,latData); % 10S-10N,170E-90W
[ts3_zs ts3] = areamean(sstMMMa,151:201,40:75,latData); % 50S-15S,150E-160W
TPI1 = ts2-(ts1+ts3)/2; % unfiltered IPO index///
dT = 1;  cf = 1/8;
TPIf_LE = lanczosfilter(TPI1,dT,cf,[],'low'); % 8 year filtered IPO index

%% internal variability
clc,clear,close all;
addpath G:\1_matlab\help;
filename = 'H:\CESM-post\LE\Temperature\his_ensmean.nc';
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
TemMMM = ncread(filename,'TEMP');
lonData = ncread(filename,'lon');
latData = ncread(filename,'lat');
depthData = ncread(filename,'z_t')/100; % unit: cm -> m
save('MatFile/lonData.mat','lonData');
save('MatFile/latData.mat','latData');
save('MatFile/depthData.mat','depthData');% clear filename
TemMMMyr_cli = mean(TemMMM,4);
datadir='H:\CESM-post\LE\Temperature\his_sub_yr1x1\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
for s=1:k
    s 
    clearvars -except TemMMM TemMMMyr_cli datadir filelist k s lonData latData
    filename=[datadir,filelist(s).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    Tem = ncread(filename,'TEMP');
    Temd = Tem - TemMMM; % detrend
    Temda = Temd - mean(Temd,4); % anomaly
    % 8-yr filter
    dT = 1;  cf = 1/8;
    Temsubar = reshape(permute(Temda,[4 1 2 3]),[size(Temda,4),length(lonData)*length(latData)*size(Temda,3)]);
    parfor i = 1:size(Temsubar,2);
        Temsubadf(:,i) = lanczosfilter(Temsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    Temsubadf = permute(reshape(Temsubadf,[size(Temda,4),length(lonData),length(latData),size(Temda,3)]),[2 3 4 1]);  % -5~707m, 8-year filtered

    save(['MatFile/0_700m_South/Temsubadf',num2str(s),'.mat'],'Temsubadf');
end
%% sst
clc,clear,close all;
addpath G:\1_matlab\help;
filename = 'E:\CESM-post\JC-data\LE\Temperature\his_sst_ensmean.nc';
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
sstMMM = ncread(filename,'TEMP');
lonData = ncread(filename,'lon');
latData = ncread(filename,'lat');
% clear filename
sstMMMyr_cli = mean(sstMMM,4);
datadir='E:\CESM-post\JC-data\LE\Temperature\his_sst_yr1x1\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
for s=1:k
    s 
    clearvars -except sstMMM datadir filelist k s lonData latData
    filename=[datadir,filelist(s).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    Tem = ncread(filename,'TEMP');
    Temd = Tem - sstMMM; % minus ensmean (detrend)
    Temda = Temd - mean(Temd,4); % anomaly  
    % 8-yr filter
    dT = 1;  cf = 1/8;
    Temsubar = reshape(permute(Temda,[4 1 2 3]),[size(Temda,4),length(lonData)*180]);
    parfor i = 1:size(Temsubar,2);
        Temsubadf(:,i) = lanczosfilter(Temsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    sstadf = permute(reshape(Temsubadf,[size(Temda,4),length(lonData),180]),[2 3 1]);  % sst, 8-year filtered
    save(['MatFile/sst_Global/sstadf',num2str(s),'.mat'],'sstadf');
end
%% taux
clc,clear,close all;
addpath G:\1_matlab\help;
filename = 'H:\CESM-post\LE\TAU\historical\his_taux_ensmean.nc';
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
tauxMMM = ncread(filename,'TAUX');
lonData = ncread(filename,'lon');
latData = ncread(filename,'lat');
% clear filename
tauxMMMyr_cli = mean(tauxMMM,4);
datadir='H:\CESM-post\LE\TAU\historical\taux_yr_1x1\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
for s=1:k
    s 
    filename=[datadir,filelist(s).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    Taux = ncread(filename,'TAUX');
    Tauxd = Taux - tauxMMM; % minus ensmean (detrend)
    Tauxda = Tauxd - mean(Tauxd,3); % anomaly  
    [d1 d2 d3] = size(Tauxda);
    % 8-yr filter
    dT = 1;  cf = 1/8;
    Tauxsubar = reshape(permute(Tauxda,[3 1 2]),[d3,d1*d2]);
    parfor i = 1:size(Tauxsubar,2);
        Tauxsubadf(:,i) = lanczosfilter(Tauxsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    tauxadf = permute(reshape(Tauxsubadf,[d3,d1,d2]),[2 3 1]);  % sst, 8-year filtered
    save(['MatFile/taux_Global/tauxadf',num2str(s),'.mat'],'tauxadf');
end
%% tauy
clc,clear,close all;
addpath G:\1_matlab\help;
filename = 'H:\CESM-post\LE\TAU\historical\his_tauy_ensmean.nc';
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
tauyMMM = ncread(filename,'TAUY');
lonData = ncread(filename,'lon');
latData = ncread(filename,'lat');
% clear filename
tauyMMMyr_cli = mean(tauyMMM,4);
datadir='H:\CESM-post\LE\TAU\historical\tauy_yr_1x1\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
for s=1:k
    s 
    filename=[datadir,filelist(s).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    Tauy = ncread(filename,'TAUY');
    Tauyd = Tauy - tauyMMM; % minus ensmean (detrend)
    Tauyda = Tauyd - mean(Tauyd,3); % anomaly  
    [d1 d2 d3] = size(Tauyda);
    % 8-yr filter
    dT = 1;  cf = 1/8;
    Tauysubar = reshape(permute(Tauyda,[3 1 2]),[d3,d1*d2]);
    parfor i = 1:size(Tauysubar,2);
        Tauysubadf(:,i) = lanczosfilter(Tauysubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    tauyadf = permute(reshape(Tauysubadf,[d3,d1,d2]),[2 3 1]);  % sst, 8-year filtered
    save(['MatFile/tauy_Global/tauyadf',num2str(s),'.mat'],'tauyadf');
end
%% curlz
filename1 = 'H:\CESM-post\LE\TAU\historical\his_taux_ensmean.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
tauxMMM = ncread(filename1,'TAUX');
filename2 = 'H:\CESM-post\LE\TAU\historical\his_tauy_ensmean.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
tauyMMM = ncread(filename2,'TAUY');
for i = 1:size(tauxMMM,3);
    curlzMMM(:,:,i) = ra_windstrcurl(latData,lonData,tauxMMM(:,:,i)',tauyMMM(:,:,i)',1);  % curlz lat*lon
end
datadir1='H:\CESM-post\LE\TAU\historical\taux_yr_1x1\'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='H:\CESM-post\LE\TAU\historical\tauy_yr_1x1\'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist2);
for s=1:k
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Taux = ncread(filename1,'TAUX');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Tauy = ncread(filename2,'TAUY');
    for i = 1:size(Taux,3);
        Curlz(:,:,i) = ra_windstrcurl(latData,lonData,Taux(:,:,i)',Tauy(:,:,i)',1);  % curlz lat*lon
    end
    Curlzd = Curlz - curlzMMM; % minus ensmean (detrend)
    Curlzda = Curlzd - mean(Curlzd,3); % anomaly  
    [d1 d2 d3] = size(Curlzda);
    % 8-yr filter
    dT = 1;  cf = 1/8;
    Curlzsubar = reshape(permute(Curlzda,[3 1 2]),[d3,d1*d2]);
    parfor i = 1:size(Curlzsubar,2);
        Curlzsubadf(:,i) = lanczosfilter(Curlzsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    curlzadf = permute(reshape(Curlzsubadf,[d3,d1,d2]),[3 2 1]);  % sst, 8-year filtered
    save(['MatFile/curlz_Global/curlzadf',num2str(s),'.mat'],'curlzadf');
end

%% slp
clc,clear,close all;
addpath G:\1_matlab\help;
filename = 'H:\CESM-post\LE\psl\historical\his_psl_ensmean.nc';
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
slpMMM = ncread(filename,'PSL');
% lonData = ncread(filename,'lon');
% latData = ncread(filename,'lat');
% clear filename
slpMMMyr_cli = mean(slpMMM,4);
datadir='H:\CESM-post\LE\psl\historical\psl_yr_1x1\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
for s=1:k
    s 
    filename=[datadir,filelist(s).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    Psl = ncread(filename,'PSL');
    Psld = Psl - slpMMM; % minus ensmean (detrend)
    Pslda = Psld - mean(Psld,3); % anomaly  
    [d1 d2 d3] = size(Pslda);
    % 8-yr filter
    dT = 1;  cf = 1/8;
    Pslsubar = reshape(permute(Pslda,[3 1 2]),[d3,d1*d2]);
    parfor i = 1:size(Pslsubar,2);
        Pslsubadf(:,i) = lanczosfilter(Pslsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    slpadf = permute(reshape(Pslsubadf,[d3,d1,d2]),[2 3 1]);  % sst, 8-year filtered
    save(['MatFile/slp_Global/slpadf',num2str(s),'.mat'],'slpadf');
end
%% SHF
% clc,clear,close all;
addpath G:\1_matlab\help;
filename = 'H:\CESM-post\LE\SHF\historical\his_SHF_ensmean.nc';
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
shfMMM = ncread(filename,'SHF');
% lonData = ncread(filename,'lon');
% latData = ncread(filename,'lat');
% clear filename
shfMMMyr_cli = mean(shfMMM,4);
datadir='H:\CESM-post\LE\SHF\historical\SHF_yr_1x1\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
%
for s=1
    s 
    filename=[datadir,filelist(s).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    SHF = ncread(filename,'SHF');
    SHFd = SHF - shfMMM; % minus ensmean (detrend)
    SHFda = SHFd - mean(SHFd,3); % anomaly  
    [d1 d2 d3] = size(SHFda);
    % 8-yr filter
    dT = 1;  cf = 1/8;
    SHFsubar = reshape(permute(SHFda,[3 1 2]),[d3,d1*d2]);
    parfor i = 1:size(SHFsubar,2);
        SHFsubadf(:,i) = lanczosfilter(SHFsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    shfadf = permute(reshape(SHFsubadf,[d3,d1,d2]),[2 3 1]);  % sst, 8-year filtered
%     save(['MatFile/shf_Global/shfadf',num2str(s),'.mat'],'shfadf');
end
%% NHFLX
% FLNS
clc,clear,close all;
addpath G:\1_matlab\help;
filename_le = 'H:\CESM-post\LE\FLNS\historical\his_FLNS_ensmean.nc';
ncid=netcdf.open(filename_le,'NOWRITE');
ncdisp(filename_le);
flnsMMM = ncread(filename_le,'FLNS');
flnsMMMyr_cli = mean(flnsMMM,4);
datadir_l='H:\CESM-post\LE\FLNS\historical\FLNS_yr_1x1\'; %指定批量数据所在的文件夹
filelist_l=dir([datadir_l,'*.nc']); %指定批量数据的类型
ncdisp([datadir_l,filelist_l(1).name]);
k=length(filelist_l);
% FSNS
addpath G:\1_matlab\help;
filename_se = 'H:\CESM-post\LE\FSNS\historical\his_FSNS_ensmean.nc';
ncid=netcdf.open(filename_se,'NOWRITE');
ncdisp(filename_se);
fsnsMMM = ncread(filename_se,'FSNS');
fsnsMMMyr_cli = mean(fsnsMMM,4);
datadir_s='H:\CESM-post\LE\FSNS\historical\FSNS_yr_1x1\'; %指定批量数据所在的文件夹
filelist_s=dir([datadir_s,'*.nc']); %指定批量数据的类型
ncdisp([datadir_s,filelist_s(1).name]);
% LHFLX
addpath G:\1_matlab\help;
filename_lhe = 'H:\CESM-post\LE\LHFLX\historical\his_LHFLX_ensmean.nc';
ncid=netcdf.open(filename_lhe,'NOWRITE');
ncdisp(filename_lhe);
lhflxMMM = ncread(filename_lhe,'LHFLX');
lhflxMMMyr_cli = mean(lhflxMMM,4);
datadir_lh='H:\CESM-post\LE\LHFLX\historical\LHFLX_yr_1x1\'; %指定批量数据所在的文件夹
filelist_lh=dir([datadir_lh,'*.nc']); %指定批量数据的类型
ncdisp([datadir_lh,filelist_lh(1).name]);
% SHFLX
addpath G:\1_matlab\help;
filename_she = 'H:\CESM-post\LE\SHFLX\historical\his_SHFLX_ensmean.nc';
ncid=netcdf.open(filename_she,'NOWRITE');
ncdisp(filename_she);
shflxMMM = ncread(filename_she,'SHFLX');
shflxMMMyr_cli = mean(shflxMMM,4);
datadir_sh='H:\CESM-post\LE\SHFLX\historical\SHFLX_yr_1x1\'; %指定批量数据所在的文件夹
filelist_sh=dir([datadir_sh,'*.nc']); %指定批量数据的类型
ncdisp([datadir_sh,filelist_sh(1).name]);
%%
nhflxMMM = fsnsMMM-shflxMMM-lhflxMMM-flnsMMM;
map = mean(mean(nhflxMMM,3),4);
max(map,[],'all')
min(map,[],'all')
%
close all;
lonData = ncread(filename_she,'lon');
latData = ncread(filename_she,'lat');
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
    hold on
    m_contourf(lonData,latData,map',[-300:20:300],'linestyle','none');
    caxis([-100,100]);
    load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    colorbar;
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
%%
for s=1:k
    s 
    filename_l=[datadir_l,filelist_l(s).name];
    ncid=netcdf.open(filename_l,'NC_NOWRITE');
    FLNS = ncread(filename_l,'FLNS');
    filename_s=[datadir_s,filelist_s(s).name];
    ncid=netcdf.open(filename_s,'NC_NOWRITE');
    FSNS = ncread(filename_s,'FSNS');
    filename_lh=[datadir_lh,filelist_lh(s).name];
    ncid=netcdf.open(filename_lh,'NC_NOWRITE');
    LHFLX = ncread(filename_lh,'LHFLX');
    filename_sh=[datadir_sh,filelist_sh(s).name];
    ncid=netcdf.open(filename_sh,'NC_NOWRITE');
    SHFLX = ncread(filename_sh,'SHFLX');
    NHFLX(:,:,:,s) = FSNS-SHFLX-LHFLX-FLNS;
end
%  
for s = 1:k
    s
    NHFLXd = NHFLX(:,:,:,k) - nhflxMMM; % minus ensmean (detrend)
    NHFLXda = NHFLXd - mean(NHFLXd,3); % anomaly  
    [d1 d2 d3] = size(NHFLXda);
    % 8-yr filter
    dT = 1;  cf = 1/8;
    NHFLXsubar = reshape(permute(NHFLXda,[3 1 2]),[d3,d1*d2]);
    parfor i = 1:size(NHFLXsubar,2);
        NHFLXsubadf(:,i) = lanczosfilter(NHFLXsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    nhflxadf = permute(reshape(NHFLXsubadf,[d3,d1,d2]),[2 3 1]);  % sst, 8-year filtered
    save(['MatFile/nhflx_Global/nhflxadf',num2str(s),'.mat'],'nhflxadf');
end
% clear FLNS FLNSd FLNSda FLNSsubar FLNSsubadf flnsadf





function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end


