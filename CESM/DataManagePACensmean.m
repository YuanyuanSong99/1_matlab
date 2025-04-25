%% ensmean data
% TEMP
clc,clear,close all;
addpath D:\1_matlab\help;
filename1 = 'E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\TEMP_192001-200512_ensmean.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
TemMMM1 = ncread(filename1,'TEMP');
TemMMM1(:,:,:,end) = []; % remove the last year (only January)
filename2 = 'E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\TEMP_200601-201312_ensmean.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
TemMMM2 = ncread(filename2,'TEMP');
TemMMM2(:,:,:,end) = []; % remove the last year (only January)
lonData = ncread(filename1,'lon');
latData = ncread(filename1,'lat');
depthData = ncread(filename1,'z_t')/100; % unit: cm -> m
TemMMM = cat(4,TemMMM1,TemMMM2);
TemMMMsub = TemMMM(:,1:90,1:37,:); % southern hemisphere; up 707m 
% minus LE ensemble mean
filename3 = 'E:\CESM-post\JC-data\LE\Temperature\his_ensmean.nc';
ncid=netcdf.open(filename3,'NOWRITE');
ncdisp(filename3);
TemLEMsub1 = ncread(filename3,'TEMP');
TemLEMsub1(:,:,:,end) = [];
filename4 = 'E:\CESM-post\JC-data\LE\Temperature\rcp85_sub_ensmean.nc';
ncid=netcdf.open(filename4,'NOWRITE');
ncdisp(filename4);
TemLEMsub2 = ncread(filename4,'TEMP');
TemLEMsub = cat(4,TemLEMsub1,TemLEMsub2(:,1:90,:,:));
TemMMMsubd = TemMMMsub-TemLEMsub; % detrend
TemMMMsubda = TemMMMsubd-mean(TemMMMsubd,4); % anomaly
% detrend and 8-yr filter
dT = 1;  cf = 1/8;
Temda = TemMMMsubda;
[d1 d2 d3 d4] = size(Temda);
Temsubar = reshape(permute(Temda,[4 1 2 3]),[d4,d1*d2*d3]);
parfor i = 1:size(Temsubar,2);
    Temsubadf_long(:,i) = lanczosfilter(Temsubar(:,i),dT,cf,[],'low'); % 8 year filtered
end
Temsubadf_long = permute(reshape(Temsubadf_long,[d4,d1,d2,d3]),[2 3 4 1]);  % -5~707m, 8-year filtered
Temsubadf = Temsubadf_long(:,:,:,5:end-4); 
clear TemLEMsub1 TemLEMsub2 TemMMM TemMMM1 TemMMM2 Temsubar
% SST  lon*lat*1*time
addpath E:\1_matlab\help;
filename1 = 'E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\SST_192001-200512_ensmean.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
sstMMM1 = ncread(filename1,'SST');
sstMMM1(:,:,:,end) = []; % remove the last year (only January)
filename2 = 'E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\SST_200601-201312_ensmean.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
sstMMM2 = ncread(filename2,'SST');
sstMMM2(:,:,:,end) = []; % remove the last year (only January)
sstMMM = cat(4,sstMMM1,sstMMM2);
% minus LE ensemble mean
filename3 = 'E:\CESM-post\JC-data\LE\Temperature\his_sst_ensmean.nc';
ncid=netcdf.open(filename3,'NOWRITE');
ncdisp(filename3);
sstLEMsub1 = ncread(filename3,'TEMP');
sstLEMsub1(:,:,:,end) = []; % remove the last year (only January)
filename4 = 'E:\CESM-post\JC-data\LE\Temperature\rcp85_sst_ensmean.nc';
ncid=netcdf.open(filename4,'NOWRITE');
ncdisp(filename4);
sstLEMsub2 = ncread(filename4,'TEMP');
sstLEMsub = cat(4,sstLEMsub1,sstLEMsub2);

sstMMMsubd = sstMMM-sstLEMsub; % detrend
sstMMMsubda = sstMMMsubd - mean(sstMMMsubd,4); % anomaly
% detrend and 8-yr filter
dT = 1;  cf = 1/8;
sstad = sstMMMsubda;
[d1 d2 d3 d4] = size(sstad);
sstar = reshape(permute(sstad,[4 1 2 3]),[d4,d1*d2*d3]);
clear sstadf
parfor i = 1:size(sstar,2);
    sstadf_long(:,i) = lanczosfilter(sstar(:,i),dT,cf,[],'low'); % 8 year filtered
end
sstadf = permute(reshape(sstadf_long,[d4,d1,d2,d3]),[2 3 1 4]);  % -5~707m, 8-year filtered
sstadf(:,:,1:4) = []; sstadf(:,:,end-3:end) = []; 
clear sstLEM* sstMM* sstar
%% single level LE ensmean
[tauxadf] = DATAadf('TAUX','E:\CESM-post\JC-data\LE\TAUX\historical\his_taux_ensmean.nc');
[tauyadf] = DATAadf('TAUY','E:\CESM-post\JC-data\LE\TAUY\historical\his_tauy_ensmean.nc');
tauxadf(:,:,1:4) = []; tauxadf(:,:,end-3:end) = [];
taux = double(tauxadf)/10; % unit: dyn/cm2 -> N/m2
tauyadf(:,:,1:4) = []; tauyadf(:,:,end-3:end) = [];
tauy = double(tauyadf)/10; % unit: dyn/cm2 -> N/m2
% wind stress curl
parfor i = 1:size(taux,3);
    % curlz lat*lon
    curlz(:,:,i) = ra_windstrcurl(latData,lonData,taux(:,:,i)',tauy(:,:,i)',1);
end
%
[slpadf] = DATAadf('PSL','E:\CESM-post\JC-data\LE\PSL\historical\his_psl_ensmean.nc');
slpadf(:,:,1:4) = []; slpadf(:,:,end-3:end) = [];
slp = double(slpadf);
%%
[flnsadf] = DATAadf('FLNS','E:\CESM-post\JC-data\LE\FLNS\historical\his_FLNS_ensmean.nc');
flnsadf(:,:,1:4) = []; flnsadf(:,:,end-3:end) = [];
flns = double(flnsadf);

[fsnsadf] = DATAadf('FSNS','E:\CESM-post\JC-data\LE\FSNS\historical\his_FSNS_ensmean.nc');
fsnsadf(:,:,1:4) = []; fsnsadf(:,:,end-3:end) = [];
fsns = double(fsnsadf);

[lhflxadf] = DATAadf('LHFLX','E:\CESM-post\JC-data\LE\LHFLX\historical\his_LHFLX_ensmean.nc');
lhflxadf(:,:,1:4) = []; lhflxadf(:,:,end-3:end) = [];
lhflx = double(lhflxadf);

[shflxadf] = DATAadf('SHFLX','E:\CESM-post\JC-data\LE\SHFLX\historical\his_SHFLX_ensmean.nc');
shflxadf(:,:,1:4) = []; shflxadf(:,:,end-3:end) = [];
shflx = double(shflxadf);
%% Z200 Pac ensmean
filename1 = 'E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\Z3_192001-200512_ensmean.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
Z3rawPAC = ncread(filename1,'Z');
lonData = ncread(filename1,'lon');
latData = ncread(filename1,'lat');
levData = ncread(filename1,'level'); % unit: cm -> m
% minus LE ensemble mean
filename3 = 'F:\ensmean_CESM_LE_Z3_hist_interpolated.nc';
ncid=netcdf.open(filename3,'NOWRITE');
ncdisp(filename3);
Z3rawLE = ncread(filename3,'Z');
Z3pacd = Z3rawPAC-Z3rawLE; % detrend
Z3pacda = Z3pacd-mean(Z3pacd,4); % anomaly
Z200ad = permute(Z3pacda(:,:,10,:),[1 2 4 3]);
% detrend and 8-yr filter
dT = 1;  cf = 1/8;
Z200ad = Z200ad;
[d1 d2 d3] = size(Z200ad);
Z200adr = reshape(permute(Z200ad,[3 1 2]),[d3,d1*d2]);
clear Z200adf
parfor i = 1:size(Z200adr,2);
    Z200adf(:,i) = lanczosfilter(Z200adr(:,i),dT,cf,[],'low'); % 8 year filtered
end
Z200adf = permute(reshape(Z200adf,[d3,d1,d2]),[2 3 1]);  % -5~707m, 8-year filtered
Z200adf(:,:,1:4) = []; Z200adf(:,:,end-3:end) = []; 
filename1 = 'E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\Z3_192001-200512_ensmean.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
Z3rawPAC = ncread(filename1,'Z');
lonData = ncread(filename1,'lon');
latData = ncread(filename1,'lat');
levData = ncread(filename1,'level'); % unit: cm -> m
% minus LE ensemble mean
filename3 = 'F:\ensmean_CESM_LE_Z3_hist_interpolated.nc';
ncid=netcdf.open(filename3,'NOWRITE');
ncdisp(filename3);
Z3rawLE = ncread(filename3,'Z');
Z3pacd = Z3rawPAC-Z3rawLE; % detrend
Z3pacda = Z3pacd-mean(Z3pacd,4); % anomaly
Z200ad = permute(Z3pacda(:,:,10,:),[1 2 4 3]);
% detrend and 8-yr filter
dT = 1;  cf = 1/8;
Z200ad = Z200ad;
[d1 d2 d3] = size(Z200ad);
Z200adr = reshape(permute(Z200ad,[3 1 2]),[d3,d1*d2]);
clear Z200adf
parfor i = 1:size(Z200adr,2);
    Z200adf(:,i) = lanczosfilter(Z200adr(:,i),dT,cf,[],'low'); % 8 year filtered
end
Z200adf = permute(reshape(Z200adf,[d3,d1,d2]),[2 3 1]);  % -5~707m, 8-year filtered
Z200adf(:,:,1:4) = []; Z200adf(:,:,end-3:end) = []; 
%% mixed layer depth
filename1 = 'E:\CESM-post\PAC-pacemaker\HMXLensmean_1920-2005.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
HMXLrawPAC = ncread(filename1,'HMXL');
% minus LE ensemble mean
filename3 = 'E:\CESM-post\LE\HMXL\HMXLensmean_1920-2005.nc';
ncid=netcdf.open(filename3,'NOWRITE');
ncdisp(filename3);
HMXLrawLE = ncread(filename3,'HMXL');
HMXLpacd = HMXLrawPAC-HMXLrawLE; % detrend
HMXLpacda = HMXLpacd-mean(HMXLpacd,3); % anomaly
% detrend and 8-yr filter
dT = 1;  cf = 1/8;
HMXLad = HMXLpacda;
[d1 d2 d3] = size(HMXLad);
HMXLadr = reshape(permute(HMXLad,[3 1 2]),[d3,d1*d2]);
clear HMXLadf
parfor i = 1:size(HMXLadr,2);
    HMXLadf(:,i) = lanczosfilter(HMXLadr(:,i),dT,cf,[],'low'); % 8 year filtered
end
HMXLadf = permute(reshape(HMXLadf,[d3,d1,d2]),[2 3 1]);  % -5~707m, 8-year filtered
HMXLadf(:,:,1:4) = []; HMXLadf(:,:,end-3:end) = []; 



%%

function [VARsubadf] = DATAadf(varname,filename3)
% varname 
% filename3 is LE ensmean of the variable
    addpath D:\1_matlab\help;
    filename1 = ['E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\',varname,'_192001-200512_ensmean.nc'];
    ncid=netcdf.open(filename1,'NOWRITE');
    ncdisp(filename1);
    VARMMM1 = ncread(filename1,varname);
    VARMMM1(:,:,end) = []; % remove the last year (only January)
    filename2 = ['E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\',varname,'_200601-201312_ensmean.nc'];
    ncid=netcdf.open(filename2,'NOWRITE');
    ncdisp(filename2);
    VARMMM2 = ncread(filename2,varname);
    VARMMM2(:,:,end) = []; % remove the last year (only January)
    VARMMM = cat(3,VARMMM1,VARMMM2);
    % minus LE ensemble mean
    ncid=netcdf.open(filename3,'NOWRITE');
    ncdisp(filename3);
    VARLEMsub1 = ncread(filename3,varname);
    VARLEMsub1(:,:,end) = []; % remove the last year (only January)
    filename4 = ['E:\CESM-post\JC-data\LE\',varname,'\rcp85\',varname,'_rcp85_ensmean.nc'];
    ncid=netcdf.open(filename4,'NOWRITE');
    ncdisp(filename4);
    VARLEMsub2 = ncread(filename4,varname);
    VARLEMsub = cat(3,VARLEMsub1,VARLEMsub2(:,:,1:8));
    VARMMMsubd = VARMMM-VARLEMsub; % detrend
    VARMMMsubda = VARMMMsubd - mean(VARMMMsubd,3); % anomaly
    % detrend and 8-yr filter
    dT = 1;  cf = 1/8;
    VARda = VARMMMsubda;
    [d1 d2 d3] = size(VARda);
    VARsubar = reshape(permute(VARda,[3 1 2]),[d3,d1*d2]);
    parfor i = 1:size(VARsubar,2);
        VARsubadf(:,i) = lanczosfilter(VARsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    VARsubadf = permute(reshape(VARsubadf,[d3,d1,d2]),[2 3 1]);  % -5~707m, 8-year filtered
end

