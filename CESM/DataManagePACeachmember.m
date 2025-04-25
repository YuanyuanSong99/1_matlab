% TPIensm = TPIf; PC1ensm = PC1z; DIensm = DI_rawlong;
num = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
for s = 1:20;
    s
% each member data
% TEMP
% addpath(genpath('D:\1_matlab\help'));
% filename1 = ['E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\multi-level\TEMP\b.e11.B20TRLENS.f09_g16.SST.restoring.ens',num{s},'.pop.h.TEMP.192001-200512.nc'];
% ncid=netcdf.open(filename1,'NOWRITE');
% ncdisp(filename1);
% TemMMM1 = ncread(filename1,'TEMP');
% TemMMM1(:,:,:,end) = []; % remove the last year (only January)
% filename2 = ['E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\multi-level\TEMP\b.e11.BRCP85LENS.f09_g16.SST.restoring.ens',num{s},'.pop.h.TEMP.200601-201312.nc'];
% ncid=netcdf.open(filename2,'NOWRITE');
% ncdisp(filename2);
% TemMMM2 = ncread(filename2,'TEMP');
% TemMMM2(:,:,:,end) = []; % remove the last year (only January)
% lonData = ncread(filename1,'lon');
% latData = ncread(filename1,'lat');
% depthData = ncread(filename1,'z_t')/100; % unit: cm -> m
% TemMMM = cat(4,TemMMM1,TemMMM2);
% TemMMMsub = TemMMM(:,1:90,1:37,:); % southern hemisphere; up 707m 
% minus LE ensemble mean
% filename3 = 'E:\CESM-post\JC-data\LE\Temperature\his_ensmean.nc';
% ncid=netcdf.open(filename3,'NOWRITE');
% ncdisp(filename3);
% TemLEMsub1 = ncread(filename3,'TEMP');
% TemLEMsub1(:,:,:,end) = [];
% filename4 = 'E:\CESM-post\JC-data\LE\Temperature\rcp85_sub_ensmean.nc';
% ncid=netcdf.open(filename4,'NOWRITE');
% ncdisp(filename4);
% TemLEMsub2 = ncread(filename4,'TEMP');
% TemLEMsub = cat(4,TemLEMsub1,TemLEMsub2(:,1:90,:,:));
% TemMMMsubd = TemMMMsub-TemLEMsub; % detrend
% TemMMMsubda = TemMMMsubd-mean(TemMMMsubd,4); % anomaly
% % detrend and 8-yr filter
% dT = 1;  cf = 1/8;
% Temda = TemMMMsubda;
% [d1 d2 d3 d4] = size(Temda);
% Temsubar = reshape(permute(Temda,[4 1 2 3]),[d4,d1*d2*d3]);
% clear Temsubadf_long
% parfor i = 1:size(Temsubar,2);
%     Temsubadf_long(:,i) = lanczosfilter(Temsubar(:,i),dT,cf,[],'low'); % 8 year filtered
% end
% Temsubadf_long = permute(reshape(Temsubadf_long,[d4,d1,d2,d3]),[2 3 4 1]);  % -5~707m, 8-year filtered
% Temsubadf = Temsubadf_long(:,:,:,5:end-4); 
% SST  lon*lat*1*time
filename1 = ['E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\single-level\0',num{s},'\b.e11.B20TRLENS.f09_g16.SST.restoring.ens',num{s},'.pop.h.SST.192001-200512.nc'];
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
TemMMM1 = ncread(filename1,'SST');
TemMMM1(:,:,:,end) = []; % remove the last year (only January)
filename2 = ['E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\single-level\0',num{s},'\b.e11.BRCP85LENS.f09_g16.SST.restoring.ens',num{s},'.pop.h.SST.200601-201312.nc'];
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
TemMMM2 = ncread(filename2,'SST');
TemMMM2(:,:,:,end) = []; % remove the last year (only January)
sstMMM = cat(4,TemMMM1,TemMMM2);
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
clear sstadf_long
% parfor i = 1:size(sstar,2);
%     sstadf_long(:,i) = lanczosfilter(sstar(:,i),dT,cf,[],'low'); % 8 year filtered
% end
% sstadf_long = permute(reshape(sstadf_long,[d4,d1,d2,d3]),[2 3 1 4]);  % -5~707m, 8-year filtered
% sstadf = sstadf_long(:,:,5:end-4); 

% calculate PC1 DI IPO
% EOF
% lats = 1:90;
% dweit = depthData(2:37)-depthData(1:36); % depth weight
% Tsub = permute(nansum(Temsubadf(:,lats,1,:)*5+Temsubadf(:,lats,2:37,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m
% [eof_maps,pc,expvar]=eofnan(Tsub);  
% EOF1z = -eof_maps(:,:,1)*std(pc(1,:)); 
% PC1z = -pc(1,:)/std(pc(1,:)); 
% % DI index (Pac_T - IO&Atl_T)
% Tsub_long = permute(nansum(Temsubadf_long(:,lats,1,:)*5+Temsubadf_long(:,lats,2:37,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m
% lats = [29:44];
% Tsub_r = cat(1,Tsub_long(160:300,:,:),Tsub_long(301:360,:,:),Tsub_long(1:159,:,:));
% lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
% [spacz_long spac_long] = areamean(Tsub_r,1:141,lats,latData); 
% [siaz_long sia_long] = areamean(Tsub_r,142:360,lats,latData); 
% DI_rawlong = spac_long-sia_long;
% DI_long = zscore(spac_long-sia_long);
% DI = DI_long(5:end-4);
% corrcoef(PC1z,DI)
% TPI  index
sstd = permute(sstad,[1 2 4 3]);
clear TPIfz AMOf8z
% TPI index ---------------------------------------------------------------
[ts1_zs ts1] = areamean(sstd,141:216,116:136,latData); % 25N-45N,140E-145W
[ts2_zs ts2] = areamean(sstd,171:271,80:101,latData); % 10S-10N,170E-90W
[ts3_zs ts3] = areamean(sstd,151:201,40:75,latData); % 50S-15S,150E-160W
TPI1 = ts2-(ts1+ts3)/2; % unfiltered IPO index///
dT = 1; % interval
cf = 1/8;
TPIf = lanczosfilter(TPI1,dT,cf,[],'low'); % 8 year filtered IPO index
TPIfz_long = zscore(TPIf);
TPIfz = TPIfz_long(5:end-4);
%
    TPIpac(:,s) = TPIf; 
    % PC1pac(s,:) = PC1z; DIpac(:,s) = DI_rawlong;
    clearvars -except TPIensm PC1ensm DIensm TPIpac PC1pac DIpac num lonData...
        latData depthData

end




%% single level
[tauxadf] = DATAadf('TAUX',num,'pop.h','E:\data-post\CESM\LE\TAU\historical\his_taux_ensmean.nc');
[tauyadf] = DATAadf('TAUY',num,'pop.h','E:\data-post\CESM\LE\TAU\historical\his_tauy_ensmean.nc');
tauxadf(:,:,1:4) = []; tauxadf(:,:,end-3:end) = [];
taux = double(tauxadf)*10; % unit: dyn/cm2 -> N/m2
tauyadf(:,:,1:4) = []; tauyadf(:,:,end-3:end) = [];
tauy = double(tauyadf)*10; % unit: dyn/cm2 -> N/m2
%
[slpadf] = DATAadf('PSL',num,'cam.h0','E:\data-post\CESM\LE\psl\historical\his_psl_ensmean.nc');
slpadf(:,:,1:4) = []; slpadf(:,:,end-3:end) = [];
slp = double(slpadf);




function [VARsubadf] = DATAadf(varname,num,modules,filename3)
% varname 
% filename3 is LE ensmean of the variable
    addpath E:\1_matlab\help;
    filename1 = ['E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\single-level\0',num,'\b.e11.B20TRLENS.f09_g16.SST.restoring.ens',num,'.',modules,'.',varname,'.192001-200512.nc'];
    ncid=netcdf.open(filename1,'NOWRITE');
    ncdisp(filename1);
    VARMMM1 = ncread(filename1,varname);
    VARMMM1(:,:,end) = []; % remove the last year (only January)
    filename2 = ['E:\CESM-post\JC-data\PAC-pacemaker\yr_1x1\single-level\0',num,'\b.e11.BRCP85LENS.f09_g16.SST.restoring.ens',num,'.',modules,'.',varname,'.200601-201312.nc'];
    ncid=netcdf.open(filename2,'NOWRITE');
    ncdisp(filename2);
    VARMMM2 = ncread(filename2,varname);
    VARMMM2(:,:,end) = []; % remove the last year (only January)
    VARMMM = cat(3,VARMMM1,VARMMM2);
    % minus LE climotology
    % filename3 = 'E:\CESM-post\JC-data\LE\Temperature\his_ensmean.nc';
    ncid=netcdf.open(filename3,'NOWRITE');
    ncdisp(filename3);
    VARLEMsub = ncread(filename3,varname);
    VARLEMsub(:,:,end) = []; % remove the last year (only January)
    VARLEMsubcli = mean(VARLEMsub,3);
    VARMMMsuba = VARMMM-VARLEMsubcli; % anomaly
    % detrend and 8-yr filter
    dT = 1;  cf = 1/8;
    VARda = VARMMMsuba;
    [d1 d2 d3] = size(VARda);
    VARsubar = reshape(permute(VARda,[3 1 2]),[d3,d1*d2]);
    parfor i = 1:size(VARsubar,2);
        VARsubad(:,i) = detrend(VARsubar(:,i)); % detrend
        VARsubadf(:,i) = lanczosfilter(VARsubad(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    VARsubadf = permute(reshape(VARsubadf,[d3,d1,d2]),[2 3 1]);  % -5~707m, 8-year filtered
end
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end
