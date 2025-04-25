% Data 
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
filename1 = '/Users/yysong/Downloads/THETA.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
% Tempo = ncread(filename1,'THETA');
lonData = ncread(filename1,'lon_c'); % cm -> m
latData = ncread(filename1,'lat_c'); % cm -> m
depthData = ncread(filename1,'dep_c'); % cm -> m
%%
tend = ncread('/Users/yysong/Downloads/tend.nc','tend');
[pac_tend, ia_tend] = caltendency(cumsum(tend,4),lonData,latData,depthData,0);
clear tend
%%
theta = ncread('/Users/yysong/Downloads/THETA.nc','THETA');
[pac_temptend, ia_temptend] = caltendency(theta,lonData,latData,depthData,1);
clear theta
%%
forc = ncread('/Users/yysong/Downloads/forc.nc','forc');
[pac_forc, ia_forc] = caltendency(forc,lonData,latData,depthData,1);
clear forc
%%
adv = ncread('/Users/yysong/Downloads/adv.nc','adv');
[pac_adv, ia_adv] = caltendency(adv,lonData,latData,depthData,1);
clear adv
%%
dif = ncread('/Users/yysong/Downloads/dif.nc','dif');
[pac_dif, ia_dif] = caltendency(dif,lonData,latData,depthData,1);
clear dif
%% plot
close all;
Fig = figure('Position',[100 100 800 400])
plot(pac_temptend-pac_temptend(1),'k-');
hold on
plot(cumsum(pac_forc),'r-');
plot(cumsum(pac_adv),'b-');
plot(cumsum(pac_dif),'g-');
plot(cumsum(pac_dif)+cumsum(pac_adv)+cumsum(pac_forc),'color',[.06 1 1]);
legend('tendency','forcing','advection','diffusion','forc+adv+diff','location','northwest')
print(Fig,['/Users/yysong/Desktop/figures/OCCA2/monthly/20240512_tendencytest/tendency_Pacific.png'],'-dpng','-r300')

%% plot
close all;
Fig = figure('Position',[100 100 800 400])
plot(ia_temptend-ia_temptend(1),'k-');
hold on
plot(cumsum(ia_forc),'r-');
plot(cumsum(ia_adv),'b-');
plot(cumsum(ia_dif),'g-');
plot(cumsum(ia_dif)+cumsum(ia_adv)+cumsum(ia_forc),'color',[.06 1 1]);
legend('tendency','forcing','advection','diffusion','forc+adv+diff','location','northwest')
print(Fig,['/Users/yysong/Desktop/figures/OCCA2/monthly/20240512_tendencytest/tendency_Atlantic.png'],'-dpng','-r300')
 
 
function [spac_0, sia_0] = caltendency(Temp,lonData,latData,depthData,ix)
% integrated and cusum
% if ix = 0, calculate tendency; otherwise, calculate forcing and adv
    intlev = 26; intdep = 700; % 700m
    lats = [70:110]; % 55S-35S
    dweit = depthData(2:intlev-1)-depthData(1:intlev-2);
    levs = [1:intlev];
    Tempint = Temp(:,:,intlev-1,:)+(Temp(:,:,intlev,:)-Temp(:,:,intlev-1,:))*(intdep-depthData(intlev-1))/(depthData(intlev)-depthData(intlev-1));
    Tempz = cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:intlev-1,:).*(permute(dweit,[3 2 1])),Tempint*(intdep-depthData(intlev-1)));
    Tempzy = Tempz*111*1000/2; % 111 km / latitude
    clear dx
    for j = 1:length(latData);
        dx(j) = 0.5*sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km'); % 半个经度之间的距离
    end
    cp = 4096 % J (kg oC)
    ro = 1025 % kg m-3
    Tempzyx = Tempzy.*dx*1000; % m
    lonDatar = cat(1,lonData(661:720),lonData(1:660));
    Tzyxsub_r_0 = cat(1,Tempzyx(661:720,:,:,:),Tempzyx(1:660,:,:,:));
    spac_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:280,lats,levs,:),1),2),3)); % 150E-70W
    sia_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(281:720,lats,levs,:),1),2),3)); % 70W-150E
    
end
function varMMM = loadsurface(filename1,filename2,varstr)
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
varMMM1 = ncread(filename1,varstr);
varMMM1(:,:,87) = [];
lonData = ncread(filename1,'lon');
latData = ncread(filename1,'lat');

ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
varMMM2 = ncread(filename2,varstr);

varMMM = double(cat(3,varMMM1(:,:,1:86),varMMM2(:,:,1:75)));
end
function trd = trend_cal_3D(var)
% var is time*d1*d2
    x = [1:size(var,1)]';
    clear trd
    for i = 1:size(var,2);
        for j = 1:size(var,3);
            par=polyfit(x,var(:,i,j),1); % regression parameters
            trd(i,j) = par(1);
        end
    end
end
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end