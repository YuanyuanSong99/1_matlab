clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
z_t = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','z_t'); % cm -> m
DXU = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','DXU'); % cm -> m
DXT = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','DXT');
DYU = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','DYU');
DYT = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','DYT');
DZU = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','dz');
DZT = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','dz');
TLONG = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','TLONG');
TLAT = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','TLAT');
TAREA = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','TAREA');
ULONG = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','ULONG');
ULAT = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','ULAT');
% TEMP
TEMP1 = ncread('/Volumes/CESM-post2/TEMPensmean_1920-2005.nc','TEMP'); 
TEMP2 = ncread('/Volumes/CESM-post2/TEMPensmean_2006-2080.nc','TEMP');  
TEMPmon = cat(4,TEMP1(:,:,:,480+2:end),TEMP2(:,:,:,1:181));
clear TEMP1 TEMP2
TEMP = squeeze(mean(reshape(TEMPmon,320,384,60,12,732/12),4)); % yearly
clear TEMPmon
TEMP_tend = diff(TEMP,1,4)/(60*60*24*365); % second   degC/s
% save('/Volumes/Togo4T/1_matlab/MatData/CESM/tendency/TEMP_tend.mat','TEMP_tend');
clear TEMP
% adv UT
UET1 = ncread('/Volumes/CESM-post2/UETensmean_1920-2005.nc','UET'); % degC/s
UET2 = ncread('/Volumes/CESM-post2/UETensmean_2006-2080.nc','UET'); % degC/s
UETmon = cat(4,UET1(:,:,:,480+2:end),UET2(:,:,:,1:181));
clear UET1 UET2
UET = squeeze(mean(reshape(UETmon,320,384,60,12,732/12),4)); % yearly
clear UETmon
[d1 d2 d3 d4] = size(UET);
DZTr = repmat(permute(DZT,[3 2 1]),[d1 d2 1]);
DXTr = repmat(DXT,[1 1 d3]);
DYTr = repmat(DYT,[1 1 d3]);
VOL = DZTr.*DXTr.*DYTr; % volume of T grid
VOLex = cat(1,VOL,VOL(1,:,:)); 
UETe = cat(1,UET,UET(1,:,:,:)); % extend to calculate difference
clear UET
dutdx = -diff(UETe.*VOLex,1,1)./VOL;
clear UETe  
% adv VT
VNT1 = ncread('/Volumes/CESM-post2/VNTensmean_1920-2005.nc','VNT'); % degC/s
VNT2 = ncread('/Volumes/CESM-post2/VNTensmean_2006-2080.nc','VNT'); % degC/s
VNTmon = cat(4,VNT1(:,:,:,480+2:end),VNT2(:,:,:,1:181));
clear VNT1 VNT2
VNT = squeeze(mean(reshape(VNTmon,320,384,60,12,732/12),4)); % yearly
clear VNTmon
[d1 d2 d3 d4] = size(VNT);
VOLey = cat(2,VOL,VOL(:,end,:));  
VNTe = cat(2,VNT,VNT(:,end,:,:)); % extend to calculate difference
clear VNT
dvtdy = -diff(VNTe.*VOLey,1,2)./VOL;
clear VNTe
% adv WT
WTT1 = ncread('/Volumes/CESM-post2/WTTensmean_1920-2005.nc','WTT'); % degC/s
WTT2 = ncread('/Volumes/CESM-post2/WTTensmean_2006-2080.nc','WTT'); % degC/s
WTTmon = cat(4,WTT1(:,:,:,480+2:end),WTT2(:,:,:,1:181));
clear WTT1 WTT2
WTT = squeeze(mean(reshape(WTTmon,320,384,60,12,732/12),4)); % yearly
clear WTTmon
[d1 d2 d3 d4] = size(WTT);
VOLez = cat(3,VOL,VOL(:,:,end));
WTTe = cat(3,WTT,WTT(:,:,end,:)); % extend to calculate difference
clear WTT
dztdz = diff(WTTe.*VOLez,1,3)./VOL;
clear WTTe

TOT_ADV = dutdx + dvtdy+ dztdz;
% save('/Volumes/Togo4T/1_matlab/MatData/CESM/tendency/dutdx.mat','dutdx');
% save('/Volumes/Togo4T/1_matlab/MatData/CESM/tendency/dvtdy.mat','dvtdy');
% save('/Volumes/Togo4T/1_matlab/MatData/CESM/tendency/dztdz.mat','dztdz');
clear dutdx dvtdy dztdz
%% vertical mixing
DIA_IMPVF_TEMP1 = ncread('/Volumes/CESM-post2/DIA_IMPVF_TEMPensmean_1920-2005.nc','DIA_IMPVF_TEMP'); % degC/s
DIA_IMPVF_TEMP2 = ncread('/Volumes/CESM-post2/DIA_IMPVF_TEMPensmean_2006-2080.nc','DIA_IMPVF_TEMP'); % degC/s
DIA_IMPVF_TEMPmon = cat(4,DIA_IMPVF_TEMP1(:,:,:,480+2:end),DIA_IMPVF_TEMP2(:,:,:,1:181));
clear DIA_IMPVF_TEMP1 DIA_IMPVF_TEMP2
DIA_IMPVF_TEMP = squeeze(mean(reshape(DIA_IMPVF_TEMPmon,320,384,60,12,732/12),4)); % yearly
clear DIA_IMPVF_TEMPmon
DIA_IMPVF_TEMPe = cat(3,DIA_IMPVF_TEMP,DIA_IMPVF_TEMP(:,:,end,:));
clear DIA_IMPVF_TEMP
DIA = -diff(DIA_IMPVF_TEMPe.*TAREA,1,3)./VOL;
clear DIA_IMPVF_TEMPe 
% set surface flux at 0th layer
TAREA = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','TAREA');
% KMT = 100*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','KMT'); 
hflux_factor = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','hflux_factor');
SHF1 = ncread('/Volumes/CESM-post2/SHFensmean_1920-2005.nc','SHF'); % degC/s
SHF2 = ncread('/Volumes/CESM-post2/SHFensmean_2006-2080.nc','SHF'); % degC/s
SHFmon = cat(3,SHF1(:,:,480+2:end),SHF2(:,:,1:181));
clear SHF1 SHF2
SHF = squeeze(mean(reshape(SHFmon,320,384,12,732/12),3)); % yearly
clear SHFmon

SHF_QSW1 = ncread('/Volumes/CESM-post2/SHF_QSWensmean_1920-2005.nc','SHF_QSW'); % degC/s
SHF_QSW2 = ncread('/Volumes/CESM-post2/SHF_QSWensmean_2006-2080.nc','SHF_QSW'); % degC/s
SHF_QSWmon = cat(3,SHF_QSW1(:,:,480+2:end),SHF_QSW2(:,:,1:181));
clear SHF_QSW1 SHF_QSW2
SHF_QSW = squeeze(mean(reshape(SHF_QSWmon,320,384,12,732/12),3)); % yearly
clear SHF_QSWmon

SRF_TEMP_FLUX = (SHF-SHF_QSW) * hflux_factor;
DIA(:,:,1,:) = (permute(SRF_TEMP_FLUX .* TAREA,[1 2 4 3]) - DIA(:,:,1,:) .* TAREA) ./ VOL(:,:,1);
clear SHF_QSW* SRF_TEMP_FLUX

KPP_SRC_TEMP1 = ncread('/Volumes/CESM-post2/KPP_SRC_TEMPensmean_1920-2005.nc','KPP_SRC_TEMP'); % degC/s
KPP_SRC_TEMP2 = ncread('/Volumes/CESM-post2/KPP_SRC_TEMPensmean_2006-2080.nc','KPP_SRC_TEMP'); % degC/s
KPP_SRC_TEMPmon = cat(4,KPP_SRC_TEMP1(:,:,:,480+2:end),KPP_SRC_TEMP2(:,:,:,1:181));
clear KPP_SRC_TEMP1 KPP_SRC_TEMP2
KPP_SRC_TEMP = squeeze(mean(reshape(KPP_SRC_TEMPmon,320,384,60,12,732/12),4)); % yearly
clear KPP_SRC_TEMPmon

VDIF = DIA + KPP_SRC_TEMP;
clear DIA KPP_SRC_TEMP
% save('/Volumes/Togo4T/1_matlab/MatData/CESM/tendency/VDIF.mat','VDIF');
%%
% horizontal diffusion
HDIFE_TEMP1 = ncread('/Volumes/CESM-post2/HDIFE_TEMPensmean_1920-2005.nc','HDIFE_TEMP'); % degC/s
HDIFE_TEMP2 = ncread('/Volumes/CESM-post2/HDIFE_TEMPensmean_2006-2080.nc','HDIFE_TEMP'); % degC/s
HDIFE_TEMPmon = cat(4,HDIFE_TEMP1(:,:,:,480+2:end),HDIFE_TEMP2(:,:,:,1:181));
clear HDIFE_TEMP1 HDIFE_TEMP2
HDIFE_TEMP = squeeze(mean(reshape(HDIFE_TEMPmon,320,384,60,12,732/12),4)); % yearly
clear HDIFE_TEMPmon
HDIFE_TEMPe = cat(1,HDIFE_TEMP,HDIFE_TEMP(1,:,:,:));
clear HDIFE_TEMP
dhdfe = diff(HDIFE_TEMPe.*VOLex,1,1)./VOL;
clear HDIFE_TEMPe
%
HDIFN_TEMP1 = ncread('/Volumes/CESM-post2/HDIFN_TEMPensmean_1920-2005.nc','HDIFN_TEMP'); % degC/s
HDIFN_TEMP2 = ncread('/Volumes/CESM-post2/HDIFN_TEMPensmean_2006-2080.nc','HDIFN_TEMP'); % degC/s
HDIFN_TEMPmon = cat(4,HDIFN_TEMP1(:,:,:,480+2:end),HDIFN_TEMP2(:,:,:,1:181));
clear HDIFN_TEMP1 HDIFN_TEMP2
HDIFN_TEMP = squeeze(mean(reshape(HDIFN_TEMPmon,320,384,60,12,732/12),4)); % yearly
clear HDIFN_TEMPmon
HDIFN_TEMPe = cat(2,HDIFN_TEMP,HDIFN_TEMP(:,end,:,:));
clear HDIFN_TEMP
dhdfn = diff(HDIFN_TEMPe.*VOLey,1,2)./VOL;
clear HDIFN_TEMPe
HDIF = dhdfe + dhdfn;
clear dhdfe dhdfn
% save('/Volumes/Togo4T/1_matlab/MatData/CESM/tendency/HDIF.mat','HDIF');

%% Solar penetration
QSW_3D1 = ncread('/Volumes/CESM-post2/QSW_3Densmean_1920-2005.nc','QSW_3D'); % degC/s
QSW_3D2 = ncread('/Volumes/CESM-post2/QSW_3Densmean_2006-2080.nc','QSW_3D'); % degC/s
QSW_3Dmon = cat(4,QSW_3D1(:,:,:,480+2:end),QSW_3D2(:,:,:,1:181));
clear QSW_3D1 QSW_3D2
QSW_3D = squeeze(mean(reshape(QSW_3Dmon,320,384,60,12,732/12),4)); % yearly
clear QSW_3Dmon
QSW_3De = cat(3,QSW_3D,QSW_3D(:,:,end,:));
clear QSW_3D
QSW = -diff(QSW_3De * hflux_factor,1,3)./permute(DZT,[3 2 1]);
clear QSW_3De 
% save('/Volumes/Togo4T/1_matlab/MatData/CESM/tendency/QSW.mat','QSW');

% compare
allvar = QSW + HDIF + VDIF + TOT_ADV;
%%
close all
figure('Position',[100 100 400 800])
plot(squeeze(TEMP_tend(161,187,1:34,1)),-z_t(1:34)/100,'k-','LineWidth',1)
hold on
plot(squeeze(allvar(161,187,1:34,1)),-z_t(1:34)/100,'g--','LineWidth',1)
plot(squeeze(QSW(161,187,1:34,1)),-z_t(1:34)/100,'r--','LineWidth',1)
plot(squeeze(VDIF(161,187,1:34,1)),-z_t(1:34)/100,'r-','LineWidth',1)
plot(squeeze(HDIF(161,187,1:34,1)),-z_t(1:34)/100,'b--','LineWidth',1)
plot(squeeze(TOT_ADV(161,187,1:34,1)),-z_t(1:34)/100,'b-','LineWidth',1)
plot(zeros(1,47),-z_t(1:47)/100,'k-','LineWidth',1.5)
ylabel('depth (m)','FontSize',14)
xlabel('TEMP tendency')
set(gca,'ylim',[-500,0],'xlim',[-5,5]*10^-6,'fontsize',14)
legend('TEMP tend','sum','QSW','VDIF','HDIF','total adv','Location','southeast')
% saveas(gca,'/Users/yysong/Desktop/figures/CESM/monthly/20240611/ens_TEMPtend_0_140E_unweited.png')
%% 
close all
figure('Position',[100 100 400 800])
plot(squeeze(TOT_ADV(161,187,1:34,1)),-z_t(1:34)/100,'b-','LineWidth',1)
hold on
plot(squeeze(dztdz(161,187,1:34,1)),-z_t(1:34)/100,'r--','LineWidth',1)
plot((squeeze(dutdx(161,187,1:34,1))+squeeze(dvtdy(161,187,1:34,1))),-z_t(1:34)/100,'r-','LineWidth',1)
plot(zeros(1,47),-z_t(1:47)/100,'k-','LineWidth',1.5)
ylabel('depth (m)','FontSize',14)
xlabel('advection')
set(gca,'ylim',[-500,0],'xlim',[-8,8]*10^-6,'fontsize',14)
legend('total adv','adv_w','adv_u+_v','Location','southeast')
% saveas(gca,'/Users/yysong/Desktop/figures/CESM/monthly/20240611/ens_advection_0_140E.png')
%%
close all
plot(squeeze(dutdx(161,187,1:47,1)),-z_t(1:47)/100,'r-','LineWidth',1.5)
hold on
plot(squeeze(dvtdy(161,187,1:47,1)),-z_t(1:47)/100,'r--','LineWidth',1.5)
plot((squeeze(dutdx(161,187,1:47,1))+squeeze(dvtdy(161,187,1:47,1))),-z_t(1:47)/100,'b-','LineWidth',1.5)
plot(zeros(1,47),-z_t(1:47)/100,'k-','LineWidth',1.5)
ylabel('depth (m)','FontSize',14)
xlabel('advection')
set(gca,'ylim',[-2000,0],'fontsize',14)
legend('adv_u','adv_v','adv_u+adv_v','Location','southeast')
% saveas(gca,'/Users/yysong/Desktop/figures/CESM/monthly/20240611/ens_advectionUV_0_140E.png')
%%
close all
figure(1)
contourf(TEMP_tend(:,:,10,1)')
caxis([-1.5,1.5]*10^-5)
colorbar
figure(2)
contourf(allvar(:,:,10,1)')
caxis([-1.5,1.5]*10^-5)
colorbar
%%
intlev = 37; % 707 m
splon = 169; % 150 E
inlats = 46:84; % 55.1805 S :  34.8794 S
% [tdp tdia tdall] = caltend(TEMP_tend,VOL,splon,inlats,intlev);
[tdp tdia tdall] = caltend(TEMP_tend.*permute(DZT,[3 2 1])/100,VOL,splon,inlats,intlev);
[allp allia allall] = caltend(allvar,VOL,splon,inlats,intlev);
[focp focia focall] = caltend(VDIF+QSW,VOL,splon,inlats,intlev);
[vdifp vdifia vdifall] = caltend(VDIF,VOL,splon,inlats,intlev);
[qswp qswia qswall] = caltend(QSW,VOL,splon,inlats,intlev);
[advp advia advall] = caltend(TOT_ADV,VOL,splon,inlats,intlev);
[hdifp hdifia hdifall] = caltend(HDIF,VOL,splon,inlats,intlev);

%%
close all
figure(1)
plot([2:61],cumsum(tdall)/10^21,'k','LineWidth',1.5)
hold on
plot(cumsum(allall)/10^21,'k--','LineWidth',1.5)
plot(cumsum(focall)/10^21,'r-','LineWidth',1.5)
plot(cumsum(advall)/10^21,'b-','LineWidth',1.5)
plot(cumsum(hdifall)/10^21,'g-','LineWidth',1.5)
plot(zeros(1,61),'k')
% set(gca,'ylim',[-100,200],'xlim',[1,61],'fontsize',14)
set(gca,'xtick',[1:10:61],'XTickLabel',[1960:10:2020])
ylabel('ZJ')
legend('tendency','sum','VDIF+QSW','DIV','HDIF','Location','northwest')
% saveas(gca,'/Users/yysong/Desktop/figures/CESM/monthly/20240611/ens_tendency_35-55S.png')
%%
close all
figure(1)
plot([2:61],tdall,'k','LineWidth',1.5)
hold on
plot(allall,'k--','LineWidth',1.5)
plot(qswall,'r-','LineWidth',1.5)
plot(advall,'b-','LineWidth',1.5)
plot(hdifall,'g-','LineWidth',1.5)
plot(vdifall,'g--','LineWidth',1.5)
plot(zeros(1,61),'k')
set(gca,'ylim',[-7,10]*10^21,'xlim',[1,61],'fontsize',14)
set(gca,'xtick',[1:10:61],'XTickLabel',[1960:10:2020])
legend('tendency','sum','QSW','DIV','HDIF','VDIF','Location','northwest')
% saveas(gca,'/Users/yysong/Desktop/figures/CESM/monthly/20240611/ens_times_35-55S.png')
%%
figure(2)
plot([2:61],cumsum(tdp)/10^21,'k','LineWidth',1.5)
hold on
plot(cumsum(allp)/10^21,'k--','LineWidth',1.5)
plot(cumsum(focp)/10^21,'r-','LineWidth',1.5)
plot(cumsum(advp)/10^21,'b-','LineWidth',1.5)
plot(cumsum(hdifp)/10^21,'g-','LineWidth',1.5)
plot(zeros(1,61),'k')
set(gca,'ylim',[-50,120],'xlim',[1,61],'fontsize',14)
set(gca,'xtick',[1:10:61],'XTickLabel',[1960:10:2020])
ylabel('ZJ')
legend('tendency','sum','VDIF+QSW','DIV','HDIF','Location','northwest')
saveas(gca,'/Users/yysong/Desktop/figures/CESM/monthly/20240611/ens_tendency_Pac.png')
%%
close all
figure(1)
plot([2:61],tdp,'k','LineWidth',1.5)
hold on
plot(allp,'k--','LineWidth',1.5)
plot(qswp,'r-','LineWidth',1.5)
plot(advp,'b-','LineWidth',1.5)
plot(hdifp,'g-','LineWidth',1.5)
plot(vdifp,'g--','LineWidth',1.5)
plot(zeros(1,61),'k')
set(gca,'ylim',[-7,10]*10^21,'xlim',[1,61],'fontsize',14)
set(gca,'xtick',[1:10:61],'XTickLabel',[1960:10:2020])
legend('tendency','sum','QSW','DIV','HDIF','VDIF','Location','northwest')
% saveas(gca,'/Users/yysong/Desktop/figures/CESM/monthly/20240611/ens_times_35-55S.png')
%%
figure(3)
plot([2:61],cumsum(tdia)/10^21,'k','LineWidth',1.5)
hold on
plot(cumsum(allia)/10^21,'k--','LineWidth',1.5)
plot(cumsum(focia)/10^21,'r-','LineWidth',1.5)
plot(cumsum(advia)/10^21,'b-','LineWidth',1.5)
plot(cumsum(hdifia)/10^21,'g-','LineWidth',1.5)
plot(zeros(1,61),'k')
set(gca,'ylim',[-50,120],'xlim',[1,61],'fontsize',14)
set(gca,'xtick',[1:10:61],'XTickLabel',[1960:10:2020])
ylabel('ZJ')
legend('tendency','sum','VDIF+QSW','DIV','HDIF','Location','northwest')
saveas(gca,'/Users/yysong/Desktop/figures/CESM/monthly/20240611/ens_tendency_IO+Atl.png')
%%
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
% hfc = ro*cp* (SHF*hflux_factor.*TAREA- squeeze(DIA(:,:,1,:) .* TAREA))./VOL(:,:,1) .*VOL(:,:,1);
hfc = ro*cp* (SHF*hflux_factor.*TAREA);
hfc_r = cat(1,hfc(splon:320,inlats,:),hfc(1:splon-1,inlats,:));
hfcp = squeeze(nansum(nansum(hfc_r(1:125,:,:),1),2));
hfcia = squeeze(nansum(nansum(hfc_r(126:320,:,:),1),2));
hfca = squeeze(nansum(nansum(hfc_r(126:205,:,:),1),2));
hfci = squeeze(nansum(nansum(hfc_r(206:320,:,:),1),2));
hfcall = squeeze(nansum(nansum(hfc_r,1),2));

clear foc*
focp = cumsum(hfcp);
focia = cumsum(hfcia);
foca = cumsum(hfca);
foci = cumsum(hfci);
focall = cumsum(hfcall);

%% 
close all
figure(1)
plot([2:61],cumsum(tdp)/10^21,'k','LineWidth',1.5)
hold on
plot(cumsum(focp+advp)/10^21,'k--','LineWidth',1.5)
plot(cumsum(focp)/10^21,'r-','LineWidth',1.5)
plot(cumsum(advp)/10^21,'b-','LineWidth',1.5)
% plot(cumsum(vdifp)/10^21,'b--','LineWidth',1.5)
% plot(cumsum(qswp)/10^21,'r--','LineWidth',1.5)

%%
figure(2)
plot([2:61],cumsum(tdia)/10^21,'k','LineWidth',1.5)
hold on
plot(cumsum(focia+advia)/10^21,'k--','LineWidth',1.5)
plot(cumsum(focia)/10^21,'r-','LineWidth',1.5)
plot(cumsum(advia)/10^21,'b-','LineWidth',1.5)
% plot([2:61],cumsum(tdia-focia(2:61))/10^21,'k','LineWidth',1.5)
% plot(cumsum(dia)/10^21,'g-','LineWidth',1.5)

function [utop utoia utoall] = caltend(var,VOL,splon,inlats,intlev)
    cp = 4096 % J (kg oC)
    ro = 1025 % kg m-3
    uto = ro*cp* var .* VOL;
    utor = cat(1,uto(splon:end,inlats,1:intlev,:),uto(1:splon-1,inlats,1:intlev,:)); % 1960-2005
    utop = squeeze(nansum(nansum(nansum(utor(1:125,:,:,:),1),2),3)); % 150E-70W; 289.625 ~ 70W
    utoia = squeeze(nansum(nansum(nansum(utor(126:end,:,:,:),1),2),3));
    utoa = squeeze(nansum(nansum(nansum(utor(126:205,:,:,:),1),2),3));
    utoi = squeeze(nansum(nansum(nansum(utor(206:end,:,:,:),1),2),3));
    utoall = squeeze(nansum(nansum(nansum(utor,1),2),3));

end
function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
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
