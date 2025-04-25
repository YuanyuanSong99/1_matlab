% test
% Data 
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
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
filename2 = '/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_2006-2022.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
TemMMM2 = ncread(filename2,'TEMP');
TemMMMall = cat(4,TemMMM1,TemMMM2(:,:,1:47,:));
Temp = TemMMMall;
clear TemMMM1 TemMMM2
%% OHC
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
%------------------------- 0-700m -----------------------------------------
depthstr = '0-700m';
% original OHC time series
lats = [35:55]; % 55S-35S
intlev = 37; intdep = 700; % 700m
dweit = depthData(2:intlev-1)-depthData(1:intlev-2);
deps = [depthData(1:intlev-1);intdep];
dz = gradient(deps);
levs = [1:intlev];
VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
Tempint = Temp(:,:,intlev-1,:)+(Temp(:,:,intlev,:)-Temp(:,:,intlev-1,:))*(intdep-depthData(intlev-1))/(depthData(intlev)-depthData(intlev-1));
% Tempz = cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:intlev-1,:).*(permute(dweit,[3 2 1])),Tempint*(intdep-depthData(intlev-1)));
Tempz = cat(3,Temp(:,:,1:intlev-1,:),Tempint).*permute(dz,[2 3 1]);
VVint = VV(:,:,intlev-1,:)+(VV(:,:,intlev,:)-VV(:,:,intlev-1,:))*(intdep-depthData(intlev-1))/(depthData(intlev)-depthData(intlev-1));
VVz = cat(3,VV(:,:,1)*5,VV(:,:,2:intlev-1).*(permute(dweit,[3 2 1])),VVint*(intdep-depthData(intlev-1)));
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
spac_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:150,lats,levs,:),1),2),3));
spacV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:150,lats,levs),1),2),3));
sia_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(151:360,lats,levs,:),1),2),3));
siaV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(151:360,lats,levs),1),2),3));
sat_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(151:230,lats,levs,:),1),2),3));
satV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(151:230,lats,levs),1),2),3));
sio_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(231:360,lats,levs,:),1),2),3));
sioV = squeeze(nansum(nansum(nansum(Vzyxsub_r_0(231:360,lats,levs),1),2),3));
sall_0 = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(:,lats,levs,:),1),2),3));
% Tendency of OHC (minus 1960) ZJ
tdpac = (spac_0(41:101)-spac_0(41))/10^21;
tdia = (sia_0(41:101)-sia_0(41))/10^21;
tdat = (sat_0(41:101)-sat_0(41))/10^21;
tdio = (sio_0(41:101)-sio_0(41))/10^21;
tdall = (sall_0(41:101)-sall_0(41))/10^21;
%%
close all
figure(1)
plot(tdia)
set(gca,'XLim',[1,61],'ylim',[-120,120],'ytick',[-120:30:120]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
ylabel('ZJ')
legend('advection','forcing','advection+forcing','location','north')
legend('boxoff')
saveas(gca,'/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240606/IO+Atl_tendency.png')

%% heat forcing (minus 1960) ZJ
% SHF
splon = 169; % 150 E
inlats = 46:84; % 55.1805 S :  34.8794 S
TAREA = 0.01^2*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','TAREA');
KMT = 100*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','KMT'); 
hflux_factor = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','hflux_factor');
rho_sw = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','rho_sw'); 
cp_sw = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','cp_sw'); 
ro = rho_sw*10^3 ; % sea water density kg/m3 (Antonov et al. 2004)
cp = cp_sw/10^4; % specific heat of water J/kg/K (Antonov et al. 2004)
SHF1 = ncread('/Volumes/CESM-post2/SHFensmean_1920-2005.nc','SHF'); % degC/s
SHF2 = ncread('/Volumes/CESM-post2/SHFensmean_2006-2080.nc','SHF'); % degC/s
SHFmon = cat(3,SHF1(:,:,480+1:end),SHF2(:,:,1:180));
clear SHF1 SHF2
[d1 d2 d3] = size(SHFmon);
SHF = squeeze(mean(reshape(SHFmon,d1,d2,12,d3/12),3)); % yearly
clear SHFmon
%%
SHFcli = nanmean(SHF,3);
map = SHFcli;
close all;
ftsz = 20; ticks = 100;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(map,[-ticks*5:1:ticks*5],ticks,TLONG',TLAT',ftsz)
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],'W m^-^2',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
saveas(gca,'/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240613/SHFcli.png')

%% trend
startyr = 1960;
endyr = 2020;
var = permute(SHF(:,:,startyr-1959:endyr-1959),[3 1 2]);
trd = trend_cal_3D(var);
% h0 = trendtest(var,0.05); % t test trend
% h0(find(isnan(Temp(:,:,1,1)) == 1)) = nan;
% trd(find(isnan(aa) == 1)) = nan;
max(trd,[],'all'), min(trd,[],'all')

units = 'W m^-^2'; name = 'NHFLX'; 
map = trd*10; ticks = 4;

close all;
ftsz = 20; 
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(map,[-ticks*5:ticks/10:ticks*5],ticks,TLONG',TLAT',ftsz)
m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
set_colorbar([0.83 0.08 0.03 0.88],[units,' decade^-^1'],5.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots0(h0,[.4 .4 .4],lonDataE,latDataE);
saveas(gca,'/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240613/SHFtrend_1960-2020.png')
%%
latstr = '35-55S'; % heat
var_r_0 = cat(1,SHF(splon:end,:,:),SHF(1:splon-1,:,:)); 
str1 = 'Surface Heat Flux (W m^-^2)';
str2 = 'SHF';
[spac_0] = areamean(var_r_0,1:125,inlats,TLAT); % 150E-70W
[sia_0] = areamean(var_r_0,126:320,inlats,TLAT); 
ln1 = spac_0;  max(ln1)
ln2 = sia_0; min(ln2)
ymax1 = 9; ymin1 = -12; yint1 = 3;
tt01 = 12; tt02 = 7;
tt11 = 27; tt12 = 7;
ymax2 = 18; ymin2 = -3; yint2 = 3;

close all;
sdiff = smooth(ln2-ln1);
ftsz = 14;
Fig = figure('position',[700 50 600 400]);

p1 = plot(smooth(ln1),'-','color',[.04 .35 .56],'LineWidth',2);
hold on
p2 = plot(smooth(ln2),'-','color',[.64 .08 .18],'LineWidth',2);
plot(zeros(1,61),'k','LineWidth',1)
yb = polyfit([1:61],ln1,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[.04 .35 .56],'linewidth',2)
text(tt01,tt02,[num2str(roundn(yb(1),-4))],'fontsize',ftsz,'Color',[.04 .35 .56])
yb = polyfit([1:61],ln2,1);
plot([1:61],polyval(yb,[1:61]),'--','Color',[.64 .08 .18],'linewidth',2)
text(tt11,tt12,[num2str(roundn(yb(1),-4))],'fontsize',ftsz,'Color',[.64 .08 .18])
set(gca,'XLim',[1,61]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',ftsz);
set(gca,'YLim',[ymin1,ymax1],'YTick',[ymin1:yint1:ymax1],'YColor','k');
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',ftsz),ylabel(str1,'FontSize',ftsz);

yyaxis right
p3 = bar(smooth(ln2-ln1));
set(p3,'FaceColor',[0.7,0.7,0.7],'facealpha',0.4,'EdgeColor',[0.7,0.7,0.7],'BarWidth',1)
set(gca,'YLim',[ymin2,ymax2],'YTick',[ymin2:yint2:ymax2],'YColor','k');
ylabel('Difference')

legend([p1,p2,p3],'Pacific','Atlantic-Indian Ocean','Difference','Location','north','Orientation','horizontal','fontsize',ftsz)
legend('boxoff')
saveas(gca,'/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240613/SHFtimeseries_35-55S_1960-2020.png')

%%
hfc = ro*cp*SHF*hflux_factor.*TAREA/10;
hfc_r = cat(1,hfc(splon:320,inlats,:),hfc(1:splon-1,inlats,:));
hfcp = squeeze(nansum(nansum(hfc_r(1:125,:,:),1),2));
hfcia = squeeze(nansum(nansum(hfc_r(126:320,:,:),1),2));
hfca = squeeze(nansum(nansum(hfc_r(126:205,:,:),1),2));
hfci = squeeze(nansum(nansum(hfc_r(206:320,:,:),1),2));
hfcall = squeeze(nansum(nansum(hfc_r,1),2));

clear foc*
focp = cumsum(hfcp)*60*60*24*365/10^21;
focia = cumsum(hfcia)*60*60*24*365/10^21;
foca = cumsum(hfca)*60*60*24*365/10^21;
foci = cumsum(hfci)*60*60*24*365/10^21;
focall = cumsum(hfcall)*60*60*24*365/10^21;

% clear foc*
% for s = 1:length(hfcp)
%     focp(s,1) = sum(hfcp)/61*s*60*60*24*365/10^21;
%     focia(s,1) = sum(hfcia)/61*s*60*60*24*365/10^21;
%     foca(s,1) = sum(hfca)/61*s*60*60*24*365/10^21;
%     foci(s,1) = sum(hfci)/61*s*60*60*24*365/10^21;
%     focall(s,1) = sum(hfcall)/61*s*60*60*24*365/10^21;
% 
% end
%
close all
Fig = figure('position',[700 100 800 400]);
plot(focp,'r')
hold on
plot(focia,'b')
plot(focall,'k')
%%
set(gca,'XLim',[1,61],'ylim',[-100,100],'ytick',[-100:50:100]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
ylabel('ZJ')
legend('tendency','forcing','location','north')
legend('boxoff')
% print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_Warming/Tendency_IO_',depthstr,'_35S_55S.png'],'-dpng','-r300')

%%
clear Tem* Tzy* VV* Vzy* flnsMMM lhflxMMM shflxMMM flnsMMM
intlev = 37; % 707 m
% ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)

% cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
DXU = 0.01*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','DXU'); % cm -> m
DXT = 0.01*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','DXT');
DYU = 0.01*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','DYU');
DYT = 0.01*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','DYT');
DZU = 0.01*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','dz');
DZT = 0.01*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','dz');
TLONG = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','TLONG');
TLAT = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','TLAT');
TAREA = 0.01^2*ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','TAREA');
ULONG = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','ULONG');
ULAT = ncread('/Volumes/CESM-post2/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc','ULAT');
%% adv UT
UET1 = ncread('/Volumes/CESM-post2/UETensmean_1920-2005.nc','UET'); % degC/s
UET2 = ncread('/Volumes/CESM-post2/UETensmean_2006-2080.nc','UET'); % degC/s
UETmon = cat(4,UET1(:,:,:,480+1:end),UET2(:,:,:,1:180));
clear UET1 UET2
[d1 d2 d3 d4] = size(UETmon);
UET = squeeze(mean(reshape(UETmon,d1,d2,d3,12,d4/12),4)); % yearly
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
uto = ro*cp* dutdx .* VOL;
clear dutdx
utor = cat(1,uto(splon:end,inlats,1:intlev,:),uto(1:splon-1,inlats,1:intlev,:)); % 1960-2005
clear uto
ULONGr = cat(1,ULONG(splon:end,:),ULONG(1:splon-1,:));
utop = squeeze(nansum(nansum(nansum(utor(1:125,:,:,:),1),2),3)); % 150E-70W; 289.625 ~ 70W
utoia = squeeze(nansum(nansum(nansum(utor(126:end,:,:,:),1),2),3)); 
utoa = squeeze(nansum(nansum(nansum(utor(126:205,:,:,:),1),2),3)); 
utoi = squeeze(nansum(nansum(nansum(utor(206:end,:,:,:),1),2),3)); 
utoall = squeeze(nansum(nansum(nansum(utor,1),2),3)); 
%
close all
figure(1)
plot(utoall,'k')
hold on
plot(utop,'r')
plot(utoia,'b')
set(gca,'ylim',[-3,3]*10^15);

%adv VT
VNT1 = ncread('/Volumes/CESM-post2/VNTensmean_1920-2005.nc','VNT'); % degC/s
VNT2 = ncread('/Volumes/CESM-post2/VNTensmean_2006-2080.nc','VNT'); % degC/s
VNTmon = cat(4,VNT1(:,:,:,480+1:end),VNT2(:,:,:,1:180));
clear VNT1 VNT2
[d1 d2 d3 d4] = size(VNTmon);
VNT = squeeze(mean(reshape(VNTmon,d1,d2,d3,12,d4/12),4)); % yearly
clear VNTmon
[d1 d2 d3 d4] = size(VNT);
VOLey = cat(2,VOL,VOL(:,end,:));  
VNTe = cat(2,VNT,VNT(:,end,:,:)); % extend to calculate difference
clear VNT
dvtdy = -diff(VNTe.*VOLey,1,2)./VOL;
clear VNTe
vto = ro*cp* dvtdy .* VOL;
clear dvtdy
vtor = cat(1,vto(splon:end,inlats,1:intlev,:),vto(1:splon-1,inlats,1:intlev,:)); % 1960-2005
clear vto
vtop = squeeze(nansum(nansum(nansum(vtor(1:125,:,:,:),1),2),3)); % 150E-70W; 289.625 ~ 70W
vtoia = squeeze(nansum(nansum(nansum(vtor(126:end,:,:,:),1),2),3)); 
vtoa = squeeze(nansum(nansum(nansum(vtor(126:205,:,:,:),1),2),3)); 
vtoi = squeeze(nansum(nansum(nansum(vtor(206:end,:,:,:),1),2),3)); 
vtoall = squeeze(nansum(nansum(nansum(vtor,1),2),3)); 
%
figure(2)
plot(vtoall,'k')
hold on
plot(vtop,'r')
plot(vtoia,'b')
set(gca,'ylim',[-3,3]*10^15);
%%
figure(3)
plot(utoall+vtoall,'k')
hold on
plot(utop+vtop,'r')
plot(utoia+vtoia,'b')
%%
figure(4)
plot(cumsum(utoall+vtoall),'k')
hold on
plot(cumsum(utop+vtop),'r')
plot(cumsum(utoia+vtoia),'b')
%% adv WT
WTT1 = ncread('/Volumes/CESM-post2/WTTensmean_1920-2005.nc','WTT'); % degC/s
WTT2 = ncread('/Volumes/CESM-post2/WTTensmean_2006-2080.nc','WTT'); % degC/s
WTTmon = cat(4,WTT1(:,:,:,480+1:end),WTT2(:,:,:,1:180));
clear WTT1 WTT2
[d1 d2 d3 d4] = size(WTTmon);
WTT = squeeze(mean(reshape(WTTmon,d1,d2,d3,12,d4/12),4)); % yearly
clear WTTmon
[d1 d2 d3 d4] = size(WTT);
VOLez = cat(3,VOL,VOL(:,:,end));
WTTe = cat(3,WTT,WTT(:,:,end,:)); % extend to calculate difference
clear WTT
dztdz = diff(WTTe.*VOLez,1,3)./VOL;
clear WTTe
wto = ro*cp* dztdz .* VOL;
clear dztdz
wtor = cat(1,wto(splon:end,inlats,1:intlev,:),wto(1:splon-1,inlats,1:intlev,:)); % 1960-2005
clear wto
wtop = squeeze(nansum(nansum(nansum(wtor(1:125,:,:,:),1),2),3)); % 150E-70W; 289.625 ~ 70W
wtoia = squeeze(nansum(nansum(nansum(wtor(126:end,:,:,:),1),2),3)); 
wtoa = squeeze(nansum(nansum(nansum(wtor(126:205,:,:,:),1),2),3)); 
wtoi = squeeze(nansum(nansum(nansum(wtor(206:end,:,:,:),1),2),3)); 
wtoall = squeeze(nansum(nansum(nansum(wtor,1),2),3)); 
%
figure(5)
plot(wtoall,'k')
hold on
plot(wtop,'r')
plot(wtoia,'b')
%%
advall = 60*60*24*365*cumsum(utoall+vtoall+wtoall)/10^21;
advp = 60*60*24*365*cumsum(utop+vtop+wtop)/10^21;
advia = 60*60*24*365*cumsum(utoia+vtoia+wtoia)/10^21;
adva = 60*60*24*365*cumsum(utoa+vtoa+wtoa)/10^21;
advi = 60*60*24*365*cumsum(utoi+vtoi+wtoi)/10^21;
%%
close all
figure(6)
plot(advall,'b')
hold on
plot(focall,'r')
plot(focall+advall,'k')
plot([1:61],zeros(61),'k--')
set(gca,'XLim',[1,61],'ylim',[-1000,1000],'ytick',[-1000:500:1000]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
ylabel('ZJ')
legend('advection','forcing','advection+forcing','location','north')
legend('boxoff')
saveas(gca,'/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240606/ALL_Adv+forc1.png')
%% 
close all
figure(1)
plot(cumsum(diff(spac_0(40:101)))*10,'b')
hold on
plot(cumsum((wtop+vtop+utop+hfcp)*60*60*24*365))
%%
close all
figure(1)
plot(diff(sia_0(40:101))*10,'k')
hold on
plot((wtoia+vtoia+utoia+hfcia)*60*60*24*365,'b')
plot(diff(sia_0(40:101))*10-(wtoia+vtoia+utoia+hfcia)*60*60*24*365,'r--')
figure(2)
plot(wtoia+vtoia+utoia,'b')
hold on
plot(hfcia,'r')
%% vertical mixing
clear utor vtor wtor
DIA_IMPVF_TEMP1 = ncread('/Volumes/CESM-post2/DIA_IMPVF_TEMPensmean_1920-2005.nc','DIA_IMPVF_TEMP'); % degC/s
DIA_IMPVF_TEMP2 = ncread('/Volumes/CESM-post2/DIA_IMPVF_TEMPensmean_2006-2080.nc','DIA_IMPVF_TEMP'); % degC/s
DIA_IMPVF_TEMPmon = cat(4,DIA_IMPVF_TEMP1(:,:,:,480+1:end),DIA_IMPVF_TEMP2(:,:,:,1:180));
clear DIA_IMPVF_TEMP1 DIA_IMPVF_TEMP2
[d1 d2 d3 d4] = size(DIA_IMPVF_TEMPmon);
DIA_IMPVF_TEMP = squeeze(mean(reshape(DIA_IMPVF_TEMPmon,d1,d2,d3,12,d4/12),4)); % yearly
clear DIA_IMPVF_TEMPmon
DIA_IMPVF_TEMPe = cat(3,DIA_IMPVF_TEMP,DIA_IMPVF_TEMP(:,:,end,:));
clear DIA_IMPVF_TEMP
DIA = -diff(DIA_IMPVF_TEMPe.*TAREA,1,3)./VOL;
clear DIA_IMPVF_TEMPe 
%%
KPP_SRC_TEMP1 = ncread('/Volumes/CESM-post2/KPP_SRC_TEMPensmean_1920-2005.nc','KPP_SRC_TEMP'); % degC/s
KPP_SRC_TEMP2 = ncread('/Volumes/CESM-post2/KPP_SRC_TEMPensmean_2006-2080.nc','KPP_SRC_TEMP'); % degC/s
KPP_SRC_TEMPmon = cat(4,KPP_SRC_TEMP1(:,:,:,480+1:end),KPP_SRC_TEMP2(:,:,:,1:180));
clear KPP_SRC_TEMP1 KPP_SRC_TEMP2
[d1 d2 d3 d4] = size(KPP_SRC_TEMPmon);
KPP = squeeze(mean(reshape(KPP_SRC_TEMPmon,d1,d2,d3,12,d4/12),4)); % yearly

vmo = ro*cp*(DIA+KPP) .* VOL;
clear dztdz
vmor = cat(1,vmo(splon:end,inlats,1:intlev,:),vmo(1:splon-1,inlats,1:intlev,:)); % 1960-2005
clear vmo
vmop = squeeze(nansum(nansum(nansum(vmor(1:125,:,:,:),1),2),3)); % 150E-70W; 289.625 ~ 70W
vmoia = squeeze(nansum(nansum(nansum(vmor(126:end,:,:,:),1),2),3)); 
vmoa = squeeze(nansum(nansum(nansum(vmor(126:205,:,:,:),1),2),3)); 
vmoi = squeeze(nansum(nansum(nansum(vmor(206:end,:,:,:),1),2),3)); 
vmoall = squeeze(nansum(nansum(nansum(vmor,1),2),3)); 
%%
vmixall = 60*60*24*365*cumsum(vmoall)/10^21;
vmixp = 60*60*24*365*cumsum(vmop)/10^21;
vmixia = 60*60*24*365*cumsum(vmoia)/10^21;
vmixa = 60*60*24*365*cumsum(vmoa)/10^21;
vmixi = 60*60*24*365*cumsum(vmoi)/10^21;

figure(6)
plot(vmixall,'k')
hold on
plot(vmixp,'r')
plot(vmixa,'b')

%% horizontal diffusion
HDIFE_TEMP1 = ncread('/Volumes/CESM-post2/HDIFE_TEMPensmean_1920-2005.nc','HDIFE_TEMP'); % degC/s
HDIFE_TEMP2 = ncread('/Volumes/CESM-post2/HDIFE_TEMPensmean_2006-2080.nc','HDIFE_TEMP'); % degC/s
HDIFE_TEMPmon = cat(4,HDIFE_TEMP1(:,:,:,480+1:end),HDIFE_TEMP2(:,:,:,1:180));
clear HDIFE_TEMP1 HDIFE_TEMP2
[d1 d2 d3 d4] = size(HDIFE_TEMPmon);
HDIFE_TEMP = squeeze(mean(reshape(HDIFE_TEMPmon,d1,d2,d3,12,d4/12),4)); % yearly
clear HDIFE_TEMPmon
HDIFE_TEMPe = cat(1,HDIFE_TEMP,HDIFE_TEMP(1,:,:,:));
clear HDIFE_TEMP
dhdfe = diff(HDIFE_TEMPe.*VOL,1,1)./VOL;
clear HDIFE_TEMPe
%%
HDIFN_TEMP1 = ncread('/Volumes/CESM-post2/HDIFN_TEMPensmean_1920-2005.nc','HDIFN_TEMP'); % degC/s
HDIFN_TEMP2 = ncread('/Volumes/CESM-post2/HDIFN_TEMPensmean_2006-2080.nc','HDIFN_TEMP'); % degC/s
HDIFN_TEMP = cat(4,HDIFN_TEMP1(:,:,:,480+1:end),HDIFN_TEMP2(:,:,:,1:180));
clear HDIFN_TEMP1 HDIFN_TEMP2
ynan = nan(d1,1,d3,d4);
dhdfn = cat(2,ynan,diff(HDIFN_TEMP.*VOL,1,2))./VOL;
hmix = dhdfe + dhdfn;


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
function [ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(sum(sum((cos(latData(lons,lats)/180*pi)).*var1(:,:,:),1,'omitmissing'),2,'omitmissing')/sum(cos(latData(lons,lats)/180*pi).*var2,'all','omitmissing'),size(var1,3),1);
end