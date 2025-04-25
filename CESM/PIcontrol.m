clc,clear,close all;
addpath D:\1_matlab\help;
addpath D:\1_matlab\help\seawater\;
datadir='E:\CESM-post\JC-data\pi_control\TEMP\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
k=length(filelist);
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
dweit = depthData(2:37)-depthData(1:36);
lats = 1:90;
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    TEMP = ncread(filename1,'TEMP'); %读入变量
    TEMPsub(:,:,:,:,s) = TEMP(:,1:55,:,1:100);
    netcdf.close(ncid);  %关闭文件    
    Tsub0(:,:,:,s) = permute(nansum(cat(3,TEMP(:,lats,1,1:100)*5,TEMP(:,lats,2:37,1:100).*(permute(dweit,[3 2 1]))),3)/700,[1 2 4 3]); % 0-700m
    sst0(:,:,:,s) = permute(TEMP(:,:,1,1:100),[1 2 4 3]);
    % 46-61S
    Tem1m0(:,:,:,s) = latmean(TEMP(:,:,1:37,1:100),[29:44],latData);

end
TEMPsubr = reshape(TEMPsub,[360,55,37,100*18]);
clear TEMPsub
Tsub1 = reshape(Tsub0,[360,90,100*18]);
sst1 = reshape(sst0,[360,180,100*18]);
Tem1m1 = reshape(Tem1m0,[360,37,100*18]);
dT = 1; cf = 1/8;
for i = 1:size(Tsub1,1);
    for j = 1:size(Tsub1,2);
        Tsub2(:,i,j) = lanczosfilter(permute(Tsub1(i,j,:),[3 1 2]),dT,cf,[],'low'); % 8 year filtered 
    end
end
Tsub = permute(Tsub2,[2 3 1]);
for i = 1:size(Tem1m1,1);
    for j = 1:size(Tem1m1,2);
        Tem1m2(:,i,j) = lanczosfilter(permute(Tem1m1(i,j,:),[3 1 2]),dT,cf,[],'low'); % 8 year filtered 
    end
end
Tem1m = permute(Tem1m2,[2 3 1]);
%% VVEL
datadir='E:\CESM-post\pi_control\VVEL\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    VVEL = ncread(filename1,'VVEL'); %读入变量
    VVELsub(:,:,:,:,s) = VVEL(:,1:55,:,1:100);
    netcdf.close(ncid);  %关闭文件    
end
VVELsubr = reshape(VVELsub,[360,55,37,100*18]);
clear VVELsub
%% pi-control的TAUX TAUY与LE一样，符号都完全与观测相反
datadir='E:\CESM-post\pi_control\TAUX\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    TAUX = ncread(filename1,'TAUX'); %读入变量
    netcdf.close(ncid);  %关闭文件    
    taux0(:,:,:,s) = TAUX(:,:,1:100);
end
taux1 = reshape(taux0,[360,180,100*18]);
tauxa = taux1 - mean(taux1,3); % anomaly
datadir='E:\CESM-post\pi_control\TAUY\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    TAUY = ncread(filename1,'TAUY'); %读入变量
    netcdf.close(ncid);  %关闭文件    
    tauy0(:,:,:,s) = TAUY(:,:,1:100);
end
tauy1 = reshape(tauy0,[360,180,100*18]);
tauya = tauy1 - mean(tauy1,3); % anomaly
dT = 1; cf = 1/8;
for i = 1:size(tauxa,1);
    for j = 1:size(tauxa,2);
        tauxaf(:,i,j) = lanczosfilter(permute(tauxa(i,j,:),[3 1 2]),dT,cf,[],'low'); % 8 year filtered 
        tauyaf(:,i,j) = lanczosfilter(permute(tauya(i,j,:),[3 1 2]),dT,cf,[],'low'); % 8 year filtered 
    end
end
taux = permute(tauxaf,[2 3 1]);
tauy = permute(tauyaf,[2 3 1]);

%% curlz
parfor i = 1:size(taux1,3);
    % curlz lat*lon *time
    curlz0(:,:,i) = ra_windstrcurl(latData,lonData,taux1(:,:,i)',tauy1(:,:,i)',1);
end
curlz1 = permute(curlz0,[2 1 3]); % lon*lat*time
curlza = curlz1 - mean(curlz1,3); % anomaly
dT = 1; cf = 1/8;
for i = 1:size(curlza,1);
    for j = 1:size(curlza,2);
        curlzaf(:,i,j) = lanczosfilter(permute(curlza(i,j,:),[3 1 2]),dT,cf,[],'low'); % 8 year filtered 
    end
end
curlz = permute(curlzaf,[2 3 1]);
%% mixed layer depth
datadir='E:\CESM-post\pi-control\HMXL\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
k=length(filelist);
dweit = depthData(2:37)-depthData(1:36);
lats = 1:90;
for s=1:k
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    HMXL = ncread(filename1,'HMXL'); %读入变量
    HMXLsub(:,:,:,s) = HMXL(:,:,1:100);
    netcdf.close(ncid);  %关闭文件    
end
HMXLsubr = reshape(HMXLsub,[360,180,100*18]);
clear HMXLsub
dT = 1; cf = 1/8;
for i = 1:size(HMXLsubr,1);
    for j = 1:size(HMXLsubr,2);
        HMXLf0(:,i,j) = lanczosfilter(permute(HMXLsubr(i,j,:),[3 1 2]),dT,cf,[],'low'); % 8 year filtered 
    end
end
HMXLf = permute(HMXLf0,[2 3 1]);

%% TPI index
ssta = sst1 - mean(sst1,3);
[ts1_zs ts1] = areamean(ssta,141:216,116:136,latData); % 25N-45N,140E-145W
[ts2_zs ts2] = areamean(ssta,171:271,80:101,latData); % 10S-10N,170E-90W
[ts3_zs ts3] = areamean(ssta,151:201,40:75,latData); % 50S-15S,150E-160W
TPI1 = ts2-(ts1+ts3)/2; % unfiltered IPO index///
dT = 1; cf = 1/8;
TPIf = lanczosfilter(TPI1,dT,cf,[],'low'); % 8 year filtered IPO index
TPIfz = zscore(TPIf);
close all
plot(TPIfz)
%% vertical mean TEMP
[d1 d2 d3] = size(Tsub);
Tsubmr1 = reshape(permute(Tsub,[3 1 2]),[d3,d1*d2]);
index = TPIfz; pngname = 'IPO';
clear par11 h01
parfor i = 1:size(Tsubmr1,2);
    [par11(i),h01(i),t] = reg1_ttest(index,Tsubmr1(:,i),0.1,0);
end
par11r = reshape(par11,[d1,d2]);
h01r = reshape(h01,[d1,d2]);
%
close all;
map = par11r;
ticks = 0.15;  ftsz = 12;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(lats),12)
testdots(h01r,[.4 .4 .4],lonData,latData(lats));
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
print(Fig,['D:\figures\CESM\Yearly\PI-control\20230606_IPO_SouthernOcean\',pngname,'.reg.Tem.verticalmean_0_700m.png'],'-dpng','-r300')
%% meridional mean TEMP

[d1 d2 d3] = size(Tem1m);
Tem1mr = reshape(permute(Tem1m,[3 1 2]),[d3,d1*d2]);
%
index = TPIfz; pngname = 'IPO';
clear par14 h04
parfor i = 1:size(Tem1mr,2);
    [par14(i),h04(i),t] = reg1_ttest(index,Tem1mr(:,i),0.1,0);
end
par14r = reshape(par14,[d1 d2]);
max(par14r,[],'all')
h04r = reshape(h04,[d1 d2]);
%% mixed layer depth
lats = [29:44]; names = '46S-61S';
latm = latmean(HMXLf,lats,latData);
HMXLr = reshape(permute(latm,[2 1]),[1800,360]);
parfor i = 1:size(HMXLr,2);
    [paH(i),hH(i),t] = reg1_ttest(TPIfz(1:1800),HMXLr(:,i),0.1,1);
end
%
latmMMM = latmean(HMXLsubr,lats,latData);
HMXLcli = nanmean(latmMMM,2);
HMXLclir = -cat(1,HMXLcli(160:360,:),HMXLcli(1:159,:))/100; % cm -> m

%%
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
map = nanmean(par14r,3);
map_r = cat(1,map(160:300,:),map(301:360,:),map(1:159,:));
h04rr = cat(1,h04r(160:300,:),h04r(301:360,:),h04r(1:159,:));
h04rr(1:2:end) = 0;
%
close all;
Fig = figure('position',[100 100 800 400]);
ax = axes('Position',[0.08 0.1 0.72 0.87],'fontsize',14,'box','on');
caxval = 0.1;
[c h] = contourf(lonData_r,-depthData(1:37),map_r',[-1:caxval/10:1],'linestyle','none');
caxis([-caxval,caxval])
hold on
line([300 300],[0 -700],'linewidth',1.5,'color','k');
plot(lonData_r,HMXLclir,'linewidth',1.5,'color','k','LineStyle','--')
[xx,yy] = find(h04rr==1);
plot(lonData_r(xx),-depthData(yy),'.','color',[.5 .5 .5],'markersize',3);
load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-caxval:caxval/5:caxval],'position',[0.9 0.1 0.02 0.87]);
set(ch.Label,'String','K','Units','normalized','position',[4.5 0.5 0],'Rotation',-90,'fontsize',12);
set(gca,'XLim',[160,360+160]);
set(gca,'XTick',[160:60:360+160]);
set(gca,'XTickLabel',{'160^oE','140^oW','80^oW','20^oW','40^oE','100^oE','160^oE'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')

yyaxis right
maph = paH'/100; % cm->m
maphr = cat(1,maph(160:360,:),maph(1:159,:));
plot(lonData_r,-maphr,'linewidth',2.5,'color',[.5 .5 .5],'LineStyle','--')
set(gca,'YLim',[-5.5,2.2],'YTick',[-5.5:1.1:2.2],'YColor',[.5 .5 .5]);
yh = ylabel('Mixed Layer Depth (m)')
set(yh,'Rotation',270,'position',[562,-1.73,0])
print(Fig,['D:\figures\CESM\Yearly\PI-control\20240130_IPO_SO\',pngname,'.reg.Tem.meridionalmean46-61S.png'],'-dpng','-r300')
%% IPO .reg. taux tauy Curlz
tauxr = permute(reshape(taux(:,1:55,:),360*55,1800),[2 1]);
tauyr = permute(reshape(tauy(:,1:55,:),360*55,1800),[2 1]);
curlzr = permute(reshape(curlz(:,1:55,:),360*55,1800),[2 1]);
clear par14 h04
parfor i = 1:size(curlzr,2);
    [par12(i),h02(i),t] = reg1_ttest(TPIfz,tauxr(:,i),0.1,1);
    [par13(i),h03(i),t] = reg1_ttest(TPIfz,tauyr(:,i),0.1,1);
    [par14(i),h04(i),t] = reg1_ttest(TPIfz,curlzr(:,i),0.1,1);
end
par12r = -reshape(par12,[360,55]);
par13r = -reshape(par13,[360,55]);
par14r = -reshape(par14,[360,55]);
h04r = reshape(h04,[360,55]);
h04r(:,1:15) = 0;
%%
close all; ftsz = 12;
map = mean(par14r,3)*10^9;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 8; 
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(1:55),12)
set_colorbar([0.83 0.08 0.03 0.88],'Wind Stress Curl  (10^-^9 Pa m^-^1)',4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks]);
hold on
testdots(h04r,[.4 .4 .4],lonData,latData(1:55));
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
d = 5; dd = 0.0006;
[x,y] = meshgrid(lonData,latData(1:55));
umap = par12r.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par13r(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1,'MaxHeadSize',5) 
print(Fig,['D:\figures\CESM\Yearly\PI-control\20230606_IPO_SouthernOcean\IPO.reg.curlz&Tau.png'],'-dpng','-r300')
%% MOHT-Ekman
[MOHTall MOHT1 MOHT2] = MOHT_Ekman(taux1(:,1:55,:),TEMPsubr,latData(1:55),depthData);
%
close all;
Fig = figure('position',[100 100 500 500])
map = mean(MOHTall,2)/10^15;
% map(85:96,:) = nan;
plot(latData(1:55),map,'linewidth',1.5)
set(gca,'XLim',[-90,90]);
set(gca,'XTick',[-80:20:80]);
ylabel('MOHT (PW)');xlabel('Latitude')
% print(Fig,['D:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\MOHT_climotology.png'],'-dpng','-r300')
%% MOHT anomaly,detrend, and 8-yr filter
MOHT1a = (MOHT1-mean(MOHT1,2))'; % anomaly
MOHT2a = (MOHT2-mean(MOHT2,2))'; % anomaly
dT = 1; cf = 1/8;
clear MOHT1af MOHT2af
for i = 1:55
    MOHT1af(:,i) = lanczosfilter(MOHT1a(:,i),dT,cf,[],'low'); % 8 year filtered
    MOHT2af(:,i) = lanczosfilter(MOHT2a(:,i),dT,cf,[],'low'); % 8 year filtered
end
% MOHT1af(1:4,:) = []; MOHT2af(1:4,:) = [];
% MOHT1af(end-3:end,:) = []; MOHT2af(end-3:end,:) = [];
%
index = TPIfz; pngname = 'IPO';
clear par7e h07e par8e h08e 
for i = 1:size(MOHT1af,2);
    [par7e(i),h07e(i),t] = reg1_ttest(index,MOHT1af(:,i),0.1,1);
    [par8e(i),h08e(i),t] = reg1_ttest(index,MOHT2af(:,i),0.1,1);
end
%%
map7e = -par7e;
map8e = -par8e;
close all;
Fig = figure('position',[100 100 800 400]);
plot(latData(1:55),map7e/10^15,'r','linewidth',1.5)
hold on
plot(latData(1:55),map8e/10^15,'b','linewidth',1.5)
plot(latData(1:55),zeros(1,length(latData(1:55))),'k','linewidth',1.5)
plot(-61*ones(20),[-.02:.04/19:.02],'k--','linewidth',1.5);
plot(-46*ones(20),[-.02:.04/19:.02],'k--','linewidth',1.5);
scatter(latData(find(h07e == 1)),map7e(find(h07e == 1))/10^15,'*','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h08e == 1)),map8e(find(h08e == 1))/10^15,'*','MarkerEdgeColor','b','linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
set(gca,'YLim',[-.02,.02],'YTick',[-.02:.01:.02]);
set(gca,'YTick',yticks,'FontSize',14);
ylabel('MOHT (PW)')
legend('Pac','Atl&IO','Location','north','Orientation','horizontal','position',[0.67,0.85,0.22,0.06])
legend('boxoff')
% print(Fig,['D:\figures\CESM\Yearly\PI-control\20230606_IPO_SouthernOcean\IPO.reg_zonalmean.MOHT_Ekman.png'],'-dpng','-r300')
%% MOHT vT 0-700m
%
f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:55)*pi/180)/360; % distance between 2 longitudes 
levdist(1) = 4; levdist(2:37) = dweit; % z distance
Temp = TEMPsubr;
Vvel = VVELsubr/100; % cm/s -> m/s
clearvars -except f ro cp londist levdist Temp Vvel map7e map8e ...
    h07e h08e TPIfz MOHT1af MOHT2af lonData latData
% 0-700m MOHT
TemVel = Temp.*Vvel;
OHCz = permute(ro*cp*sum((TemVel).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z
OHCzy = OHCz.*londist';
OHCzy_r = cat(1,OHCzy(160:300,:,:,:),OHCzy(301:360,:,:,:),OHCzy(1:159,:,:,:));
MOHT71_p = permute(nansum(OHCzy_r(1:141,:,:,:),1),[2 3 4 1]); % Pac unit:W
MOHT72_p = permute(nansum(OHCzy_r(142:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
MOHT7all_p = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
clear OHC* TemVel Vvel filename* filelist* datadir*
%% MOHT anomaly,detrend, and 8-yr filter
MOHT71ad = (MOHT71_p-mean(MOHT71_p,2))'; % anomaly
MOHT72ad = (MOHT72_p-mean(MOHT72_p,2))'; % anomaly
dT = 1; cf = 1/8;
clear MOHT71adf MOHT72adf
for i = 1:55
    MOHT71adf(:,i) = lanczosfilter(MOHT71ad(:,i),dT,cf,[],'low'); % 8 year filtered
    MOHT72adf(:,i) = lanczosfilter(MOHT72ad(:,i),dT,cf,[],'low'); % 8 year filtered
end
% MOHT71adf(1:4,:) = []; MOHT72adf(1:4,:) = [];
% MOHT71adf(end-3:end,:) = []; MOHT72adf(end-3:end,:) = [];
%
index = TPIfz; pngname = 'IPO';
clear par7 h07 par8 h08 
for i = 1:size(MOHT71adf,2);
    [par7(i),h07(i),t] = reg1_ttest(index,MOHT71adf(:,i),0.1,1);
    [par8(i),h08(i),t] = reg1_ttest(index,MOHT72adf(:,i),0.1,1);
end
%%
close all;
Fig = figure('position',[100 100 800 400]);
p1 = plot(latData(1:55),par7/10^15,'-','color',[0.90,0.51,0.31],'linewidth',1.5)
hold on
p2 = plot(latData(1:55),par8/10^15,'-','color',[0.17,0.73,0.69],'linewidth',1.5)
p3 = plot(latData(1:55),map7e/10^15,'r-','linewidth',1.5)
p4 = plot(latData(1:55),map8e/10^15,'b-','linewidth',1.5)
plot(latData(1:55),zeros(1,length(latData(1:55))),'k','linewidth',1.5)
plot(-61*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
plot(-46*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
scatter(latData(find(h07 == 1)),par7(find(h07 == 1))/10^15,'o','MarkerEdgeColor',[0.90,0.51,0.31],'linewidth',1)
scatter(latData(find(h08 == 1)),par8(find(h08 == 1))/10^15,'o','MarkerEdgeColor',[0.17,0.73,0.69],'linewidth',1)
scatter(latData(find(h07e == 1)),map7e(find(h07e == 1))/10^15,'o','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h08e == 1)),map8e(find(h08e == 1))/10^15,'o','MarkerEdgeColor','b','linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
set(gca,'YLim',[-.025,.032],'YTick',[-.03:.01:.03]);
set(gca,'YTick',yticks,'FontSize',14);
ylabel('MOHT (PW)')
% legend([p1,p2,p3,p4],'Pac_7_0_0_m','Atl&IO_7_0_0_m','Pac_E_k_m_a_n','Atl&IO_E_k_m_a_n','Location','north','Orientation','vertical','position',[0.14,0.6,0.18,0.3])
% legend('boxoff')
print(Fig,['D:\figures\CESM\Yearly\PI-control\20230606_IPO_SouthernOcean\IPO.reg_zonalmean.MOHT_700m_&_Ekman.png'],'-dpng','-r300')















function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end
function contourfSPolar(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end
function [ts] = latmean(var,lats,latData)
% average along longitude (all latitudes)
% var is lon*lat*depth*time
% ts is lon*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,2)./nansum(cos(latData(lats)'/180*pi).*var2,2),[1 3 4 2]);
end
function testdots(h,clor,lonData,latData)
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',2,'color',clor); % F-test dots
end
function [MOHTall MOHT1 MOHT2] = MOHT_Ekman(tauxraw,TemLEMsub,latData,depthData)
    f = sw_f(latData);
    My = -tauxraw./f'; % Ekman 质量输运
    ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
    cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
    Vy = My/ro;
%     Vy(:,85:96,:) = nan; % 赤道不适用ekman输送
    londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
    dweit = depthData(2:37) - depthData(1:36);
    levdist(1) = 1; levdist(2:37) = dweit; % z distance
    clear OHC OHCT MOHT
    OHC = permute(ro*cp*sum((TemLEMsub(:,1:55,1:37,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3])/700;
    OHCT = OHC.*Vy;
    OHCTw = OHCT.*londist';
    MOHTall = permute(nansum(OHCTw,1),[2 3 1]); % zonal sum unit:W
    OHCTw_r = cat(1,OHCTw(160:300,:,:),OHCTw(301:360,:,:),OHCTw(1:159,:,:));
    MOHT1 = permute(nansum(OHCTw_r(1:141,:,:),1),[2 3 1]); % Pac unit:W
    MOHT2 = permute(nansum(OHCTw_r(142:end,:,:),1),[2 3 1]); % IO+Atl unit:W
end