% Data 
% clc,clear,close all;
addpath D:\1_matlab\help;
addpath D:\1_matlab\help\seawater\;
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = 37;
dweit = depthData(2:37)-depthData(1:37-1);
ncid=netcdf.open('E:\CESM-post\JC-data\LE\Temperature\his_sst_ensmean.nc','NOWRITE');
ncdisp('E:\CESM-post\JC-data\LE\Temperature\his_sst_ensmean.nc');
tsMMM = ncread('E:\CESM-post\JC-data\LE\Temperature\his_sst_ensmean.nc','TEMP'); 
sstMMM = permute(tsMMM(:,:,1,:),[1 2 4 3]);
datadir='E:\CESM-post\JC-data\LE\Temperature\his_sst_yr1x1\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); 
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
%
for s = 1:k
    s
    filename=[datadir,filelist(s).name];    
    ncid=netcdf.open(filename,'NC_NOWRITE');
    Tem = ncread(filename,'TEMP');
% sst minus MMM & anomaly (without filter)
    sst0 = permute(Tem(:,:,1,:),[1,2,4,3]);
    sst1 = sst0 - sstMMM; % minus MMM (trend)
    sstd(:,:,:,s) = sst1 - mean(sst1,3); % anomaly

%     load(['MatFile/0_700m_South/Temsubadf',num2str(s),'.mat']);
%     Temsubadf(:,:,:,1:4) = []; Temsubadf(:,:,:,end-3:end) = [];
% % sst minus MMM & anomaly & (8-yr filter) 
%     load(['MatFile/sst_Global/sstadf',num2str(s),'.mat']);
%     sstadf(:,:,1:4) = []; sstadf(:,:,end-3:end) = [];
%     sstdf(:,:,:,s) = sstadf; 
% % 0-707m mean
%     Tsubm(:,:,:,s) = permute(nansum(cat(3,Temsubadf(:,:,1,:)*5,Temsubadf(:,:,2:37,:).*permute(dweit,[3 2 1 4])),3)/707,[1 2 4 3]);
% meridional mean
%     names = '46S-61S';
%     Tem1m(:,:,:,s) = latmean(Temsubadf(:,:,1:37,:),[29:44],latData);
% end
% %
% for s = 1:k
%     load(['MatFile/nhflx_Global/nhflxadf',num2str(s),'.mat']);
%     nhflxadf(:,:,1:4) = []; nhflxadf(:,:,end-3:end) = [];
%     nhflxdf(:,:,:,s) = nhflxadf;
% end
% for s = 1:40
%     load(['MatFile/slp_Global/slpadf',num2str(s),'.mat']);
%     slpadf(:,:,1:4) = []; slpadf(:,:,end-3:end) = [];
%     slp(:,:,:,s) = double(slpadf);
end
% nyr = size(Tsubm,3);
clear Tem sst0 sst1 Temsubadf sstadf
%% TPI and AMO index
nmem = 40; % number of member
clear TPIfz AMOf8z
for s=1:nmem
    s
% TPI index ---------------------------------------------------------------
    [ts1_zs ts1] = areamean(sstd(:,:,:,s),141:216,116:136,latData); % 25N-45N,140E-145W
    [ts2_zs ts2] = areamean(sstd(:,:,:,s),171:271,80:101,latData); % 10S-10N,170E-90W
    [ts3_zs ts3] = areamean(sstd(:,:,:,s),151:201,40:75,latData); % 50S-15S,150E-160W
    TPI1 = ts2-(ts1+ts3)/2; % unfiltered IPO index///
    dT = 1; % interval 
    cf = 1/8;
    TPIf_long = lanczosfilter(TPI1,dT,cf,[],'low'); % 8 year filtered IPO index
    TPIf(:,s) = TPIf_long(5:end-4);
    TPIfz(:,s) = zscore(TPIf(:,s));  
    Pi{s} = find(TPIfz(:,s) > 0.5);    Ni{s} = find(TPIfz(:,s) < -0.5); 
% AMO index ---------------------------------------------------------------
%     [AMO_zs AMO] = areamean(sstd(:,:,:,s),290:360,91:161,latData); % 70W-0,0-70N
%     dT = 1; % interval
%     cf = 1/8;
%     AMOf8 = lanczosfilter(AMO,dT,cf,[],'low'); % 8 year filtered AMO index
%     AMOf8z(:,s) = zscore(AMOf8);
%     close all;
% positive TPI com SST (each member)
%     sstcom(:,:,s) = nanmean(sstdf(:,:,Pi{s},s),3);
%     Fig = figure('position',[100 100 700 700]);
%     subplot(2,1,1)
%     plot(TPIfz(:,s));
%     hold on
%     plot(zeros(87))
%     set(gca,'YLim',[-3,3]);
% 
%     subplot(2,1,2)
%     caxval = 0.25;
%     contourVARra(sstcom(:,:,s),[-6:caxval/10:6],caxval,0,360,-90,90,lonData,latData)
%     ch = colorbar;
%     set(ch,'Ticks',[-caxval:caxval/5:caxval]);
%     print(Fig,['D:\figures\CESM\Yearly\IPO.TPI.sst\sst_po',num2str(s),'.png'],'-dpng','-r300')
    
end  
%% TPI -> upper ocean tem
Tsubmr1 = reshape(permute(Tsubm(:,1:55,:,:),[3 4 1 2]),[nyr,nmem,360*55]);
clear par11 h01
for s = 1:nmem
    s
parfor i = 1:size(Tsubmr1,3);
    [par11(i,s),h01(i,s),t] = reg1_ttest(TPIfz(:,s),Tsubmr1(:,s,i),0.1,1);
end
end
par11r = reshape(par11,[360,55,nmem]);
h01r = reshape(h01,[360,55,nmem]);
%%
for s = 1:nmem
close all;
    map = par11r(:,:,s);
    Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
    ticks = 0.2;
    contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(1:55),12)
    set_colorbar([0.83 0.08 0.03 0.88],[],4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
    m_line([0:1:360],-46,'linewidth',2,'color','k');
    m_line([0:1:360],-61,'linewidth',2,'color','k');
    m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
    m_line(160,[-61:1:-46],'linewidth',2,'color','k');
    print(Fig,['D:\figures\CESM\Yearly\IPO_reg_Tem0_707m\reg',num2str(s),'.png'],'-dpng','-r300')
end
%%
close all;
map = mean(par11r,3);
h01rm = sum(h01r,3);
h01rm(find(h01rm < 40*0.75)) = 0;
h01rm(find(h01rm >= 40*0.75)) = 1;
ftsz = 12;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 0.15;
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(1:55),12)
testdots(h01rm,[.4 .4 .4],lonData,latData(1:55));
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_Tem0_707m\ComReg_46_61_tested.png'],'-dpng','-r300')

%% meridional mean
% meridional mean
Tem1mr = reshape(permute(Tem1m,[3 4 1 2]),[nyr,nmem,360*37]);
clear par14 h04
for s = 1:nmem
    s
parfor i = 1:size(Tem1mr,3);
    [par14(i,s),h04(i,s),t] = reg1_ttest(TPIfz(:,s),Tem1mr(:,s,i),0.1,1);
end
end
par14r = reshape(par14,[360,nlev,nmem]);
h04r = reshape(h04,[360,nlev,nmem]);
h04rm = sum(h04r,3);
h04rm(find(h04rm < 40*0.75)) = 0;  h04rm(find(h04rm >= 40*0.75)) = 1;
%%
% mixed layer
ncid=netcdf.open('E:\CESM-post\LE\HMXL\HMXLensmean_1920-2005.nc','NC_NOWRITE');
HMXLmmm = ncread('E:\CESM-post\LE\HMXL\HMXLensmean_1920-2005.nc','HMXL');
datadir='E:\CESM-post\LE\HMXL\yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*B20TRC5CNBDRD*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
k=length(filelist);
for s=1:k
    s 
    filename=[datadir,filelist(s).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    HMXL = ncread(filename,'HMXL');
    HMXLd = HMXL - HMXLmmm; % detrend
    HMXLda = HMXLd - mean(HMXLd(:,:,61:86),3); % anomaly
    % 8-yr filter
    dT = 1;  cf = 1/8;
    HMXLsubar = reshape(permute(HMXLda,[3 1 2]),[size(HMXLda,3),length(lonData)*length(latData)]);
    clear HMXLsubadf
    parfor i = 1:size(HMXLsubar,2);
        HMXLsubadf(:,i) = lanczosfilter(HMXLsubar(:,i),dT,cf,[],'low'); % 8 year filtered
    end
    HMXLsubadf = permute(reshape(HMXLsubadf,[size(HMXLda,3),length(lonData),length(latData)]),[2 3 1]);  % -5~707m, 8-year filtered

    save(['MatFile/hmxl_Global/HMXLsubadf',num2str(s),'.mat'],'HMXLsubadf');
end
%% meridional mixed layer
lats = [29:44]; names = '46S-61S'; 
for s = 1:nmem
    s
    load(['MatFile/hmxl_Global/HMXLsubadf',num2str(s),'.mat']);
    HMXLsubadf(:,:,1:4) = []; HMXLsubadf(:,:,end-3:end) = [];
    latm = latmean(HMXLsubadf,lats,latData);
    HMXLr = reshape(permute(latm,[2 1]),[78,360]);
parfor i = 1:size(HMXLr,2);
    [paH(i,s),hH(i,s),t] = reg1_ttest(TPIfz(1:78,s),HMXLr(:,i),0.1,1);
end
end
%%
latmMMM = latmean(HMXLmmm,lats,latData);
HMXLcli = nanmean(latmMMM,2);
HMXLclir = -cat(1,HMXLcli(160:360,:),HMXLcli(1:159,:))/100; % cm -> m
%% com 40
map = nanmean(par14r,3);
map_r = cat(1,map(160:300,:),map(301:360,:),map(1:159,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
h04rm_r = cat(1,h04rm(160:300,:),h04rm(301:360,:),h04rm(1:159,:));
h04rm_r(1:2:end) = 0;
close all;
Fig = figure('position',[100 100 800 400]);
ax = axes('Position',[0.08 0.1 0.72 0.87],'fontsize',14,'box','on');
caxval = 0.1;
[c h] = contourf(lonData_r,-depthData(1:nlev),map_r',[-1:caxval/10:1],'linestyle','none');
caxis([-caxval,caxval])
hold on
line([300 300],[0 -700],'linewidth',1.5,'color','k');
plot(lonData_r,HMXLclir,'linewidth',1.5,'color','k','LineStyle','--')
[xx,yy] = find(h04rm_r==1);
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
maph = nanmean(paH,2)/100; % cm->m
maphr = cat(1,maph(160:360,:),maph(1:159,:));
plot(lonData_r,-maphr,'linewidth',2.5,'color',[.5 .5 .5],'LineStyle','--')
set(gca,'YLim',[-5.5,2.2],'YTick',[-5.5:1.1:2.2],'YColor',[.5 .5 .5]);
yh = ylabel('Mixed Layer Depth (m)')
set(yh,'Rotation',270,'position',[562,-1.73,0])
% print(Fig,['D:\figures\CESM\Yearly\LE\20240130_IPO_SO\IPO_reg_meridionalMean_46S_61S.png'],'-dpng','-r300')
%%
for s = 1:nmem
    map = par14r(:,:,s);
    map_r = cat(1,map(160:300,:),map(301:360,:),map(1:159,:));
    close all;
    Fig = figure('position',[100 100 800 400]);
    caxval = 0.8;
    lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
    [c h] = contourf(lonData_r,-depthData(1:nlev),map_r',[-1:caxval/10:1],'linestyle','none');
    caxis([-caxval,caxval])
    hold on
    % line([220 220 340 340 220],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
    % [ac ah] = contour(lonData,-depthData(1:23),(nanmean(latmc(:,1:23,:),3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
    load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
    bl_re5(18,:) = [0.85 0.18 0.2]; bl_re5(19,:) = [0.8 0.16 0.2];
    bl_re5(20,:) = [0.75 0.15 0.2];
    colormap(bl_re5);
    ch = colorbar;
    set(ch,'Ticks',[-caxval:caxval/5:caxval]);
%     set(gca,'XLim',[150,360+149]);
%     set(gca,'XTick',[150:60:360+149]);
%     set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
    set(gca,'YLim',[-700,0]);
    set(gca,'YTick',[-700:100:0],'FontSize',14);
    set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
    ylabel('Depth (m)')
    % title([' Com.meridional mean temperature (',names,')'])
    print(Fig,['D:\figures\CESM\Yearly\IPO_corr_meridionalMean_46S_61S\IPO.Corr',num2str(s),'.png'],'-dpng','-r300')
end
%% corr Pac with IO+Atl
Tsubr1 = cat(1,Tsubm(160:300,:,:,:),Tsubm(301:360,:,:,:),Tsubm(1:159,:,:,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
%
for s = 1:40
    [spacz(:,s) spac] = areamean(Tsubr1(:,:,:,s),1:141,29:44,latData);
    [siaz(:,s) sia] = areamean(Tsubr1(:,:,:,s),142:360,29:44,latData);
    [r1(s),p1(s),n_eff1] = corr_eff(spacz(:,s),siaz(:,s),0.1);  %  90% confidence
    DI(:,s) = zscore(spac-sia);
    [r2(s),p2(s),n_eff2] = corr_eff(DI(:,s),TPIfz(:,s),0.1);  %  90% confidence    
end
mean(r1)
mean(r2)
%% DI&IPO scatter
close all;
figure('position',[2700 100 500 500])
for s = 1:40
    num_pc1p = find(DI(:,s)>0);
    num_pc1n = find(DI(:,s)<0);
    numhb = [num_pc1n;num_pc1p];
    ch = scatter(DI(numhb,s),TPIfz(numhb,s),[],DI(numhb,s),'filled')
    hold on
    yb0(s,:) = polyfit(DI(:,s),TPIfz(:,s),1);
end
yb = mean(yb0,1);
plot(DI,polyval(yb,DI),'k','linewidth',3)
% scatter(AMOfz(find(DI>0 & DI<0.5)),TPIfz(find(DI>0 & DI<0.5)),'ro')
% scatter(AMOfz(find(DI>-0.5 & DI<0)),TPIfz(find(DI>-0.5 & PC1z<0)),'bo')
set(gca,'XLim',[-3,3]);
set(gca,'XTick',[-3:1:3]);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
text(1,-3,['Slope = ',num2str(roundn(yb(1),-2)),'*'],'fontsize',12)
text(3.2,-0.3,'DI','FontSize',12)
text(-0.6,3.4,'IPO','FontSize',12)
%
new_fig_handle = convert_to_std_coordinate_system(gca,0);
caxis([-2,2]);
load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
cb = colorbar;
set(cb,'Ticks',[-2:0.4:2],'TickLabels',[-2:0.4:2],'Position',[0.88,0.1,0.03,0.8],'FontSize',12);
set(cb.Label,'String','DI','Rotation',270,'position',[4 0 0],'FontSize',12)
saveas(new_fig_handle,['D:\figures\CESM\Yearly\LE\scatter_IPO_DI_46_61S.png'],'png')


%% CESM LE taux tauy
% 有 wind stress curl
% 很奇怪，回归风场符号完全与理论相反，可能是数据出了问题？目前还没
% 找到原因（2023.1.10）.但pacemaker试验结果正常。
% 且回归 taux tauy curlz比观测大了大概3倍，奇怪（2023.5.18）
for s = 1:40
    load(['MatFile/taux_Global/tauxadf',num2str(s),'.mat']);
    tauxadf(:,:,1:4) = []; tauxadf(:,:,end-3:end) = [];
    taux(:,:,:,s) = double(tauxadf);
    load(['MatFile/tauy_Global/tauyadf',num2str(s),'.mat']);
    tauyadf(:,:,1:4) = []; tauyadf(:,:,end-3:end) = [];
    tauy(:,:,:,s) = double(tauyadf);
end
%
lats = 1:55;
tauxr1 = reshape(permute(taux(:,lats,:,:),[3 4 1 2]),[nyr,nmem,nlon*length(lats)]);
tauyr1 = reshape(permute(tauy(:,lats,:,:),[3 4 1 2]),[nyr,nmem,nlon*length(lats)]);
clear par12 h02 par13 h03
for s = 1:nmem
parfor i = 1:size(tauxr1,3);
    [par12(i,s),h02(i,s),t] = reg1_ttest(TPIfz(:,s),tauxr1(:,s,i),0.1,1);
    [par13(i,s),h03(i,s),t] = reg1_ttest(TPIfz(:,s),tauyr1(:,s,i),0.1,1);
end
end
par12r = -reshape(par12,[nlon,length(lats),nmem]);
h02r = reshape(h02,[nlon,length(lats),nmem]);
par13r = -reshape(par13,[nlon,length(lats),nmem]);
h03r = reshape(h03,[nlon,length(lats),nmem]);
mapu = mean(par12r,3);
mapv = mean(par13r,3);
%% taux
close all;
ftsz = 12; ticks1 = 5;
h02rm = sum(h02r,3);
h02rm(find(h02rm < 40*0.75)) = 0; h02rm(find(h02rm >= 40*0.75)) = 1;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(mapu*10^3,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData(lats),12)
testdots(h02rm,[.4 .4 .4],lonData,latData(lats));
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'TAUX  (10^-^3 Pa)',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_taux\ComReg_46_61S_tested.png'],'-dpng','-r300')
%% tauy
close all;
ftsz = 12; ticks1 = 5;
h03rm = sum(h03r,3);
h03rm(find(h03rm < 40*0.75)) = 0; h03rm(find(h03rm >= 40*0.75)) = 1;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(mapv*10^3,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData(lats),12)
testdots(h03rm,[.4 .4 .4],lonData,latData(lats));
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'TAUY  (10^-^3 Pa)',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_tauy\ComReg_46_61S_tested.png'],'-dpng','-r300')
%% for s = 1:110
close all;
Fig = figure('position',[100 100 800 400]);
ftsz = 12;
% m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
caxval = 0.4;
[c h] = m_contourf(lonData,latData(lats),mapu',[-caxval*10:caxval/10:caxval*10],'linestyle','none');
caxis([-caxval,caxval])
hold on
load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-caxval:caxval/5:caxval],'position',[0.78 0.09 0.03 0.85],'fontsize',12);
m_coast('linewidth',1,'color','k');
m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
d = 8; dd = 0.02;
[x,y] = meshgrid(lonData,latData(lats));
umap = mapu.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,mapv(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
% print(Fig,['D:\figures\CESM\Yearly\IPO_reg_Tau\comReg.png'],'-dpng','-r300')
%% wind stress curl
for s = 1:40
    load(['MatFile/curlz_Global/curlzadf',num2str(s),'.mat']);
    curlzadf(:,:,1:4) = []; curlzadf(:,:,end-3:end) = [];
    curlz(:,:,:,s) = double(curlzadf);
end
lats = 1:55;
curlzr1 = reshape(permute(curlz(:,lats,:,:),[3 4 1 2]),[nyr,nmem,nlon*length(lats)]);
clear par14 h04
for s = 1:nmem
parfor i = 1:size(curlzr1,3);
    [par14(i,s),h04(i,s),t] = reg1_ttest(TPIfz(:,s),curlzr1(:,s,i),0.1,1); 
end
end
par14r = -reshape(par14,[nlon,55,nmem]);
h04r = reshape(h04,[nlon,length(lats),nmem]);
%%
for s = 1:nmem
close all;
    map = par14r(:,:,s)*10^8;
    Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
    ticks = 2;
    contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(1:55),12)
    set_colorbar([0.83 0.08 0.03 0.88],'10^-^8',4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
    m_line([0:1:360],-46,'linewidth',2,'color','k');
    m_line([0:1:360],-61,'linewidth',2,'color','k');
    m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
    m_line(160,[-61:1:-46],'linewidth',2,'color','k');
    print(Fig,['D:\figures\CESM\Yearly\IPO_reg_curlz\reg',num2str(s),'.png'],'-dpng','-r300')
end
%%
close all; ftsz = 12;
map = mean(par14r,3)*10^9;
h04rm = sum(h04r,3);
h04rm(find(h04rm < 40*0.75)) = 0; h04rm(find(h04rm >= 40*0.75)) = 1;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 8; 
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(1:55),12)
testdots(h04rm,[.4 .4 .4],lonData,latData(lats));
set_colorbar([0.83 0.08 0.03 0.88],'Wind Stress Curl  (10^-^9 Pa m^-^1)',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
hold on
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
d = 5; dd = 0.0006;
[x,y] = meshgrid(lonData,latData(lats));
umap = mapu.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,mapv(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1,'MaxHeadSize',5) 
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_curlz\ComReg_46_61_tested.png'],'-dpng','-r300')
%% IPO reg. slp
lats = 1:55;
slpr1 = reshape(permute(slp(:,lats,:,:),[3 4 1 2]),[nyr,nmem,nlon*length(lats)]);
clear par15 h02 h05
for s = 1:nmem
parfor i = 1:size(slpr1,3);
    [par15(i,s),h05(i,s),t] = reg1_ttest(TPIfz(:,s),slpr1(:,s,i),0.1,1)
end
end
par15r = reshape(par15,[nlon,length(lats),nmem]);
h05r = reshape(h05,[nlon,length(lats),nmem]);
%%
close all;
map = mean(par15r,3);
h05rm = sum(h05r,3);
h05rm(find(h05rm < 40*0.75)) = 0; h05rm(find(h05rm >= 40*0.75)) = 1;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 60;
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(lats),12)
testdots(h05rm,[.4 .4 .4],lonData,latData(lats));
set_colorbar([0.83 0.08 0.03 0.88],'SLP  (Pa)',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
hold on
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
d = 5; dd = 0.0006;
[x,y] = meshgrid(lonData,latData(lats));
umap = mapu.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,mapv(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1,'MaxHeadSize',5) 
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_slp\ComReg_46_61_tested.png'],'-dpng','-r300')
% 
%% IPO reg. nhflx
%
lats = 1:120;
nhflxr1 = reshape(permute(nhflxdf(:,lats,:,:),[3 4 1 2]),[nyr,nmem,nlon*length(lats)]);
clear par16 h16 
for s = 1:nmem
parfor i = 1:size(nhflxr1,3);
    [par16(i,s),h16(i,s),t] = reg1_ttest(TPIfz(:,s),nhflxr1(:,s,i),0.1,1)
end
end
par16r = reshape(par16,[nlon,length(lats),nmem]);
h16r = reshape(h16,[nlon,length(lats),nmem]);
%%
close all; ftsz = 12;
map = mean(par16r,3);
max(map,[],'all')
min(map,[],'all')
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 6;
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(lats),12)
set_colorbar([0.83 0.08 0.03 0.88],'NHFLX  (W m^-^2)',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
hold on
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_nhflx\ComReg_46_61S_tested.png'],'-dpng','-r300')
%%
close all;    
for s = 1:40
map = mean(par16r(:,:,s),3);
ticks = 4;
max(map,[],'all')
Fig = figure('position',[700 100 600 300]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
    hold on
    m_contourf(lonData,latData(lats),map',[-ticks*5:ticks/10:ticks*5],'linestyle','none');
    caxis([-ticks,ticks]);
    load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
    print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_nhflx\Reg',num2str(s),'.png'],'-dpng','-r300')
end
%% IPO reg. sst
lats = 1:120;
sstr1 = reshape(permute(sstdf(:,lats,:,:),[3 4 1 2]),[nyr,nmem,nlon*length(lats)]);
clear par17 h17 
for s = 1:nmem
parfor i = 1:size(sstr1,3);
    [par17(i,s),h17(i,s),t] = reg1_ttest(TPIfz(:,s),sstr1(:,s,i),0.1,1)
end
end
par17r = reshape(par17,[nlon,length(lats),nmem]);
h17r = reshape(h17,[nlon,length(lats),nmem]);
max(par17r,[],'all')
min(par17r,[],'all')
%
close all;
map = mean(par17r,3);
h17rm = sum(h17r,3);
h17rm(find(h17rm < 40*0.75)) = 0; h17rm(find(h17rm >= 40*0.75)) = 1;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 0.2;
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(lats),12)
testdots(h17rm,[.4 .4 .4],lonData,latData(lats));
set_colorbar([0.83 0.08 0.03 0.88],'SST (K)',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
hold on
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_sst\ComReg.46_61S_tested.png'],'-dpng','-r300')
%% IPO reg zonal mean (SST nhflx curlz slp taux)
clear var par11* par12* map1 map2
for s = 1:40
    [par11s(s,:) par12s(s,:) h11s(s,:) h12s(s,:)] = reg_zonalmean(sstdf(:,:,:,s),TPIfz(:,s));
    [par11n(s,:) par12n(s,:) h11n(s,:) h12n(s,:)] = reg_zonalmean(nhflxdf(:,:,:,s),TPIfz(:,s));
    [par11c(s,:) par12c(s,:) h11c(s,:) h12c(s,:)] = reg_zonalmean(curlz(:,:,:,s),TPIfz(:,s));
    [par11t(s,:) par12t(s,:) h11t(s,:) h12t(s,:)] = reg_zonalmean(taux(:,:,:,s),TPIfz(:,s));
    [par11p(s,:) par12p(s,:) h11p(s,:) h12p(s,:)] = reg_zonalmean(slp(:,:,:,s),TPIfz(:,s));
end
%% sst
h11ss = sum(h11s,1); h11ss(find(h11ss < 40*0.75)) = 0; h11ss(find(h11ss >= 40*0.75)) = 1; 
h12ss = sum(h12s,1); h12ss(find(h12ss < 40*0.75)) = 0; h12ss(find(h12ss >= 40*0.75)) = 1; 
close all; map1 = mean(par11s,1); map2 = mean(par12s,1);
zonalmean_plot(map1,map2,h11ss,h12ss,latData,[1:90],[-0.1,0.15],[-0.1:0.05:0.15],'Sea Surface Temperature (K)','SST');
% nhflx
h11ns = sum(h11n,1); h11ns(find(h11ns < 40*0.75)) = 0; h11ns(find(h11ns >= 40*0.75)) = 1; 
h12ns = sum(h12n,1); h12ns(find(h12ns < 40*0.75)) = 0; h12ns(find(h12ns >= 40*0.75)) = 1; 
close all; map1 = mean(par11n,1); map2 = mean(par12n,1);
zonalmean_plot(map1,map2,h11ns,h12ns,latData,[1:90],[-1.5,1.5],[-1.5:0.5:1.5],'Net Heat Flux (W/m^-^2)','NHFLX');
%% curlz
h11cs = sum(h11c,1); h11cs(find(h11cs < 40*0.75)) = 0; h11cs(find(h11cs >= 40*0.75)) = 1; 
h12cs = sum(h12c,1); h12cs(find(h12cs < 40*0.75)) = 0; h12cs(find(h12cs >= 40*0.75)) = 1; 
close all; map1 = -mean(par11c,1); map2 = -mean(par12c,1);
zonalmean_plot(map1,map2,h11cs,h12cs,latData,[1:90],[-8,8]*10^-9,[-8:2:8]*10^-9,'Wind Stress Curl (Pa/m)','CURLZ');
%% taux 
h11ts = sum(h11t,1); h11ts(find(h11ts < 40*0.75)) = 0; h11ts(find(h11ts >= 40*0.75)) = 1; 
h12ts = sum(h12t,1); h12ts(find(h12ts < 40*0.75)) = 0; h12ts(find(h12ts >= 40*0.75)) = 1; 
close all; map1 = mean(-par11t,1); map2 = mean(-par12t,1);
zonalmean_plot(map1,map2,h11ts,h12ts,latData,[1:90],[-8,4]*10^-3,[-8:2:4]*10^-3,'Zonal Wind Stress (Pa)','TAUX');
% slp
h11ps = sum(h11p,1); h11ps(find(h11ps < 40*0.75)) = 0; h11ps(find(h11ps >= 40*0.75)) = 1; 
h12ps = sum(h12p,1); h12ps(find(h12ps < 40*0.75)) = 0; h12ps(find(h12ps >= 40*0.75)) = 1; 
close all; map1 = mean(par11p,1); map2 = mean(par12p,1);
zonalmean_plot(map1,map2,h11ps,h12ps,latData,[1:90],[-32,55],[-30:10:50],'Sea Level Pressure (Pa)','SLP');
%% MOHT calculate (Ekman layer)
f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:37)-depthData(1:36);
levdist(1) = 4; levdist(2:37) = dweit; % z distance
datadir1='E:\CESM-post\JC-data\LE\Temperature\his_sub_yr1x1\'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='E:\CESM-post\JC-data\LE\TAUX\historical\taux_yr_1x1\'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1);
clear OHC Vy
for s=1:k
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Taux = ncread(filename2,'TAUX');
    % Ekman layer MOHT
    OHC(:,:,:,s) = permute(ro*cp*sum((Temp(:,1:90,1:37,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3])/700;
    My = -Taux(:,1:90,:)./f'; % Ekman 质量输运
    Vy(:,:,:,s) = My/ro;
end
clear Temp Taux My filename1 filename2 datadir1 datadir2 filelist1 filelist2
%%
OHCT = OHC.*Vy;
OHCTw = OHCT.*londist';
OHCTw_r = cat(1,OHCTw(160:300,:,:,:),OHCTw(301:360,:,:,:),OHCTw(1:159,:,:,:));
MOHTe1 = permute(nansum(OHCTw_r(1:141,:,:,:),1),[2 3 4 1]); % Pac unit:W
MOHTe2 = permute(nansum(OHCTw_r(142:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
%
MOHT = permute(nansum(mean(OHCTw_r,4),1),[2 3 1]); % Pac unit:W
MOHT(85:end,:) = nan;
close all;
Fig = figure('position',[100 100 600 300])
map = -mean(MOHT,2)/10^15; % 都是taux符号的错
% map(85:96,:) = nan;
plot(latData(1:90),map,'linewidth',1.5)
set(gca,'XLim',[-90,90],'YLim',[-2,2],'XGrid','on','YGrid','on');
set(gca,'XTick',[-80:20:80],'YTick',[-2:0.5:2]);
ylabel('MOHT (PW)');xlabel('Latitude')
% print(Fig,['D:\figures\CESM\Yearly\LE\MOHT_Ekman_climotology.png'],'-dpng','-r300') 
% 比观测大几倍，但又不到十倍，很奇怪,可能是因为taux
%% MOHT anomaly,detrend, and 8-yr filter
clear MOHT1a MOHT2a MOHT1ad MOHT1adf MOHT2ad MOHT2adf
for s = 1:40
MOHT1a(:,:,s) = MOHTe1(:,:,s)-mean(MOHTe1(:,:,s),2);
MOHT2a(:,:,s) = MOHTe2(:,:,s)-mean(MOHTe2(:,:,s),2);
dT = 1; cf = 1/8;
for i = 1:90
    MOHT1ad(:,i,s) = detrend(MOHT1a(i,:,s)'); % linear detrend
    MOHT1adf(:,i,s) = lanczosfilter(MOHT1ad(:,i,s),dT,cf,[],'low'); % 8 year filtered
    MOHT2ad(:,i,s) = detrend(MOHT2a(i,:,s)'); % linear detrend
    MOHT2adf(:,i,s) = lanczosfilter(MOHT2ad(:,i,s),dT,cf,[],'low'); % 8 year filtered
end
end
MOHT1adf(1:4,:,:) = []; MOHT2adf(1:4,:,:) = [];
MOHT1adf(end-3:end,:,:) = []; MOHT2adf(end-3:end,:,:) = [];
%
clear par7 h07 par8 h08 
for s = 1:40;
index = TPIfz(:,s);
for i = 1:size(MOHT1adf,2);
    [par7(i,s),h07(i,s),t] = reg1_ttest(index,MOHT1adf(:,i,s),0.1,1);
    [par8(i,s),h08(i,s),t] = reg1_ttest(index,MOHT2adf(:,i,s),0.1,1);
end
end
%% 
% for s = 1:40
close all;
h07s = sum(h07,2); h07s(find(h07s < 40*0.75)) = 0; h07s(find(h07s >= 40*0.75)) = 1; 
h08s = sum(h08,2); h08s(find(h08s < 40*0.75)) = 0; h08s(find(h08s >= 40*0.75)) = 1; 
map7 = -mean(par7,2)/10^15;
map8 = -mean(par8,2)/10^15;
Fig = figure('position',[100 100 800 400]);
plot(latData(1:90),map7,'r-','linewidth',1.5)
hold on
plot(latData(1:90),map8,'b-','linewidth',1.5)
plot(latData(1:90),zeros(1,length(latData(1:90))),'k','linewidth',1.5)
plot(-61*ones(20),[-.02:.04/19:.02],'k--','linewidth',1.5);
plot(-46*ones(20),[-.02:.04/19:.02],'k--','linewidth',1.5);
scatter(latData(find(h07s == 1)),map7(find(h07s == 1)),'o','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h08s == 1)),map8(find(h08s == 1)),'o','MarkerEdgeColor','b','linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
set(gca,'YLim',[-.02,.02]);
set(gca,'YTick',yticks,'FontSize',14);
ylabel('MOHT (PW)')

% yyaxis right
% map11 = -mean(par11t,1);
% map12 = -mean(par12t,1);
% plot(latData(1:90),map11,'r--','linewidth',1.5)
% plot(latData(1:90),map12,'b--','linewidth',1.5)
% scatter(latData(find(h11ts == 1)),map11(find(h11ts == 1)),'*','MarkerEdgeColor','r','linewidth',1)
% scatter(latData(find(h12ts == 1)),map12(find(h12ts == 1)),'*','MarkerEdgeColor','b','linewidth',1)
% set(gca,'YLim',[-8,8]*10^-3,'ycolor','k');
% set(gca,'YTick',[-8:2:8]*10^-3,'FontSize',14);
% ylabel('TAUX (Pa)','color','k','rotation',270,'position',[-26.5 0 -1]);

% end
% print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_MOHT_Ekman\Com_tested.png'],'-dpng','-r300')
%% MOHT calculate (0-700m) 
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:37)-depthData(1:36);
levdist(1) = 4; levdist(2:37) = dweit; % z distance
datadir1='E:\CESM-post\JC-data\LE\Temperature\his_sub_yr1x1\'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*.nc']); %指定批量数据的类型
ncdisp([datadir1,filelist1(1).name]);
datadir2='E:\CESM-post\JC-data\LE\VVEL\his_vvel_yr1x1\'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*.nc']); %指定批量数据的类型
ncdisp([datadir2,filelist2(1).name]);
k=length(filelist1);
clear OHCz*
%
for s=1:40
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp = ncread(filename1,'TEMP');
    filename2=[datadir2,filelist2(s).name];
    ncid=netcdf.open(filename2,'NC_NOWRITE');
    Vvel = ncread(filename2,'VVEL')/100; % cm/s -> m/s
    % 0-700m MOHT
    TemVel = Temp.*Vvel(:,1:90,:,:);
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,1:37,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(160:300,:,:,:),OHCzy(301:360,:,:,:),OHCzy(1:159,:,:,:));
    MOHTp1v(:,:,s) = permute(nansum(OHCzy_r(1:141,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHTia1v(:,:,s) = permute(nansum(OHCzy_r(142:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHTall1v(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
save("MatFile\MOHT_160E_700m_VVEL\MOHTall1v",'MOHTall1v');
% save("MatFile\MOHT_160E_700m_VVEL\MOHTa1v",'MOHTa1v');
save("MatFile\MOHT_160E_700m_VVEL\MOHTia1v",'MOHTia1v');
save("MatFile\MOHT_160E_700m_VVEL\MOHTp1v",'MOHTp1v');
%%
close all;
Fig = figure('position',[100 100 600 300])
map = mean(mean(MOHTall1v,2),3)/10^15;
% map(85:96,:) = nan;
plot(latData(1:90),map,'linewidth',1.5)
set(gca,'XLim',[-90,90],'YLim',[-2,2],'XGrid','on','YGrid','on');
set(gca,'XTick',[-80:20:80],'YTick',[-2:0.5:2]);
ylabel('MOHT (PW)');xlabel('Latitude')
% print(Fig,['D:\figures\CESM\Yearly\LE\MOHT_Tv_climotology.png'],'-dpng','-r300')  
% 跟实际差距较大
%% MOHT anomaly,detrend, and 8-yr filter
load("MatFile\MOHT_160E_700m_VVEL\MOHTall1v",'MOHTall1v');
load("MatFile\MOHT_160E_700m_VVEL\MOHTa1v",'MOHTa1v');
load("MatFile\MOHT_160E_700m_VVEL\MOHTia1v",'MOHTia1v');
load("MatFile\MOHT_160E_700m_VVEL\MOHTp1v",'MOHTp1v');
%%
clear MOHT71a MOHT72a MOHT71ad MOHT71adf MOHT72ad MOHT72adf
for s = 1:40
MOHT71a(:,:,s) = MOHTp1v(:,:,s)-mean(MOHTp1v(:,:,s),2);
MOHT72a(:,:,s) = MOHTia1v(:,:,s)-mean(MOHTia1v(:,:,s),2);
dT = 1; cf = 1/8;
for i = 1:90
    MOHT71ad(:,i,s) = detrend(MOHT71a(i,:,s)'); % linear detrend
    MOHT71adf(:,i,s) = lanczosfilter(MOHT71ad(:,i,s),dT,cf,[],'low'); % 8 year filtered
    MOHT72ad(:,i,s) = detrend(MOHT72a(i,:,s)'); % linear detrend
    MOHT72adf(:,i,s) = lanczosfilter(MOHT72ad(:,i,s),dT,cf,[],'low'); % 8 year filtered
end
end
MOHT71adf(1:4,:,:) = []; MOHT72adf(1:4,:,:) = [];
MOHT71adf(end-3:end,:,:) = []; MOHT72adf(end-3:end,:,:) = [];
%
clear par77 h077 par87 h087 
for s = 1:40
index = TPIfz(:,s);
for i = 1:size(MOHT71adf,2);
    [par77(i,s),h077(i,s),t] = reg1_ttest(index,MOHT71adf(:,i,s),0.1,1);
    [par88(i,s),h088(i,s),t] = reg1_ttest(index,MOHT72adf(:,i,s),0.1,1);
end
end
%%
% for s = 1:40
close all;
h077s = sum(h077,2); h077s(find(h077s < 40*0.75)) = 0; h077s(find(h077s >= 40*0.75)) = 1; 
h088s = sum(h088,2); h088s(find(h088s < 40*0.75)) = 0; h088s(find(h088s >= 40*0.75)) = 1; 
Fig = figure('position',[100 100 800 400]);
map77 = mean(par77,2)/10^15;
map88 = mean(par88,2)/10^15;
p1 = plot(latData(1:90),map77,'-','color',[0.90,0.51,0.31],'linewidth',1.5)
hold on
p2 = plot(latData(1:90),map88,'-','color',[0.17,0.73,0.69],'linewidth',1.5)
p3 = plot(latData(1:90),map7,'r-','linewidth',1.5)
p4 = plot(latData(1:90),map8,'b-','linewidth',1.5)
plot(latData(1:90),zeros(1,length(latData(1:90))),'k','linewidth',1.5)
plot(-61*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
plot(-46*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
scatter(latData(find(h07s == 1)),map7(find(h07s == 1)),'o','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h08s == 1)),map8(find(h08s == 1)),'o','MarkerEdgeColor','b','linewidth',1)
scatter(latData(find(h077s == 1)),map77(find(h077s == 1)),'o','MarkerEdgeColor',[0.90,0.51,0.31],'linewidth',1)
scatter(latData(find(h088s == 1)),map88(find(h088s == 1)),'o','MarkerEdgeColor',[0.17,0.73,0.69],'linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
set(gca,'YLim',[-.025,.032],'YTick',[-.03:.01:.03]);
set(gca,'YTick',[-.03:.01:.03],'FontSize',14);
ylabel('MOHT (PW)')
% end
% print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_MOHT_0_700m\Com_MOHT_160E_700m_Ekman_tested.png'],'-dpng','-r300')
%%
close all;
Fig = figure('position',[100 100 800 400]);
p1 = plot(latData(1:90),zeros(90,1),'-','color',[0.90,0.51,0.31],'linewidth',1.5)
hold on
p2 = plot(latData(1:90),-0.015*ones(90,1),'-','color',[0.17,0.73,0.69],'linewidth',1.5)
p3 = plot(latData(1:90),-0.01*ones(90,1),'r-','linewidth',1.5)
p4 = plot(latData(1:90),-0.02*ones(90,1),'b-','linewidth',1.5)
set(gca,'XLim',[-70,-35]);set(gca,'XTick',[-70:10:-35]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',[-.025,.032],'YTick',[-.03:.01:.03]);
ylabel('MOHT (PW)')
legend([p1,p2,p3,p4],'Pacific Eulerian mean','Atlantic & Indian Ocean Eulerian mean','Pacific Ekman',...
    'Atlantic & Indian Ocean Ekman','Location','north','Orientation','vertical','position',[0.27,0.6,0.18,0.3])
legend('boxoff')
print(Fig,['D:\figures\CESM\Yearly\LE\20240130_IPO_SO\Com_legend.png'],'-dpng','-r300')
%% dMOHT/dy



%% scatter time series
for s = 1:40;
    lats = [29:44];
    var = nhflxdf(:,:,:,s); varname = 'NHFLX'; titlename = 'Net Heat Flux';
    var_r = cat(1,var(160:300,:,:),var(301:360,:,:),var(1:159,:,:));
    lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
    [tz1 t1] = areamean(var_r,1:141,lats,latData); % Pac
    [tz2 t2] = areamean(var_r,142:360,lats,latData); % Atl & IO
    %
    % tx1(:,s) = DI(1:72); tx2(:,s) = zscore(t1-t2); pocname = 'DI';
    tx1(:,s) = TPIfz(:,s); tx2(:,s) = zscore(t1-t2); pocname = 'IPO';
    [r(s),p(s),n_eff2] = corr_eff(tx1(:,s),tx2(:,s),0.1);  %  90% confidence
end
%%
close all;
Fig = figure('position',[700 100 400 400]);
for s = 1:40
    xlimv = [-4,4];
    ylimv = [-5,5];
    scatter(tx1(:,s),tx2(:,s),'o','markeredgeColor',[0,46,166]/255,'LineWidth',1);
    hold on
    set(gca,'XLim',xlimv,'xtick',[-4:1:4]);
    set(gca,'FontSize',12);
    set(gca,'YLim',ylimv,'ytick',[-5:1:5]);
    xlabel('IPO','FontSize',12),ylabel(titlename,'FontSize',12);
end
%
for s = 1:40
    yb(s,:) = polyfit(tx1(:,s),tx2(:,s),1);
    plot(tx1(:,s),polyval(yb(s,:),tx1(:,s)),'color',[.7 .7 .7],'linewidth',1)
end
for s = 1:40
    plot(tx1(:,s),polyval(mean(yb,1),tx1(:,s)),'k','linewidth',4)
end
text(1.3,-3.5,['Corr = ',num2str(roundn(mean(r),-2))],'fontsize',12)
text(1.3,-4,['Slope = ',num2str(roundn(mean(yb(:,1),1),-2))],'fontsize',12)
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_',varname,'_scatter_46_61S_tested.png'],'-dpng','-r300')





function contourVARra(var,val,ctr,lon1,lon2,lat1,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[lat1 lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
end
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
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
function [Fig] = diffdepth(Temcom,nlevel,k,lonData,latData,names,val,ticks)
    Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
    contourVARra(Temcom(:,:,nlevel),[-val:ticks/10:val],1,0,360,-90,90,lonData,latData);
    caxis([-ticks,ticks]);
    ch = colorbar;
    set(ch,'Ticks',[-ticks:ticks/5:ticks]);
    % testdots(h05,'k',lonData,latData)
    title([' Com.T',' ',names,'m']);
end
function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
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
function [par11 par12 h01 h02] = reg_zonalmean(nhflx,index)
lats = [1:90];
lonmh1 = permute(nanmean(nhflx(160:300,lats,:),1),[2 3 1]);
nhflx_r = cat(1,nhflx(300:360,:,:),nhflx(1:299,:,:));
lonmh2 = permute(nanmean(nhflx_r(1:220,lats,:),1),[2 3 1]);
% corr with PC1
var1 = lonmh1';
size(var1)
size(index)
var2 = lonmh2';
clear par11 h01 par12 h02
for i = 1:size(var1,2);
    [par11(i),h01(i),t] = reg1_ttest(index,var1(:,i),0.1,1); % filtered
    [par12(i),h02(i),t] = reg1_ttest(index,var2(:,i),0.1,1); % filtered
end
max(par11,[],'all')
min(par11,[],'all')
end

function zonalmean_plot(par11,par12,h01,h02,latData,lats,ylimval,yticks,titlename,pngname)
%
close all;
Fig = figure('position',[100 100 800 400]);
plot(latData(lats),par11,'r','linewidth',1.5)
hold on
plot(latData(lats),par12,'b','linewidth',1.5)
plot(latData(lats),zeros(1,length(lats)),'k','linewidth',1.5)
plot(-61*ones(20),[ylimval(1):(ylimval(2)-ylimval(1))/19:ylimval(2)],'k--','linewidth',1.5);
plot(-46*ones(20),[ylimval(1):(ylimval(2)-ylimval(1))/19:ylimval(2)],'k--','linewidth',1.5);
scatter(latData(find(h01 == 1)),par11(find(h01 == 1)),'*','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h02 == 1)),par12(find(h02 == 1)),'*','MarkerEdgeColor','b','linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
set(gca,'YLim',ylimval);
set(gca,'YTick',yticks,'FontSize',14);
ylabel(titlename)
print(Fig,['D:\figures\CESM\Yearly\LE\IPO_reg_zonalmean_',pngname,'\com_46_61S_tested.png'],'-dpng','-r300')
end
function testdots(h,clor,lonData,latData)
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',3,'color',clor); % F-test dots
end