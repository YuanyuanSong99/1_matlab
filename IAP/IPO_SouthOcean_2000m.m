%% EOF
clc,clear,close all;
addpath G:\1_matlab\help;
load('MatFile/Tempadlf');
Tempadlf(:,:,:,1:4) = [];
Tempadlf(:,:,:,end-3:end) = []; % 1944-2016
load("MatFile\TPI_filtered.mat");
TPI_filtered(1:4) = []; TPI_filtered(end-3:end) = [];
TPIfz = zscore(TPI_filtered);
load("MatFile\AMOf8.mat");
AMOf8(1:4) = []; AMOf8(end-3:end) =[];
AMOfz = zscore(AMOf8);
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
%%  load u v 
filename = 'G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\10mSI-MO\uwnd.10m.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
uwnd = ncread(filename,'uwnd');
uwnd = permute(mean(reshape(uwnd(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
ua = double(uwnd - mean(uwnd(:,:,42:71),3)); % remove climatology from 1981 to 2010
filename = 'G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\10mSI-MO\vwnd.10m.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
vwnd = ncread(filename,'vwnd');
vwnd = permute(mean(reshape(vwnd(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
va = double(vwnd - mean(vwnd(:,:,42:71),3)); % remove climatology from 1981 to 2010
uar = permute(ua,[3 1 2]);
var = permute(va,[3 1 2]);
dT = 1; cf = 1/8;
for i = 1:360
    for j = 1:181
    uad(:,i,j) = detrend(uar(:,i,j)); % linear detrend
    uadf(:,i,j) = lanczosfilter(uad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    vad(:,i,j) = detrend(var(:,i,j)); % linear detrend
    vadf(:,i,j) = lanczosfilter(vad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
uadf = permute(uadf,[2 3 1]);
uadf(:,:,1:4) = []; 
vadf = permute(vadf,[2 3 1]);
vadf(:,:,1:4) = [];
%% EOF
lats = 1:90;
dweit = depthData(2:end)-depthData(1:end-1); % depth weight
Tsub = permute(nansum(Tempadlf(:,lats,1,:)+Tempadlf(:,lats,2:end,:).*(permute(dweit,[3 2 1])),3)/2000,[1 2 4 3]); % 0-2000m
[eof_maps,pc,expvar]=eofnan(Tsub); 
%
close all;
ftsz = 12; ticks = 0.2;
% Fig = figure('position',[10 50 550 400]);
% ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
% contourfS(-eof_maps(:,:,1)*10,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
% m_line([0:1:360],-45,'linewidth',2,'color','k'); 
% m_line([0:1:360],-60,'linewidth',2,'color','k'); 
% m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
% m_line(150,[-60:1:-45],'linewidth',2,'color','k'); 
% set_colorbar([0.83 0.08 0.03 0.88],[],4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
Fig = figure('position',[100 100 600 300]);
contourVARra(-eof_maps(:,:,1)*10,[-0.4:0.01:0.4],ticks,0,360,-90,90,lonData,latData(lats))
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% print(Fig,['G:\figures\IAP\Yearly\20221124_IPO_SouthernOcean_2000m\Tem0_2000m_EOFmap_90S_90N.png'],'-dpng','-r300')
%% PC1 .reg. zonal mean 150E-60W
lonm1 = permute(nanmean(Tempadlf(150:300,lats,:,:),1),[2 3 4 1]);
% lonm1m = permute(nanmean(Temp(150:300,lats,1:27,:),1),[2 3 4 1]); % climotology
% corr with PC1
var = permute(lonm1,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(pc(1,:),var(:,i,j),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
%
% close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData,par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% [ac ah] = contour(latData(lats),-depthData(1:27),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:100:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.PC1.zonal mean temperature (South Pac)')
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\Tem0-700m_EOFmap_0_90S_Pac.png'],'-dpng','-r300')
%% PC1 .reg. zonal mean 30W-150E 
Temp_r = cat(1,Tempadlf(300:360,:,:,:),Tempadlf(1:299,:,:,:));
lonm1 = permute(nanmean(Temp_r(1:210,lats,1:27,:),1),[2 3 4 1]);
% corr with PC1
var = permute(lonm1,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(pc(1,:),var(:,i,j),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
%
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:27),par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% [ac ah] = contour(latData,-depthData(1:27),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.PC1.zonal mean temperature (South IO+Atl)')
print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\Tem0-700m_EOFmap_0_90S_IO+Atl.png'],'-dpng','-r300')
%% PC1 .reg. zonal mean 0-360
lonm1 = permute(nanmean(Tempadlf(:,lats,1:27,:),1),[2 3 4 1]);
% lonm1m = permute(nanmean(Temp(150:300,lats,1:27,:),1),[2 3 4 1]); % climotology
% corr with PC1
var = permute(lonm1,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(pc(1,:),var(:,i,j),0.05); % filtered
    end
end
max(par11,[],'all')
min(par11,[],'all')
%
% close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:27),par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% [ac ah] = contour(latData(lats),-depthData(1:27),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.PC1.zonal mean temperature (Southern Ocean)')
print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\Tem0-700m_EOFmap_0_90S_zm_0_360.png'],'-dpng','-r300')
%% PC1 .reg. meridional mean (45S-60S) 
lats = [1:30]; names = '45S-60S'; index = pc(1,:);
latm = latmean(Tempadlf,lats,latData);
% latmc = latmean(Temp,lats,latData); % climotology
% corr with TPI
var = permute(latm,[3 1 2]);
clear par11 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = corr_eff(index,var(:,i,j),0.05);
    end
end
max(par11,[],'all')
min(par11,[],'all')
%
close all;
par11_r = cat(1,par11(150:300,:),par11(301:360,:),par11(1:149,:));
lonData_r = [lonData(150:300);lonData(301:360);lonData(1:149)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,par11_r',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% line([300 300],[0 -700],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[150,360+149]);
set(gca,'XTick',[150:60:360+149]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-2000,0]);
set(gca,'YTick',[-2000:500:0],'FontSize',14);
set(gca,'YTickLabel',[2000:-500:0],'FontSize',14);
ylabel('Depth (m)')
title(['Corr.PC1.meridional mean temperature (',names,')'])
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\Tem0-700m_EOFmap_0_90S_lonm.png'],'-dpng','-r300')

%% PC1>0.5 u v 
ucom1 = nanmean(uadf(:,:,num_pc1p),3);
vcom1 = nanmean(vadf(:,:,num_pc1p),3);
ucom1(:,181)=[]; vcom1(:,181)=[];
ucom1(:,1:10) = nan; ucom1(:,170:180)=nan;
close all;
ticks = 0.6; ftsz = 12;
Fig = figure('position',[100 100 800 400]);
% m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
hold on
[c h] = m_contourf(lonData,latData,ucom1(:,1:180)',[-0.9:0.03:0.9],'linestyle','none');
caxis([-ticks,ticks])
hold on
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks],'position',[0.78 0.09 0.03 0.85],'fontsize',12);
m_coast('linewidth',1,'color','k');
m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
d = 6; dd = 0.03;
[x,y] = meshgrid(lonData,latData);
umap = ucom1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,vcom1(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\PC1p_comUV.png'],'-dpng','-r300')
%%  time series 
pc1z = zscore(pc(1,:));
num_pc1p = find(pc1z>0.5);
num_pc1n = find(pc1z<-0.5);
close all;
Fig = figure('position',[2700 100 800 400])
plot(zscore(pc(1,:)/10),'-','color',[0,46,166]/255,'LineWidth',1.5)
hold on
plot(TPIfz,'color',[252,198,48]/255,'LineWidth',1.5)
plot(AMOfz,'color','k','LineWidth',1.5)
plot(zeros(1,81),'k','LineWidth',1)
plot(zeros(1,81)-0.5,'k--','LineWidth',1)
plot(zeros(1,81)+0.5,'k--','LineWidth',1)
r1 = corrcoef(pc1z,TPIfz);
text(5,-2.2,['R1 = ',num2str(roundn(r1(1,2),-2)),'*'],'fontsize',14)
r2 = corrcoef(pc1z,AMOfz);
text(55,-2.2,['R2 = ',num2str(roundn(r2(1,2),-2)),'**'],'fontsize',14)
scatter(num_pc1p,pc1z(num_pc1p),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
scatter(num_pc1n,pc1z(num_pc1n),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
set(gca,'XLim',[1,73]);
set(gca,'XTick',[7:10:73]);
set(gca,'XTickLabel',[1950:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
legend('PC1','IPO','AMO','Location','north','Orientation','horizontal')
legend('boxoff');
% print(Fig,['G:\figures\IAP\Yearly\20221124_IPO_SouthernOcean_2000m\TimeSeries_0_2000m_PC1.png'],'-dpng','-r300')
%% scatter plot
close all;
num_pc1p = find(pc1z>0);
num_pc1n = find(pc1z<0);
figure('position',[2700 100 500 500])
pc1p1 = find(pc1z >= 0.5 & pc1z < 0.9); pc1p2 = find(pc1z >= 0.9 & pc1z < 1.3); 
pc1p3 = find(pc1z >= 1.3 & pc1z < 1.7); pc1p4 = find(pc1z >= 1.7); 
numhb = [num_pc1n,num_pc1p];
ch = scatter(pc1z(numhb),TPIfz(numhb),[],pc1z(numhb),'filled')
% hold on
% scatter(AMOfz(find(pc1z>0 & pc1z<0.5)),TPIfz(find(pc1z>0 & pc1z<0.5)),'ro')
% scatter(AMOfz(find(pc1z>-0.5 & pc1z<0)),TPIfz(find(pc1z>-0.5 & pc1z<0)),'bo')
set(gca,'XLim',[-2,2]);
set(gca,'XTick',[-2:1:2]);
set(gca,'YLim',[-2,2]);
set(gca,'YTick',[-2:1:2]);
text(2.15,-0.2,'PC1','FontSize',12)
text(-0.4,2.25,'IPO','FontSize',12)
%
new_fig_handle = convert_to_std_coordinate_system(gca,0);
caxis([-2,2]);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
cb = colorbar;
set(cb,'Ticks',[-2:0.4:2],'TickLabels',[-2:0.4:2],'Position',[0.88,0.1,0.03,0.8],'FontSize',12);
set(cb.Label,'String','PC1','Rotation',270,'position',[4 0 0],'FontSize',12)
saveas(new_fig_handle,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\scatter_PC1_IPO.png'],'png')
%% IPO+
num_ipop = find(TPIfz>0.5);
num_ipon = find(TPIfz<-0.5);
% different depth
Temcom1 = nanmean(Tempadlf(:,:,:,num_ipop),4);
Temcom1d = nansum(Temcom1(:,:,1)+Temcom1(:,:,2:27).*permute(dweit,[3 2 1]),3)/700; % 补充 weight
close all;
Fig = diffdepth(Temcom1d(:,lats),1,lonData,latData,lats,1.2,0.2)
m_line(-60,[-60:1:-45],'linewidth',2,'color','k'); 
m_line(150,[-70:1:-45],'linewidth',2,'color','k'); 
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\IPOp_comTem0-700m.png'],'-dpng','-r300')
%% zonal mean Pac 150E-60W
Temcom1z = permute(nanmean(Temcom1(150:300,lats,:),1),[2 3 4 1]);
close all;
ticks = 0.1;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData,Temcom1z',[-1:0.02:1],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' IPO+ zonal mean Tem (South Pac)'])
print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\IPOp_zonalmean_Tem_150E_60W.png'],'-dpng','-r300')
%% zonal mean IO+Atl 60W-150E
lats = 1:45;
Temcom1_r = cat(1,Temcom1(300:360,:,:,:),Temcom1(1:299,:,:,:));
Temcom1z = permute(nanmean(Temcom1_r(1:210,lats,1:27,:),1),[2 3 4 1]);
% Temcom1z = permute(nanmean(Temcom1(150:300,lats,:),1),[2 3 4 1]);
close all;
ticks = 0.1;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:27),Temcom1z',[-1:0.02:1],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
set(gca,'XLim',[-75,-45]);
set(gca,'XTick',[-75:15:-45]);
set(gca,'XTickLabel',{'75^oS','60^oS','45^oS'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' IPO+ zonal mean Tem (South IO+Atl)'])
print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\IPOp_zonalmean_Tem_60W_150E.png'],'-dpng','-r300')
%% meridional mean
Temcom1m = latmean(Temcom1(:,:,1:41),[1:30],latData);
close all;
Fig = figure('position',[100 100 800 400]);
ticks = 0.1;
Temcom1m_r = cat(1,Temcom1m(150:300,:),Temcom1m(301:360,:),Temcom1m(1:149,:));
[c h] = contourf(lonData_r,-depthData,Temcom1m_r',[-1:0.02:1],'linestyle','none');
caxis([-ticks,ticks])
hold on
% line([220 220 340 340 220],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData(1:23),(nanmean(latmc(:,1:23,:),3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
set(gca,'XLim',[150,360+149]);
set(gca,'XTick',[150:60:360+149]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-2000,-800]);
set(gca,'YTick',[-2000:100:-800],'FontSize',14);
set(gca,'YTickLabel',[2000:-100:-800],'FontSize',14);
ylabel('Depth (m)')
title([' IPO+ meridional mean Tem 45-60S'])
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\IPOp_meridionalmean_Tem_45-60S.png'],'-dpng','-r300')
%% IPO+-  u v 
ucom1 = nanmean(uadf(:,:,num_ipon),3);
vcom1 = nanmean(vadf(:,:,num_ipon),3);
ucom1(:,181)=[]; vcom1(:,181)=[];
ucom1(:,1:10) = nan; ucom1(:,170:180)=nan;
close all;
ticks = 0.6;
Fig = figure('position',[100 100 800 400]);
% m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
    m_proj('stereographic','lon',0,'lat',-90,'radius',45,'rotangle',0);
hold on
[c h] = m_contourf(lonData,latData,ucom1(:,1:180)',[-0.9:0.03:0.9],'linestyle','none');
caxis([-ticks,ticks])
hold on
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks],'position',[0.78 0.09 0.03 0.85],'fontsize',12);
m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
d = 6; dd = 0.03;
[x,y] = meshgrid(lonData,latData);
umap = ucom1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,vcom1(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\IPOn_comUV.png'],'-dpng','-r300')






function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function contourfSPolar(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',45,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end
function contourVARra(var,val,ctr,lon1,lon2,lat1,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[lat1 lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
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
function [Fig] = diffdepth(Temcom,nlevel,lonData,latData,lats,val,ticks)
    Fig = figure('position',[100 100 600 300]);
%     contourVARra(Temcom(:,:,nlevel),[-val:ticks/10:val],1,0,360,-90,90,lonData,latData);
    contourfSPolar(Temcom(:,:,nlevel),[-val:ticks/10:val],ticks,lonData,latData(lats),12)
    caxis([-ticks,ticks]);
    ch = colorbar;
    set(ch,'Ticks',[-ticks:ticks/5:ticks],'position',[0.78 0.09 0.03 0.85]);
    % testdots(h05,'k',lonData,latData)
end
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end