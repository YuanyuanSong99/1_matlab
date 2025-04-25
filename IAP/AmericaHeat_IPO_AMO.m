%% EOF
clc,clear,close all;
addpath E:\1_matlab\help;
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
%%
dweit = depthData(2:27)-depthData(1:26); % depth weight
Tsub = permute(nansum(Tempadlf(:,:,1,:)+Tempadlf(:,:,2:27,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m
[eof_maps,pc,expvar]=eofnan(Tsub); 
%%
close all;
ticks = 0.1;
Fig = figure('position',[100 100 600 300]);
contourVARra(eof_maps(:,:,1)*10,[-0.4:0.01:0.4],0.1,0,360,-90,90,lonData,latData)
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
m_line([0 360 ],[-45 -45],'linewidth',1.5,'color','k')
m_line([0 360 ],[45 45],'linewidth',1.5,'color','k')
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\Tem0-700m_EOFmap.png'],'-dpng','-r300')
%%  load u v 
filename = 'E:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\10mSI-MO\uwnd.10m.mon.mean.nc'; % 1836-2015
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
uwnd = ncread(filename,'uwnd');
uwnd = permute(mean(reshape(uwnd(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
ua = double(uwnd - mean(uwnd(:,:,42:71),3)); % remove climatology from 1981 to 2010
filename = 'E:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\10mSI-MO\vwnd.10m.mon.mean.nc'; % 1836-2015
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
%% u v climotology
ucom1 = nanmean(uwnd(:,:,:),3);
vcom1 = nanmean(vwnd(:,:,:),3);
ucom1(:,181)=[]; vcom1(:,181)=[];
ucom1(:,1:10) = nan; ucom1(:,170:180)=nan;
close all;
Fig = figure('position',[100 100 800 400]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
d = 6; dd = 0.6;
[x,y] = meshgrid(lonData,latData);
umap = ucom1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,vcom1(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
m_text(10,85,'10 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,10./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\uvwnd.climo.png'],'-dpng','-r300')
%%
pc1z = zscore(pc(1,:));
num_pc1p = find(pc1z>0.5);
num_pc1n = find(pc1z<-0.5);
num_ipop = find(TPIfz>0.5);
num_ipon = find(TPIfz<-0.5);
num_amop = find(AMOfz>0.5);
num_amon = find(AMOfz<-0.5);
num_ipop_amon = intersect(num_ipop,num_amon);
num_ipop_amop = intersect(num_ipop,num_amop);
num_ipon_amon = intersect(num_ipon,num_amon);
num_ipon_amop = intersect(num_ipon,num_amop);
%% PAT sst time series 30S-30N,140W-20W Pacific-Atlantic
[Pat_zs Pat] = areamean(Tsub,220:340,60:121,latData); 
close all;
Fig = figure('position',[2700 100 800 400])
plot(zscore(pc(1,:)/10),'-','color',[0,46,166]/255,'LineWidth',1.5)
hold on
% plot(Pat_zs,'--','color',[0,46,166]/255,'LineWidth',1.5)
plot(TPIfz,'color',[252,198,48]/255,'LineWidth',1.5)
plot(AMOfz,'color',[255,106,106]/255,'LineWidth',1.5)
plot(zeros(1,81),'k','LineWidth',1)
plot(zeros(1,81)-0.5,'k--','LineWidth',1)
plot(zeros(1,81)+0.5,'k--','LineWidth',1)
scatter(num_pc1p,pc1z(num_pc1p),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
scatter(num_pc1n,pc1z(num_pc1n),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
% scatter(num_ipop_amon,pc1z(num_ipop_amon),'o','MarkerEdgeColor',[255,106,106]/255,'LineWidth',1.5);
% scatter(num_ipon_amop,pc1z(num_ipon_amop),'o','MarkerEdgeColor',[255,106,106]/255,'LineWidth',1.5);
% scatter(num_ipop_amop,pc1z(num_ipop_amop),'o','MarkerEdgeColor','g','LineWidth',1.5);
% scatter(num_ipon_amon,pc1z(num_ipon_amon),'o','MarkerEdgeColor','g','LineWidth',1.5);
set(gca,'XLim',[1,73]);
set(gca,'XTick',[7:10:73]);
set(gca,'XTickLabel',[1950:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
legend('PC1','TPI','AMO','Location','north','Orientation','horizontal')
legend('boxoff');
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\Tem100-500m_EOFpc1.png'],'-dpng','-r300')
%% scatter plot
close all;
figure('position',[2700 100 500 500])
pc1p1 = find(pc1z >= 0.5 & pc1z < 0.9); pc1p2 = find(pc1z >= 0.9 & pc1z < 1.3); 
pc1p3 = find(pc1z >= 1.3 & pc1z < 1.7); pc1p4 = find(pc1z >= 1.7); 
numhb = [num_pc1n,num_pc1p];
ch = scatter(AMOfz(numhb),TPIfz(numhb),[],pc1z(numhb),'filled')
% hold on
% scatter(AMOfz(find(pc1z>0 & pc1z<0.5)),TPIfz(find(pc1z>0 & pc1z<0.5)),'ro')
% scatter(AMOfz(find(pc1z>-0.5 & pc1z<0)),TPIfz(find(pc1z>-0.5 & pc1z<0)),'bo')
set(gca,'XLim',[-2,2]);
set(gca,'XTick',[-2:1:2]);
set(gca,'YLim',[-2,2]);
set(gca,'YTick',[-2:1:2]);
text(2.15,-0.2,'AMO','FontSize',12)
text(-0.4,2.25,'IPO','FontSize',12)
%
new_fig_handle = convert_to_std_coordinate_system(gca,0);
caxis([-2,2]);
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
cb = colorbar;
set(cb,'Ticks',[-2:0.4:2],'TickLabels',[-2:0.4:2],'Position',[0.88,0.1,0.03,0.8],'FontSize',12);
set(cb.Label,'String','PC1','Rotation',270,'position',[4 0 0],'FontSize',12)
saveas(new_fig_handle,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\scatter_IPO_AMO_PC1.png'],'png')
%% PC1 .reg. zonal mean 0-180 east
lonm1 = permute(nanmean(Tempadlf(1:180,:,1:27,:),1),[2 3 4 1]);
lonm1m = permute(nanmean(Temp(1:180,:,1:27,:),1),[2 3 4 1]); % climotology
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
contourf(latData,-depthData(1:27),par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
line([-45 -45],[0 -700],'linewidth',2.5,'color','k')
line([45 45],[0 -700],'linewidth',2.5,'color','k')
[ac ah] = contour(latData,-depthData(1:27),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.PC1.zonal mean temperature (East)')
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\PC1.corr.zonalmean.0-180east.png'],'-dpng','-r300')
%% PC1 .reg. zonal mean 0-180 west
lonm1m = permute(nanmean(Temp(181:360,:,1:27,:),1),[2 3 4 1]); % climotology
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
contourf(latData,-depthData(1:27),par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
line([-45 -45],[0 -700],'linewidth',2.5,'color','k')
line([45 45],[0 -700],'linewidth',2.5,'color','k')
[ac ah] = contour(latData,-depthData(1:27),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title('Corr.PC1.zonal mean temperature (West)')
% print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\PC1.corr.zonalmean.180-360west.png'],'-dpng','-r300')
%% PC1 .reg. meridional mean (30S-30N) 
lats = [60:121]; names = '30S-30N'; index = pc(1,:);
latm = latmean(Tempadlf,lats,latData);
latmc = latmean(Temp,lats,latData); % climotology
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
Fig = figure('position',[100 100 800 400]);
contourf(lonData,-depthData,par11',[-1:0.1:1],'linestyle','none')
caxis([-1,1])
hold on
% line([220 220 340 340 220],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
[ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-1:0.2:1]);
set(gca,'XLim',[0,360]);
set(gca,'XTick',[0:60:360]);
set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title(['Corr.PC1.meridional mean temperature (',names,')'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\PC1.corr.meridionalmean',names,'.tem.png'],'-dpng','-r300')
%% PC1>0.5 u v 
ucom1 = nanmean(uadf(:,:,num_pc1p),3);
vcom1 = nanmean(vadf(:,:,num_pc1p),3);
ucom1(:,181)=[]; vcom1(:,181)=[];
ucom1(:,1:10) = nan; ucom1(:,170:180)=nan;
close all;
Fig = figure('position',[100 100 800 400]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
[c h] = m_contourf(lonData,latData,ucom1(:,1:180)',[-0.9:0.03:0.9],'linestyle','none');
caxis([-0.6,0.6])
hold on
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.6:0.12:0.6]);
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
d = 6; dd = 0.03;
[x,y] = meshgrid(lonData,latData);
umap = ucom1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,vcom1(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\uvwnd.pc1p.png'],'-dpng','-r300')
%% IPO+ & AMO- composition
% different depth
Temcom1 = nanmean(Tempadlf(:,:,:,num_ipop_amon),4);
Temcom1d = nansum(Temcom1(:,:,1)+Temcom1(:,:,2:27).*permute(dweit,[3 2 1]),3)/700; % 补充 weight
close all;
Fig = diffdepth(Temcom1d,1,lonData,latData,'0-700',1.2,0.2)
% m_line([220 340 340 220 220 ],[-30 -30 30 30 -30 ],'linewidth',1.5,'color','k')
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\diffdepth.0-700m.ipop_amon.tem.png'],'-dpng','-r300')
%% sst
Temcom1 = nanmean(Tempadlf(:,:,:,num_amon),4);
close all;
Fig = diffdepth(Temcom1(:,:,1),1,lonData,latData,'00',1.2,0.4)
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\sst.amon.png'],'-dpng','-r300')
%% zonal mean
Temcom1z = permute(nanmean(Temcom1(220:340,:,:),1),[2 3 4 1]);
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData,-depthData,Temcom1z',[-1:0.02:1],'linestyle','none')
caxis([-0.4,0.4])
hold on
line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.4:0.08:0.4]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' Com.zonal mean temperature 140W-20W'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\zonalmean.140W_20W.ipop_amon.png'],'-dpng','-r300')
% meridional mean
names = '45S-45N';
Temcom1m = latmean(Temcom1(:,:,1:27),[45:136],latData);
close all;
Fig = figure('position',[100 100 800 400]);
[c h] = contourf(lonData,-depthData(1:27),Temcom1m',[-1:0.02:1],'linestyle','none');
caxis([-0.4,0.4])
hold on
% line([220 220 340 340 220],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData(1:23),(nanmean(latmc(:,1:23,:),3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.4:0.08:0.4]);
set(gca,'XLim',[0,360]);
set(gca,'XTick',[0:60:360]);
set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' Com.meridional mean temperature (',names,')'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\meridionalmean.45S-45N.ipop_amon.png'],'-dpng','-r300')
% u v
ucom1 = nanmean(uadf(:,:,num_ipop_amon),3);
vcom1 = nanmean(vadf(:,:,num_ipop_amon),3);
ucom1(:,181)=[]; vcom1(:,181)=[];
ucom1(:,1:10) = nan; ucom1(:,170:180)=nan;
close all;
Fig = figure('position',[100 100 800 400]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
[c h] = m_contourf(lonData,latData,ucom1(:,1:180)',[-0.9:0.03:0.9],'linestyle','none');
caxis([-0.6,0.6])
hold on
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.6:0.12:0.6]);
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
d = 6; dd = 0.03;
[x,y] = meshgrid(lonData,latData);
umap = ucom1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,vcom1(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\uvwnd.ipop_amon.png'],'-dpng','-r300')
%% IPO+ & AMO+ composition
% different depth
Temcom1 = nanmean(Tempadlf(:,:,:,num_ipop_amop),4);
Temcom1d = nanmean(Temcom1(:,:,12:23),3);
close all;
Fig = diffdepth(Temcom1d,1,k,lonData,latData,'100-500',1.2,0.4)
m_line([220 340 340 220 220 ],[-30 -30 30 30 -30 ],'linewidth',1.5,'color','k')
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\diffdepth.100-500m.ipop_amop.tem.png'],'-dpng','-r300')
% zonal mean
Temcom1z = permute(nanmean(Temcom1(220:340,:,:),1),[2 3 4 1]);
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData,-depthData,Temcom1z',[-1:0.02:1],'linestyle','none')
caxis([-0.4,0.4])
hold on
line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.4:0.08:0.4]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' Com.zonal mean temperature 140W-20W'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\zonalmean.140W_20W.ipop_amop.png'],'-dpng','-r300')
% meridional mean
names = '30S-30N';
Temcom1m = latmean(Temcom1(:,:,1:23),[60:121],latData);
close all;
Fig = figure('position',[100 100 800 400]);
[c h] = contourf(lonData,-depthData(1:23),Temcom1m',[-1:0.02:1],'linestyle','none');
caxis([-0.4,0.4])
hold on
line([220 220 340 340 220],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData(1:23),(nanmean(latmc(:,1:23,:),3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.4:0.08:0.4]);
set(gca,'XLim',[0,360]);
set(gca,'XTick',[0:60:360]);
set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' Com.meridional mean temperature (',names,')'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\meridionalmean.30S-30N.ipop_amop.png'],'-dpng','-r300')
% u v
ucom1 = nanmean(uadf(:,:,num_ipop_amop),3);
vcom1 = nanmean(vadf(:,:,num_ipop_amop),3);
ucom1(:,181)=[]; vcom1(:,181)=[];
ucom1(:,1:10) = nan; ucom1(:,170:180)=nan;
close all;
Fig = figure('position',[100 100 800 400]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
[c h] = m_contourf(lonData,latData,ucom1(:,1:180)',[-1.8:0.03:1.8],'linestyle','none');
caxis([-0.6,0.6])
hold on
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.6:0.12:0.6]);
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
d = 6; dd = 0.03;
[x,y] = meshgrid(lonData,latData);
umap = ucom1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,vcom1(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\uvwnd.ipop_amop.png'],'-dpng','-r300')
%% IPO- & AMO- composition
% different depth
Temcom1 = nanmean(Tempadlf(:,:,:,num_ipon_amon),4);
Temcom1d = nanmean(Temcom1(:,:,12:23),3);
close all;
Fig = diffdepth(Temcom1d,1,k,lonData,latData,'100-500',1.2,0.4)
m_line([220 340 340 220 220 ],[-30 -30 30 30 -30 ],'linewidth',1.5,'color','k')
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\diffdepth.100-500m.ipon_amon.tem.png'],'-dpng','-r300')
% zonal mean
Temcom1z = permute(nanmean(Temcom1(220:340,:,:),1),[2 3 4 1]);
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData,-depthData,Temcom1z',[-1:0.02:1],'linestyle','none')
caxis([-0.4,0.4])
hold on
line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.4:0.08:0.4]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' Com.zonal mean temperature 140W-20W'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\zonalmean.140W_20W.ipon_amon.png'],'-dpng','-r300')
% meridional mean
names = '30S-30N';
Temcom1m = latmean(Temcom1(:,:,1:23),[60:121],latData);
close all;
Fig = figure('position',[100 100 800 400]);
[c h] = contourf(lonData,-depthData(1:23),Temcom1m',[-1:0.02:1],'linestyle','none');
caxis([-0.4,0.4])
hold on
line([220 220 340 340 220],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData(1:23),(nanmean(latmc(:,1:23,:),3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.4:0.08:0.4]);
set(gca,'XLim',[0,360]);
set(gca,'XTick',[0:60:360]);
set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' Com.meridional mean temperature (',names,')'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\meridionalmean.30S-30N.ipon_amon.png'],'-dpng','-r300')
% u v
ucom1 = nanmean(uadf(:,:,num_ipon_amon),3);
vcom1 = nanmean(vadf(:,:,num_ipon_amon),3);
ucom1(:,181)=[]; vcom1(:,181)=[];
ucom1(:,1:10) = nan; ucom1(:,170:180)=nan;
close all;
Fig = figure('position',[100 100 800 400]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
[c h] = m_contourf(lonData,latData,ucom1(:,1:180)',[-0.9:0.03:0.9],'linestyle','none');
caxis([-0.6,0.6])
hold on
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.6:0.12:0.6]);
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
d = 6; dd = 0.03;
[x,y] = meshgrid(lonData,latData);
umap = ucom1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,vcom1(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\uvwnd.ipon_amon.png'],'-dpng','-r300')
%% IPO- & AMO+ composition
% different depth
Temcom1 = nanmean(Tempadlf(:,:,:,num_ipon_amop),4);
Temcom1d = nansum(Temcom1(:,:,1)+Temcom1(:,:,2:27).*permute(dweit,[3 2 1]),3)/700; % 补充 weight
close all;
Fig = diffdepth(Temcom1d,1,lonData,latData,'100-500',1.2,0.2)
% m_line([220 340 340 220 220 ],[-30 -30 30 30 -30 ],'linewidth',1.5,'color','k')
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\diffdepth.0-700m.ipon_amop.tem.png'],'-dpng','-r300')
%% zonal mean
Temcom1z = permute(nanmean(Temcom1(220:340,:,:),1),[2 3 4 1]);
close all;
Fig = figure('position',[100 100 800 400]);
contourf(latData,-depthData,Temcom1z',[-1:0.02:1],'linestyle','none')
caxis([-0.4,0.4])
hold on
line([-30 -30 30 30 -30],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(latData,-depthData(1:23),(nanmean(lonm1m,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240); 
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.4:0.08:0.4]);
set(gca,'XLim',[-70,65]);
set(gca,'XTick',[-60:20:65]);
set(gca,'XTickLabel',{'60^oS','40^oS','20^oS','0^o','20^oN','40^oN','60^oN'},'FontSize',14);
set(gca,'YLim',[-500,0]);
set(gca,'YTick',[-500:100:0],'FontSize',14);
set(gca,'YTickLabel',[500:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' Com.zonal mean temperature 140W-20W'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\zonalmean.140W_20W.ipon_amop.png'],'-dpng','-r300')
%% meridional mean
names = '45S-45N';
Temcom1m = latmean(Temcom1(:,:,1:27),[45:136],latData);
close all;
Fig = figure('position',[100 100 800 400]);
[c h] = contourf(lonData,-depthData(1:27),Temcom1m',[-1:0.02:1],'linestyle','none');
caxis([-0.4,0.4])
hold on
% line([220 220 340 340 220],[-100 -500 -500 -100 -100],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData(1:23),(nanmean(latmc(:,1:23,:),3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.4:0.08:0.4]);
set(gca,'XLim',[0,360]);
set(gca,'XTick',[0:60:360]);
set(gca,'XTickLabel',{'0^o','60^oE','120^oE','180^o','120^oW','60^oW','0^o'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
title([' Com.meridional mean temperature (',names,')'])
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\meridionalmean.45S-45N.ipon_amop.png'],'-dpng','-r300')
%% u v
ucom1 = nanmean(uadf(:,:,num_ipon_amop),3);
vcom1 = nanmean(vadf(:,:,num_ipon_amop),3);
ucom1(:,181)=[]; vcom1(:,181)=[];
ucom1(:,1:10) = nan; ucom1(:,170:180)=nan;
close all;
Fig = figure('position',[100 100 800 400]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
[c h] = m_contourf(lonData,latData,ucom1(:,1:180)',[-0.9:0.03:0.9],'linestyle','none');
caxis([-0.6,0.6])
hold on
load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-0.6:0.12:0.6]);
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
d = 6; dd = 0.03;
[x,y] = meshgrid(lonData,latData);
umap = ucom1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,vcom1(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
print(Fig,['E:\figures\IAP\Yearly\220922EOF0_700m_IPOAMO\uvwnd.ipon_amop.png'],'-dpng','-r300')















function contourVARra(var,val,ctr,lon1,lon2,lat1,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[lat1 lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('E:/1_matlab/help/colorbar_mat/bl_re5.mat');
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
function [Fig] = diffdepth(Temcom,nlevel,lonData,latData,names,val,ticks)
    Fig = figure('position',[100 100 600 300]);
    contourVARra(Temcom(:,:,nlevel),[-val:ticks/10:val],1,0,360,-90,90,lonData,latData);
    caxis([-ticks,ticks]);
    ch = colorbar;
    set(ch,'Ticks',[-ticks:ticks/5:ticks]);
    % testdots(h05,'k',lonData,latData)
    title([' Com.T',' ',names,'m']);
end



