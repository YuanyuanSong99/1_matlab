% trend
lats = 1:55;
dweit = depthData(2:27)-depthData(1:26); % depth weight
Tsubraw = permute(nansum(TemMMMsub(:,lats,1,:)+TemMMMsub(:,lats,2:27,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m raw
Tsubiv = permute(nansum(TemMMMsubd(:,lats,1,:)+TemMMMsubd(:,lats,2:27,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m internal variability
Tsubex = permute(nansum(TemLEMsub(:,lats,1,:)+TemLEMsub(:,lats,2:27,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m external forcing
startyr = 2005;  endyr = 2010;
str = 'raw';
var = permute(Tsubraw(:,:,startyr-1919:endyr-1919),[3 1 2]);
x = [1:size(var,1)]';
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1);
    end
end
h0 = trendtest(var,0.05); % t test trend
max(trd,[],'all')
min(trd,[],'all')
close all;
ftsz = 12; ticks = 0.6;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd*10,[-ticks*5:0.01:ticks*5],ticks,lonData,latData(lats),12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K decade^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots(h0,'k',lonData,latData(lats));
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\Tem0_700m_',num2str(startyr),'_',num2str(endyr),'trend_',str,'.png'],'-dpng','-r300')

%%
addpath D:\1_matlab\help;
addpath D:\1_matlab\help\seawater\;
% EOF
lats = 1:90;
dweit = depthData(2:37)-depthData(1:36); % depth weight
Tsub = permute(nansum(cat(3,Temsubadf(:,lats,1,:)*5,Temsubadf(:,lats,2:37,:).*(permute(dweit,[3 2 1]))),3)/700,[1 2 4 3]); % 0-700m
[eof_maps,pc,expvar]=eofnan(Tsub);  
EOF1z = -eof_maps(:,:,1)*std(pc(1,:)); 
PC1z = -pc(1,:)/std(pc(1,:)); 
% DI index (Pac_T - IO&Atl_T)
Tsub_long = permute(nansum(cat(3,Temsubadf_long(:,lats,1,:)*5,Temsubadf_long(:,lats,2:37,:).*(permute(dweit,[3 2 1]))),3)/700,[1 2 4 3]); % 0-700m
lats = [29:44];
Tsub_r = cat(1,Tsub_long(160:300,:,:),Tsub_long(301:360,:,:),Tsub_long(1:159,:,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
[spacz_long spac_long] = areamean(Tsub_r,1:141,lats,latData); 
[siaz_long sia_long] = areamean(Tsub_r,142:360,lats,latData); 
DI_rawlong = spac_long-sia_long;
DI_long = zscore(spac_long-sia_long);
DI = DI_long(5:end-4);
corrcoef(PC1z,DI)
%%
close all;
ftsz = 12; ticks = 0.1;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(EOF1z,[-0.4:0.01:0.4],ticks,lonData,latData(1:90),12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\Tem0-700m_EOFmap_0_90S_zscore.png'],'-dpng','-r300')
%% IPO index
% TPI and AMO index
[d1 d2 d3] = size(sstadf);
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
Pi = find(TPIfz > 0.5);    Ni = find(TPIfz < -0.5);
%%  PC1
pc1z = PC1z;
len = length(pc1z);
num_pc1p = find(pc1z>0.5);
num_pc1n = find(pc1z<-0.5);
close all;
Fig = figure('position',[2700 100 800 400])
plot(pc1z,'-','color',[0,46,166]/255,'LineWidth',1.5)
hold on
% plot(TPIfz,'color',[252,198,48]/255,'LineWidth',1.5)
plot(zeros(1,len),'k','LineWidth',1)
plot(zeros(1,len)-0.5,'k--','LineWidth',1)
plot(zeros(1,len)+0.5,'k--','LineWidth',1)
% rr = corrcoef(pc1z,TPIfz);
% text(5,-2.2,['Corr = ',num2str(roundn(rr(1,2),-2)),'**'],'fontsize',14)
scatter(num_pc1p,pc1z(num_pc1p),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
scatter(num_pc1n,pc1z(num_pc1n),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
set(gca,'XLim',[1,len]);
set(gca,'XTick',[7:10:len]);
set(gca,'XTickLabel',[1930:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
% print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230327_IPO_SouthernOcean\ensmean\Tem0-700m_EOFmap_PC1.png'],'-dpng','-r300')
%% PC1/DI/IPO reg meridional mean
Tem1m = latmean(Temsubadf(:,:,1:37,:),[29:44],latData);
[d1 d2 d3] = size(Tem1m);
Tem1mr = reshape(permute(Tem1m,[3 1 2]),[d3,d1*d2]);
%
index = TPIfz; pngname = 'IPO';
clear par14 h04
parfor i = 1:size(Tem1mr,2);
    [par14(i),h04(i),t] = reg1_ttest(index,Tem1mr(:,i),0.1,1);
end
par14r = reshape(par14,[d1 d2]);
max(par14r,[],'all')
h04r = reshape(h04,[d1 d2]);
map = nanmean(par14r,3);
%% mixed layer depth
lats = [29:44]; names = '46S-61S';
latm = latmean(HMXLadf,lats,latData);
HMXLr = reshape(permute(latm,[2 1]),[78,360]);
parfor i = 1:size(HMXLr,2);
    [paH(i),hH(i),t] = reg1_ttest(TPIfz(1:78),HMXLr(:,i),0.1,1);
end
%
latmMMM = latmean(HMXLrawLE,lats,latData);
HMXLcli = nanmean(latmMMM,2);
HMXLclir = -cat(1,HMXLcli(160:360,:),HMXLcli(1:159,:))/100; % cm -> m
%%
map_r = cat(1,map(160:300,:),map(301:360,:),map(1:159,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
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
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20240130_IPO_SO\',pngname,'.reg.Tem.meridionalmean46-61S.png'],'-dpng','-r300')
%% positive TPI com SST (each member)
sstcom = nanmean(sstadf(:,:,Pi),3);
Fig = figure('position',[100 100 700 700]);
subplot(2,1,1)
plot(TPIfz);
hold on
plot(zeros(96))
set(gca,'YLim',[-3,3]);
set(gca,'XLim',[1,len]);
set(gca,'XTick',[7:10:len]);
set(gca,'XTickLabel',[1930:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
subplot(2,1,2)
caxval = 0.25;
contourVARra(sstcom,[-6:caxval/10:6],caxval,0,360,-90,90,lonData,latData)
ch = colorbar;
set(ch,'Ticks',[-caxval:caxval/5:caxval]);
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230327_IPO_SouthernOcean\ensmean\IPO.com.po.sst.png'],'-dpng','-r300')
%% IPO/DI reg upper ocean Tem
[d1 d2 d3] = size(Tsub);
Tsubmr1 = reshape(permute(Tsub,[3 1 2]),[d3,d1*d2]);
index = TPIfz; pngname = 'IPO';
clear par11 h01
parfor i = 1:size(Tsubmr1,2);
    [par11(i),h01(i),t] = reg1_ttest(index,Tsubmr1(:,i),0.1,1);
end
par11r = reshape(par11,[d1,d2]);
h01r = reshape(h01,[d1,d2]);
%%
close all;
map = par11r; ftsz = 12;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 0.15; 
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(1:90),12)
testdots(h01r,[.4 .4 .4],lonData,latData(1:90));
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\',pngname,'.reg.Tem0_700.png'],'-dpng','-r300')
%% corr Pac with IO+Atl
% Tsubr1 = cat(1,Tsubm(151:301,:,:,:),Tsubm(302:360,:,:,:),Tsubm(1:150,:,:,:));
% lonData_r = [lonData(151:301);lonData(302:360);lonData(1:150)+360];
% %
% for s = 1:110
%     [spacz(:,s) spac] = areamean(Tsubr1(:,:,:,s),1:141,29:44,latData);
%     [siaz(:,s) sia] = areamean(Tsubr1(:,:,:,s),142:360,29:44,latData);
%     [r(s),p(s),n_eff2] = corr_eff(spacz,siaz,0.1);  %  90% confidence
% end
% mean(r)
%% IPO reg Taux Tauy
lats = 1:180;
[d1 d2 d3] = size(taux);
tauxr1 = reshape(permute(taux(:,lats,:),[3 1 2]),[d3,d1*d2]);
tauyr1 = reshape(permute(tauy(:,lats,:),[3 1 2]),[d3,d1*d2]);
clear par12 h02 par13 h03
parfor i = 1:size(tauxr1,2);
    [par12(i),h02(i),t] = reg1_ttest(TPIfz,tauxr1(:,i),0.1,1);
    [par13(i),h03(i),t] = reg1_ttest(TPIfz,tauyr1(:,i),0.1,1);
end
par12r = reshape(par12,[d1,d2]);
par13r = reshape(par13,[d1,d2]);
h02r = reshape(h02,[d1,d2]);
h03r = reshape(h03,[d1,d2]);
mapu = nanmean(par12r,3);
mapv = nanmean(par13r,3);
%%
close all;
Fig = figure('position',[100 100 800 400]);
ftsz = 12; caxval = 0.006;
% m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
[c h] = m_contourf(lonData,latData(lats),mapu',[-caxval*10:caxval/10:caxval*10],'linestyle','none');
caxis([-caxval,caxval])
hold on
load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-caxval:caxval/5:caxval],'position',[0.9 0.09 0.03 0.85],'fontsize',12);
m_coast('linewidth',1,'color','k');
m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
d = 8; dd = 0.0005;
[x,y] = meshgrid(lonData,latData(lats));
umap = mapu.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,mapv(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.0005 Pa','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.05./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
% print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230327_IPO_SouthernOcean\ensmean\IPO.reg.Taux&Tauy.png'],'-dpng','-r300')
%% zonal wind u
close all; ftsz = 12; ticks1 = 5;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par12r*10^3,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData(lats),12)
hold on
testdots(h02r,[.4 .4 .4],lonData,latData);
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'TAUX  (10^-^3 Pa)',4.5,12,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1]);
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\IPO.reg.Taux.png'],'-dpng','-r300')
%%
curlzr1 = permute(curlz,[3 2 1]); % 88*360*180
curlzr2 = reshape(curlzr1,d3,d1*d2);
clear par14 h04
parfor i = 1:size(curlzr2,2);
    [par14(i),h04(i),t] = reg1_ttest(TPIfz,curlzr2(:,i),0.1,1);
end
par14r = reshape(par14,[d1,d2]);
h04r = reshape(h04,[d1,d2]);
%%
close all; ftsz = 12;
map = mean(par14r,3)*10^9;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 8;
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData,12)
set_colorbar([0.83 0.08 0.03 0.88],'Wind Stress Curl  (10^-^9 Pa m^-^1)',4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks]);
hold on
testdots(h04r,[.4 .4 .4],lonData,latData);
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
d = 5; dd = 0.0006;
[x,y] = meshgrid(lonData,latData(lats));
umap = par12r.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par13r(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1,'MaxHeadSize',5) 
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\IPO.reg.curlz&Tau.png'],'-dpng','-r300')
%% IPO reg. slp
[d1 d2 d3] = size(slp);
slpr1 = reshape(permute(slp,[3 1 2]),[d3,d1*d2]);
clear par15 h05
parfor i = 1:size(slpr1,2);
    [par15(i),h05(i),t] = reg1_ttest(TPIfz,slpr1(:,i),0.1,1)
end
par15r = reshape(par15,[d1,d2]);
h05r = reshape(h05,[d1,d2]);
%%
close all;
map = mean(par15r,3);
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 60;
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(lats),12)
set_colorbar([0.83 0.08 0.03 0.88],'SLP  (Pa)',4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks]);
hold on
testdots(h05r,[.4 .4 .4],lonData,latData(lats));
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
d = 5; dd = 0.0006;
[x,y] = meshgrid(lonData,latData(lats));
umap = par12r.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par13r(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1,'MaxHeadSize',5) 
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\IPO.reg.slp&Tau.png'],'-dpng','-r300')
%% IPO reg. nhflx
nhflxadf = fsnsadf-shflxadf-lhflxadf-flnsadf;
lats = 1:120;
nhflxr1 = reshape(permute(nhflxadf(:,lats,:),[3 1 2]),[86,360*length(lats)]);
clear par16 h16 
parfor i = 1:size(nhflxr1,2);
    [par16(i),h16(i),t] = reg1_ttest(TPIfz,nhflxr1(:,i),0.1,1);
end
par16r = reshape(par16,[360,length(lats)]);
h16r = reshape(h16,[360,length(lats)]);
max(par16r,[],'all')
min(par16r,[],'all')
%%
close all;
map = mean(par16r,3); ftsz = 12;
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 6;
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(lats),12)
set_colorbar([0.83 0.08 0.03 0.88],'NHFLX (W m^-^2)',4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks]);
hold on
testdots(h16r,[.4 .4 .4],lonData,latData(lats));
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\IPO.reg.NHFLX.png'],'-dpng','-r300')
%% IPO reg. sst
lats = 1:120;
sstr1 = reshape(permute(sstadf(:,lats,:),[3 1 2]),[86,360*length(lats)]);
clear par17 h17 
parfor i = 1:size(sstr1,2);
    [par17(i),h17(i),t] = reg1_ttest(TPIfz,sstr1(:,i),0.1,1);
end
par17r = reshape(par17,[360,length(lats)]);
h17r = reshape(h17,[360,length(lats)]);
max(par17r,[],'all')
min(par17r,[],'all')
%
close all;
map = mean(par17r,3);
Fig = figure('position',[10 50 550 400]); ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
ticks = 0.2;
contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(lats),12)
set_colorbar([0.83 0.08 0.03 0.88],'SST (K)',4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks]);
hold on
testdots(h17r,[.4 .4 .4],lonData,latData(lats));
m_line([0:1:360],-46,'linewidth',2,'color','k');
m_line([0:1:360],-61,'linewidth',2,'color','k');
m_line(-60,[-61:1:-46],'linewidth',2,'color','k');
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\IPO.reg.SST.png'],'-dpng','-r300')
%% IPO reg zonal mean (SST)
reg_zonalmean(sstadf,TPIfz,[-0.1,0.15],[-0.1:0.05:0.15],'SST (K)','IPO.reg.zonalmean_SST',latData);
%% IPO reg zonal mean (net heat flux)
reg_zonalmean(nhflxadf,TPIfz,[-1.5,1.5],[-1.5:0.5:1.5],'Net Heat Flux (W m^-^2)','IPO.reg.zonalmean_NetHeatFlux',latData);
%% IPO reg zonal mean (zonal wind)
[paru1 paru2 hu1 hu2] = reg_zonalmean(taux,TPIfz,[-4,4]*10^-3,[-4:1:4]*10^-3,'Zonal Wind Stress (Pa)','IPO.reg.zonalmean_Taux',latData);
%% IPO reg zonal mean (curlz)
reg_zonalmean(permute(curlz,[2 1 3]),TPIfz,[-0.8,0.8]*10^-8,[-0.8:0.2:0.8]*10^-8,'Wind Stress Curl (Pa/m)','IPO.reg.zonalmean_curlz',latData);
%% IPO reg zonal mean (slp)
reg_zonalmean(slpadf,TPIfz,[-50,50],[-50:10:50],'Sea Level Pressure (Pa)','IPO.reg.zonalmean_SLP',latData);
%% PC1 corr other time series (scatter)
lats = [29:44];
var = taux; varname = 'TAUX'; titlename = 'Zonal Wind Stress';
var_r = cat(1,var(160:300,:,:),var(301:360,:,:),var(1:159,:,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
[tz1 t1] = areamean(var_r,1:141,lats,latData); % Pac
[tz2 t2] = areamean(var_r,142:360,lats,latData); % Atl & IO
%
tx1 = TPIfz; tx2 = zscore(t1-t2); pocname = 'IPO';
%
xlimv = [-3,3];
ylimv = [-4,3];
[r,p,n_eff2] = corr_eff(tx1,tx2,0.05)  %  90% confidence  
close all;
Fig = figure('position',[700 100 400 400]);
scatter(tx1,tx2,'*','Color',[0,46,166]/255,'LineWidth',1);
hold on
yb = polyfit(tx1,tx2,1);
plot(tx1,polyval(yb,tx1),'k','linewidth',3)
text(-2.5,-2.8,['Corr = ',num2str(roundn(r,-2)),'**'],'fontsize',12)
text(-2.5,-3.4,['Slope = ',num2str(roundn(yb(1),-2)),'**'],'fontsize',12)
set(gca,'XLim',xlimv);
set(gca,'FontSize',12);
set(gca,'YLim',ylimv);
xlabel(pocname,'FontSize',12),ylabel(titlename,'FontSize',12);
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\',pocname,'.regwith.',varname,'.scatter.png'],'-dpng','-r300')
%% MOHT-Ekman
filename1 = ['H:\CESM-post\PAC-pacemaker\yr_1x1\TAUX_192001-200512_ensmean.nc'];
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
TAUXMMM1 = ncread(filename1,'TAUX');
TAUXMMM1(:,:,end) = []; % remove the last year (only January)
filename2 = ['H:\CESM-post\PAC-pacemaker\yr_1x1\TAUX_200601-201312_ensmean.nc'];
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
TAUXMMM2 = ncread(filename2,'TAUX');
TAUXMMM2(:,:,end) = []; % remove the last year (only January)
TAUXMMM = cat(3,TAUXMMM1,TAUXMMM2)/10; % unit: dyn/cm2 -> N/m2
clear TAUXMMM1 TAUXMMM2
[MOHTall_p MOHT1_p MOHT2_p] = MOHT_Ekman(TAUXMMM,TemLEMsub,latData,depthData);
%
datadiru1='H:\CESM-post\LE\TAUX\historical\taux_yr_1x1\'; %指定批量数据所在的文件夹
filelistu1=dir([datadiru1,'*.nc']); %指定批量数据的类型
ncdisp([datadiru1,filelistu1(1).name]);
datadiru2='H:\CESM-post\LE\TAUX\rcp85\TAUX_yr_1x1\'; %指定批量数据所在的文件夹
filelistu2=dir([datadiru2,'*.nc']); %指定批量数据的类型
datadirt1='H:\CESM-post\LE\Temperature\his_sub_yr1x1\'; %指定批量数据所在的文件夹
filelistt1=dir([datadirt1,'*.nc']); %指定批量数据的类型
datadirt2='H:\CESM-post\LE\Temperature\rcp85_sub_yr1x1_2006-2013\'; %指定批量数据所在的文件夹
filelistt2=dir([datadirt2,'*.nc']); %指定批量数据的类型
for s=1:40
    s 
    filenameu1=[datadiru1,filelistu1(s).name];
    ncid=netcdf.open(filenameu1,'NC_NOWRITE');
    TAUXle1 = -ncread(filenameu1,'TAUX');
    TAUXle1(:,:,end) = [];
    filenameu2=[datadiru2,filelistu2(s).name];
    ncid=netcdf.open(filenameu2,'NC_NOWRITE');
    TAUXle2 = -ncread(filenameu2,'TAUX');
    TAUXle = cat(3,TAUXle1,TAUXle2);
    filenamet1=[datadirt1,filelistt1(s).name];
    ncid=netcdf.open(filenamet1,'NC_NOWRITE');
    Temle1 = ncread(filenamet1,'TEMP');
    Temle1(:,:,:,end) = [];    
    filenamet2=[datadirt2,filelistt2(s).name];
    ncid=netcdf.open(filenamet2,'NC_NOWRITE');
    Temle2 = ncread(filenamet2,'TEMP');
    Temle = cat(4,Temle1,Temle2(:,1:90,:,:));
    [MOHTall_le(:,:,s)  MOHT1_le(:,:,s) MOHT2_le(:,:,s)] = MOHT_Ekman(TAUXle,Temle,latData,depthData);
end
MOHTall_lem = mean(MOHTall_le,3);
MOHT1_lem = mean(MOHT1_le,3);
MOHT2_lem = mean(MOHT2_le,3);
%
close all;
Fig = figure('position',[100 100 500 500])
map = mean(MOHTall_lem,2)/10^15;
% map(85:96,:) = nan;
plot(latData(1:90),map,'linewidth',1.5)
set(gca,'XLim',[-90,90]);
set(gca,'XTick',[-80:20:80]);
ylabel('MOHT (PW)');xlabel('Latitude')
% print(Fig,['D:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\MOHT_climotology.png'],'-dpng','-r300')
%% MOHT anomaly,detrend, and 8-yr filter
MOHT1d = MOHT1_p - MOHT1_lem; % minus external forcing
MOHT2d = MOHT2_p - MOHT2_lem; % minus external forcing
MOHT1ad = (MOHT1d-mean(MOHT1d,2))'; % anomaly
MOHT2ad = (MOHT2d-mean(MOHT2d,2))'; % anomaly
dT = 1; cf = 1/8;
clear MOHT1adf MOHT2adf
for i = 1:90
    MOHT1adf(:,i) = lanczosfilter(MOHT1ad(:,i),dT,cf,[],'low'); % 8 year filtered
    MOHT2adf(:,i) = lanczosfilter(MOHT2ad(:,i),dT,cf,[],'low'); % 8 year filtered
end
MOHT1adf(1:4,:) = []; MOHT2adf(1:4,:) = [];
MOHT1adf(end-3:end,:) = []; MOHT2adf(end-3:end,:) = [];
%
index = TPIfz; pngname = 'IPO';
clear par7e h07e par8e h08e 
for i = 1:size(MOHT1adf,2);
    [par7e(i),h07e(i),t] = reg1_ttest(index,MOHT1adf(:,i),0.1,1);
    [par8e(i),h08e(i),t] = reg1_ttest(index,MOHT2adf(:,i),0.1,1);
end
%%
close all;
Fig = figure('position',[100 100 800 400]);
plot(latData(1:90),par7e/10^15,'r','linewidth',1.5)
hold on
plot(latData(1:90),par8e/10^15,'b','linewidth',1.5)
plot(latData(1:90),zeros(1,length(latData(1:90))),'k','linewidth',1.5)
plot(-61*ones(20),[-.015:.03/19:.015],'k--','linewidth',1.5);
plot(-46*ones(20),[-.015:.03/19:.015],'k--','linewidth',1.5);
scatter(latData(find(h07e == 1)),par7e(find(h07e == 1))/10^15,'*','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h08e == 1)),par8e(find(h08e == 1))/10^15,'*','MarkerEdgeColor','b','linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:10:-35]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',[-.015,.015],'YTick',[-.015:.005:.015]);
set(gca,'YTick',yticks,'FontSize',14);
ylabel('MOHT (PW)')
legend('Pac','Atl&IO','Location','north','Orientation','horizontal','position',[0.67,0.85,0.22,0.06])
legend('boxoff')
% print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\IPO.reg_zonalmean.MOHT_Ekman.png'],'-dpng','-r300')
%% MOHT vT 0-700m
% LE 
f = sw_f(latData(1:90));
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
londist = 2*pi*6.371393*10^6*cos(latData(1:90)*pi/180)/360; % distance between 2 longitudes 
dweit = depthData(2:37)-depthData(1:36);
levdist(1) = 4; levdist(2:37) = dweit; % z distance
datadirt1='H:\CESM-post\LE\Temperature\his_sub_yr1x1\'; %指定批量数据所在的文件夹
filelistt1=dir([datadirt1,'*.nc']); %指定批量数据的类型
datadirt2='H:\CESM-post\LE\Temperature\rcp85_sub_yr1x1_2006-2013\'; %指定批量数据所在的文件夹
filelistt2=dir([datadirt2,'*.nc']); %指定批量数据的类型
datadirv1='H:\CESM-post\LE\VVEL\his_vvel_yr1x1\'; %指定批量数据所在的文件夹
filelistv1=dir([datadirv1,'*.nc']); %指定批量数据的类型
ncdisp([datadirv1,filelistv1(1).name]);
datadirv2='H:\CESM-post\LE\VVEL\rcp85_vvel_yr1x1\'; %指定批量数据所在的文件夹
filelistv2=dir([datadirv2,'*.nc']); %指定批量数据的类型
k=length(filelistt1);
clear OHCz* MOHT7*
%
for s=1:40
    s 
    filenamet1=[datadirt1,filelistt1(s).name];
    ncid=netcdf.open(filenamet1,'NC_NOWRITE');
    Temle1 = ncread(filenamet1,'TEMP');
    Temle1(:,:,:,end) = [];    
    filenamet2=[datadirt2,filelistt2(s).name];
    ncid=netcdf.open(filenamet2,'NC_NOWRITE');
    Temle2 = ncread(filenamet2,'TEMP');
    Temle = cat(4,Temle1,Temle2(:,1:90,:,:));
    filenamev1=[datadirv1,filelistv1(s).name];
    ncid=netcdf.open(filenamev1,'NC_NOWRITE');
    Vvel1 = ncread(filenamev1,'VVEL')/100; % cm/s -> m/s
    Vvel1(:,:,:,end) = [];   
    filenamev2=[datadirv2,filelistv2(s).name];
    ncid=netcdf.open(filenamev2,'NC_NOWRITE');
    Vvel2 = ncread(filenamev2,'VVEL')/100; % cm/s -> m/s
    Vvel = cat(4,Vvel1,Vvel2);
   % 0-700m MOHT
    TemVel = Temle.*Vvel(:,1:90,:,:);
    OHCz = permute(ro*cp*sum((TemVel(:,1:90,1:37,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z 
    OHCzy = OHCz.*londist';
    OHCzy_r = cat(1,OHCzy(160:300,:,:,:),OHCzy(301:360,:,:,:),OHCzy(1:159,:,:,:));
    MOHT71_le(:,:,s) = permute(nansum(OHCzy_r(1:141,:,:,:),1),[2 3 4 1]); % Pac unit:W
    MOHT72_le(:,:,s) = permute(nansum(OHCzy_r(142:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
    MOHT7all_le(:,:,s) = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
end
clear OHC* TemVel Vvel filename* filelist* datadir*
MOHT71_lem = mean(MOHT71_le,3);
MOHT72_lem = mean(MOHT72_le,3);
%%
filename1 = ['H:\CESM-post\PAC-pacemaker\yr_1x1\VVEL_192001-200512_ensmean.nc'];
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
VVELMMM1 = ncread(filename1,'VVEL');
VVELMMM1(:,:,:,end) = []; % remove the last year (only January)
filename2 = ['H:\CESM-post\PAC-pacemaker\yr_1x1\VVEL_200601-201312_ensmean.nc'];
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
VVELMMM2 = ncread(filename2,'VVEL');
VVELMMM2(:,:,:,end) = []; % remove the last year (only January)
VVELMMM = cat(4,VVELMMM1,VVELMMM2);
clear VVELMMM1 VVELMMM2
Temp = TemLEMsub;
Vvel = VVELMMM(:,1:90,1:37,:)/100; % cm/s -> m/s
% 0-700m MOHT
TemVel = Temp(:,1:90,:,:).*Vvel(:,1:90,:,:);
OHCz = permute(ro*cp*sum((TemVel(:,1:90,1:37,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3]); % integrated along z
OHCzy = OHCz.*londist';
OHCzy_r = cat(1,OHCzy(160:300,:,:,:),OHCzy(301:360,:,:,:),OHCzy(1:159,:,:,:));
MOHT71_p = permute(nansum(OHCzy_r(1:141,:,:,:),1),[2 3 4 1]); % Pac unit:W
MOHT72_p = permute(nansum(OHCzy_r(142:end,:,:,:),1),[2 3 4 1]); % IO+Atl unit:W
MOHT7all_p = permute(nansum(OHCzy_r,1),[2 3 4 1]); % all
clear OHC* TemVel Vvel filename* filelist* datadir*
%% MOHT anomaly,detrend, and 8-yr filter
MOHT71d = MOHT71_p - MOHT71_lem(:,:); % minus external forcing
MOHT72d = MOHT72_p - MOHT72_lem(:,:); % minus external forcing
MOHT71ad = (MOHT71d-mean(MOHT71d,2))'; % anomaly
MOHT72ad = (MOHT72d-mean(MOHT72d,2))'; % anomaly
dT = 1; cf = 1/8;
clear MOHT71adf MOHT72adf
for i = 1:90
    MOHT71adf(:,i) = lanczosfilter(MOHT71ad(:,i),dT,cf,[],'low'); % 8 year filtered
    MOHT72adf(:,i) = lanczosfilter(MOHT72ad(:,i),dT,cf,[],'low'); % 8 year filtered
end
MOHT71adf(1:4,:) = []; MOHT72adf(1:4,:) = [];
MOHT71adf(end-3:end,:) = []; MOHT72adf(end-3:end,:) = [];
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
p1 = plot(latData(1:90),par7/10^15,'-','color',[0.90,0.51,0.31],'linewidth',1.5)
hold on
p2 = plot(latData(1:90),par8/10^15,'-','color',[0.17,0.73,0.69],'linewidth',1.5)
p3 = plot(latData(1:90),par7e/10^15,'r-','linewidth',1.5)
p4 = plot(latData(1:90),par8e/10^15,'b-','linewidth',1.5)
plot(latData(1:90),zeros(1,length(latData(1:90))),'k','linewidth',1.5)
plot(-61*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
plot(-46*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
scatter(latData(find(h07 == 1)),par7(find(h07 == 1))/10^15,'o','MarkerEdgeColor',[0.90,0.51,0.31],'linewidth',1)
scatter(latData(find(h08 == 1)),par8(find(h08 == 1))/10^15,'o','MarkerEdgeColor',[0.17,0.73,0.69],'linewidth',1)
scatter(latData(find(h07e == 1)),par7e(find(h07e == 1))/10^15,'o','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h08e == 1)),par8e(find(h08e == 1))/10^15,'o','MarkerEdgeColor','b','linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
set(gca,'YLim',[-.025,.032],'YTick',[-.03:.01:.03]);
set(gca,'YTick',yticks,'FontSize',14);
ylabel('MOHT (PW)')
% legend([p1,p2,p3,p4],'Pac_7_0_0_m','Atl&IO_7_0_0_m','Pac_E_k_m_a_n','Atl&IO_E_k_m_a_n','Location','north','Orientation','vertical','position',[0.14,0.6,0.18,0.3])
% legend('boxoff')
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\IPO.reg_zonalmean.MOHT_vt0_700m&Ekman.png'],'-dpng','-r300')
%% 1979-2012 trend
lats = 1:90;
dweit = depthData(2:37)-depthData(1:36); % depth weight
Tsubraw = permute(nansum(TemLEMsub(:,lats,1,:)*5+TemLEMsub(:,lats,2:37,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m
% Tsubraw = permute(nansum(TemMMMsubda(:,lats,1,:)*5+TemMMMsubda(:,lats,2:37,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); % 0-700m
%
Tano = Tsubraw-mean(Tsubraw(:,:,1981-1919:2010-1919),3); % minus climotology 1981-2010
startyr = 1979;
endyr = 2010;
var = permute(Tano(:,:,startyr-1919:endyr-1919),[3 1 2]);
x = [1:size(var,1)]';
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1);
    end
end
max(trd,[],'all')
min(trd,[],'all')
%%
close all;
ftsz = 12; ticks = 0.2;
map = trd;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(map*10,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(180,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K decade^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230327_IPO_SouthernOcean\ensmean\Tem_700m_trend_1979_2010_EX.png'],'-dpng','-r300')
%% (teleconnection) IPO reg Z200 zonal deviation
varz = permute(Z200adf,[3 1 2]);
index = TPIfz; pngname = 'IPO';
clear par1 h1 par2 h2 par3 h3 par4 h4
for i = 1:size(varz,2);
    for j = 1:size(varz,3);
        [par1(i,j),h1(i,j),t] = reg1_ttest(index(1:78),varz(:,i,j),0.05,1);
    end
end
par1(:,1:10) = nan; par1(:,170:180)=nan;
h1(1:2:end,1:2:end) = 0;
%%
close all;
ticks1 = 10; ftsz = 12;
Fig = figure('position',[100 100 650 300])
contourVARra(par1,[-14:1:14],ticks1,0,360,-90,90,lonData,latData)
testdots(h1,[.4 .4 .4],lonData,latData)
set(gca,'position',[0.08,0.11,0.775,0.815],'fontsize',ftsz);
set_colorbar([0.87 0.095 0.026 0.832],'Z200 (gpm)',3.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\IPO.reg.Z200.png'],'-dpng','-r300')





function contourVARra(var,val,ctr,lon1,lon2,lat1,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[lat1 lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
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
    Fig = figure('position',[100 100 600 300]);
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
function [par11 par12 h01 h02] = reg_zonalmean(nhflx,index,ylimval,yticks,titlename,pngname,latData)
lats = [1:90];
lonmh1 = permute(nanmean(nhflx(160:300,lats,:),1),[2 3 1]);
nhflx_r = cat(1,nhflx(300:360,:,:),nhflx(1:299,:,:));
lonmh2 = permute(nanmean(nhflx_r(1:210,lats,:),1),[2 3 1]);
% corr with PC1
var1 = lonmh1';
var2 = lonmh2';
clear par11 h01 par12 h02
for i = 1:size(var1,2);
    [par11(i),h01(i),t] = reg1_ttest(index,var1(:,i),0.1,1); % filtered
    [par12(i),h02(i),t] = reg1_ttest(index,var2(:,i),0.1,1); % filtered
end
max(par11,[],'all')
min(par11,[],'all')
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
set(gca,'XTick',[-70:10:-35]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',ylimval);
set(gca,'YTick',yticks,'FontSize',14);
ylabel(titlename)
legend('Pac','Atl&IO','Location','north','Orientation','horizontal','position',[0.67,0.85,0.22,0.06])
legend('boxoff')
print(Fig,['D:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\',pngname,'.png'],'-dpng','-r300')
end
function [MOHTall MOHT1 MOHT2] = MOHT_Ekman(tauxraw,TemLEMsub,latData,depthData)
    f = sw_f(latData);
    My = -tauxraw./f'; % Ekman 质量输运
    ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
    cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
    Vy = My/ro;
    Vy(:,85:96,:) = nan; % 赤道不适用ekman输送
    londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
    dweit = depthData(2:37) - depthData(1:36);
    levdist(1) = 1; levdist(2:37) = dweit; % z distance
    clear OHC OHCT MOHT
    OHC = permute(ro*cp*sum((TemLEMsub(:,1:90,1:37,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3])/700;
    OHCT = OHC(:,:,:).*Vy(:,1:90,:);
    OHCTw = OHCT.*londist(1:90)';
    MOHTall = permute(nansum(OHCTw,1),[2 3 1]); % zonal sum unit:W
    OHCTw_r = cat(1,OHCTw(160:300,:,:),OHCTw(301:360,:,:),OHCTw(1:159,:,:));
    MOHT1 = permute(nansum(OHCTw_r(1:141,:,:),1),[2 3 1]); % Pac unit:W
    MOHT2 = permute(nansum(OHCTw_r(142:end,:,:),1),[2 3 1]); % IO+Atl unit:W
end
function h0 = trendtest(var,alp)
% var is time*lon*lat
d1 = size(var,1);
d2 = size(var,2);
d3 = size(var,3);
h0 = zeros(d2,d3);
x = [1:d1]';
for i=1:d2;
    for j=1:d3;
        m3=var(:,i,j);
        par(i,j,:)=polyfit(x,m3,1); % regression parameters
        y1=polyval(permute(par(i,j,:),[3 2 1]),x);
        Q = sum((m3-y1).^2); % t test
        c = 1/sum((x-mean(x)).^2);
        T = par(i,j,1)/sqrt(c)/sqrt(Q/(d1-2));
        if abs(T) > tinv(1-alp/2,d1-2); % 双侧 t test
            h0(i,j) = 1;
        end
    end
end
end
function testdots(h,clor,lonData,latData)
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',2,'color',clor); % F-test dots
end