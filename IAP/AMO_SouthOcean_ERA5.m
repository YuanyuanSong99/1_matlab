% Trend 0-700m Tem150
clc,clear,close all;
addpath G:\1_matlab\help;
datadir='G:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% anomaly 
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
Tempa = double(Temp - nanmean(Temp(:,:,:,42:71),4)); % remove climatology from 1981-2010 
sst_iap = permute(Tempa(:,:,1,:),[1 2 4 3]);
%%
lats = 1:60;
% 0-700m
dweit = depthData(2:27)-depthData(1:26); % depth weight
Tsubraw = permute(nansum(Tempa(:,lats,1,:)+Tempa(:,lats,2:27,:).*(permute(dweit,[3 2 1])),3)/700,[1 2 4 3]); 
% 700-2000m
dweit = depthData(28:41)-depthData(27:40); % depth weight
Tsubraw = permute(nansum(Tempa(:,lats,41,:)+Tempa(:,lats,27:40,:).*(permute(dweit,[3 2 1])),3)/1300,[1 2 4 3]); 
clear Tsubraw_d Tsubraw_df
Tsubraw_d = detrend3d(Tsubraw);
Tsubraw_df = filter3d(Tsubraw_d,1,1/8);
%% trend
startyr = 1940;
endyr = 2020;
var = permute(Tsubraw(:,:,startyr-1939:endyr-1939),[3 1 2]);
x = [1:size(var,1)]';
clear trd
for i = 1:size(var,2);
    for j = 1:size(var,3);
        par=polyfit(x,var(:,i,j),1); % regression parameters
        trd(i,j) = par(1); 
    end
end
h0 = trendtest(var,0.05); % t test trend
max(trd,[],'all')
min(trd,[],'all')
%
close all;
ftsz = 12; ticks = 1.5/10;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd*10,[-ticks*5:0.01:ticks*5],ticks,lonData,latData(lats),12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K decade^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
testdots(h0,[.4 .4 .4],lonData,latData(lats));
% print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\Tem0-700m_trend_',num2str(startyr),'_',num2str(endyr),'.png'],'-dpng','-r300')
%%
close all
ftsz = 14;
Fig = figure('position',[610 50 1550 450]);
ax = axes('Position',[0.005 0.3 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(trd*10,[-ticks*5:0.01:ticks*5],ticks,lonData,latData(lats),12)
    hb1 = colorbar('location','southoutside');
    set(hb1,'Units','normalized','position',[0.05 0.15 0.88 0.05],'fontsize',ftsz,'Ticks',[-ticks:ticks/5:ticks],'TickLabels',[-ticks:ticks/5:ticks]);
    set(hb1.Label,'String','(K decade^-^1)','Units','normalized','fontsize',ftsz);
print(Fig,['G:\figures\IAP\Yearly\20230911_IPO_SouthernOcean_46_61S\Tem0-700m_trend_colorbar.png'],'-dpng','-r300')
%% raw non-detrended
lats = [29:44];
Tsub_r_0 = cat(1,Tsubraw(160:300,:,:),Tsubraw(301:360,:,:),Tsubraw(1:159,:,:));
[spacz_0 spac_0] = areamean(Tsub_r_0,1:141,lats,latData); 
[siaz_0 sia_0] = areamean(Tsub_r_0,142:360,lats,latData); 
DI_0 = zscore(spac_0-sia_0);
[sallz_0 sall_0] = areamean(Tsub_r_0,1:360,lats,latData); 
%%
close all;
Fig = figure('position',[700 100 800 400]);
plot(spac_0,'-','color','r','LineWidth',2)
hold on
plot(sia_0,'-','color','b','LineWidth',2)
plot(sall_0,'-','color','k','LineWidth',2)
plot(zeros(1,81),'k','LineWidth',1)
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[-.3,.3],'YTick',[-.3:.1:.3]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel('T_7_0_0  (K)','FontSize',14);
legend('Pacific','Atlantic-Indian Ocean','Southern Ocean','Location','northwest')
legend('boxoff')
title('Original')
% print(Fig,['G:\figures\IAP\Yearly\20230911_IPO_SouthernOcean_46_61S\time_non-detrended.png'],'-dpng','-r300')
%% raw detrended
lats = [29:44];
Tsub_rd_0 = cat(1,Tsubraw_d(160:300,:,:),Tsubraw_d(301:360,:,:),Tsubraw_d(1:159,:,:));
[spaczd_0 spacd_0] = areamean(Tsub_rd_0,1:141,lats,latData); 
[siazd_0 siad_0] = areamean(Tsub_rd_0,142:360,lats,latData); 
DIzd_0 = zscore(spacd_0-siad_0);
[sallzd_0 salld_0] = areamean(Tsub_rd_0,1:360,lats,latData); 
%%
close all;
Fig = figure('position',[700 100 800 400]);
plot(spacd_0,'-','color','r','LineWidth',2)
hold on
plot(siad_0,'-','color','b','LineWidth',2)
plot(salld_0,'-','color','k','LineWidth',2)
plot(zeros(1,81),'k','LineWidth',1)
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[-.2,.2],'YTick',[-.2:.1:.2]);
set(gca,'XGrid','on','YGrid','on');
xlabel('YEAR','FontSize',12),ylabel('T_7_0_0  (K)','FontSize',14);
legend('Pacific','Atlantic-Indian Ocean','Southern Ocean','Location','northwest')
legend('boxoff')
title('Detrended')
% print(Fig,['G:\figures\IAP\Yearly\20230911_IPO_SouthernOcean_46_61S\time_detrended.png'],'-dpng','-r300')


%%
lats = 1:55;
com1 = nanmean(Tsubraw(:,:,1979-1939:1993-1939),3);
com2 = nanmean(Tsubraw(:,:,1994-1939:2012-1939),3);
close all;
ftsz = 12; ticks = 0.25;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(com1,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_48_63S\Com_raw_1979_1993.png'],'-dpng','-r300')

Fig = figure('position',[610 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(com2,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k');
set_colorbar([0.83 0.08 0.03 0.88],'K',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['G:\figures\IAP\Yearly\20230321_IPO_SouthernOcean_48_63S\Com_raw_1994_2012.png'],'-dpng','-r300')

%% load Temperature
% clc,clear,close all;
addpath G:\1_matlab\help;
addpath G:\1_matlab\help\seawater\;
load('MatFile/Tempadlf');
Tempadlf_long = Tempadlf; % 1940-2020
Tempadlf(:,:,:,1:4) = [];
Tempadlf(:,:,:,end-3:end) = []; % 1944-2016
load("MatFile\TPI_filtered8_hadley.mat");
TPIf_long = TPI_filtered8_hadley(51:151); % 1920-2020
TPIf_raw = TPI_filtered8_hadley(71:151); % 1940-2020
TPIfz = zscore(TPIf_raw(5:end-4)); % 1944-2016
load("MatFile\AMOf8_hadley.mat");
AMOf_long = AMOf8_hadley(51:151); % 1920-2020
AMOf_raw = AMOf8_hadley(71:151); % 1940-2020
AMOfz = zscore(AMOf_raw(5:end-4)); % 1944-2016
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
[nlon nlat nlev nyr] = size(Tempadlf);
%%  load u v ；taux tauy curlz
datadir='G:\data\ERA5_1x1\10m_wind\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE');
lonDataE = ncread(filetemp,'longitude'); % 0-360
lonEmap = lonDataE; lonEmap(361) = lonEmap(1);
latDataE = ncread(filetemp,'latitude'); % 90 - -90
%
k=length(filelist);
for s=1:k-2;  % 1940-2020
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    u10 = ncread(filename1,'u10');     uwnd(:,:,s) = mean(u10,3);   % yearly mean
    v10 = ncread(filename1,'v10');     vwnd(:,:,s) = mean(v10,3);   % yearly mean
    netcdf.close(ncid);  %关闭文件    
    [taux_raw(:,:,s) tauy_raw(:,:,s)]= ra_windstr(uwnd(:,:,s),vwnd(:,:,s)); % lon*lat
    curlz0(:,:,s) = ra_windstrcurl(latDataE,lonDataE,uwnd(:,:,s)',vwnd(:,:,s)',0); % lat*lon
end
curlz_raw = permute(curlz0,[2 1 3]);
[uadf] = caladf(uwnd);
uadf(:,:,1:4) = []; uadf(:,:,end-3:end) = [];
[vadf] = caladf(vwnd);
vadf(:,:,1:4) = []; vadf(:,:,end-3:end) = [];
[taux] = caladf(taux_raw);
taux(:,:,1:4) = []; taux(:,:,end-3:end) = [];
[tauy] = caladf(tauy_raw);
tauy(:,:,1:4) = []; tauy(:,:,end-3:end) = [];
[curlz] = caladf(curlz_raw);
curlz(:,:,1:4) = []; curlz(:,:,end-3:end) = [];
%%  load slp
datadir='G:\data\ERA5_1x1\sea_level_pressure\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE');
k=length(filelist);
for s=1:k-2;  % 1940-2020
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    slp0 = ncread(filename1,'msl');     slp(:,:,s) = mean(slp0,3);   % yearly mean
    netcdf.close(ncid);  %关闭文件    
end
[slpadf] = caladf(slp);
slpadf(:,:,1:4) = []; slpadf(:,:,end-3:end) = [];
% load sst
filename = 'G:\data\hadley\HadISST_sst.nc'; % 1870-2021
ncid=netcdf.open(filename,'NOWRITE');
ncdisp(filename);
sst = ncread(filename,'sst');
sst(find(sst == -1000)) = nan;
lonDataH = ncread(filename,'longitude'); % -179.5-179.5
latDataH = ncread(filename,'latitude');  % 89.5- -89.5 
sst = permute(mean(reshape(sst(:,:,841:1812),[360 180 12 972/12]),3),[1 2 4 3]); % yearly 1940-2020
ssta = double(sst - mean(sst(:,:,42:71),3)); % remove climatology from 1981 to 2010
% 转换为相同lonData latData
ssta1 = fliplr(ssta);
ssta2 = cat(1,ssta1(182:end,:,:,:),ssta1(1:181,:,:,:));
sstar = permute(ssta2,[3 1 2]);
dT = 1; cf = 1/8;
clear sstad sstadf
for i = 1:360
    for j = 1:180
    sstad(:,i,j) = detrend(sstar(:,i,j)); % linear detrend
    sstadf(:,i,j) = lanczosfilter(sstad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
sstadf = permute(sstadf,[2 3 1]);
sstadf(:,:,1:4) = []; sstadf(:,:,end-3:end) = [];
% load surface heat data
datadir='G:\data\ERA5_1x1\heat\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
k=length(filelist);
for s=1:k-2;  % 1940-2020
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    slhf = ncread(filename1,'slhf');     lhf(:,:,s) = mean(slhf,3);   % latent heat flux  upward  J m-2 
    ssr = ncread(filename1,'ssr');       swf(:,:,s) = mean(ssr,3);    % short wave flux   downward
    str = ncread(filename1,'str');       lwf(:,:,s) = mean(str,3);    % long wave flux    upward
    sshf = ncread(filename1,'sshf');     shf(:,:,s) = mean(sshf,3);   % sensible heat flux upward    
    netcdf.close(ncid);  %关闭文件    
end
%% climo
nhflx = (swf+lhf+lwf+shf)/3600; % positive downward W m-2
[nhflxadf] = caladf(nhflx);
nhflxadf(:,:,1:4) = []; nhflxadf(:,:,end-3:end) = [];
mapcli = nanmean(nhflx,3);
close all;
Fig = figure('position',[10 50 550 400]);
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
m_contourf(lonDataE,latDataE,mapcli',[-400:20:400],'linestyle','none');
caxis([-100,100]);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
colorbar
m_coast('linewidth',1,'color','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
% print(Fig,['G:\figures\IAP\Yearly\nhflx_ERA5.png'],'-dpng','-r300')
%% load z200
datadir='G:\data\ERA5_1x1\Geopotential\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE');
level = ncread(filetemp,'level');
k=length(filelist);
for s=1:k-2;  % 1940-2020
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    z = ncread(filename1,'z');       z200(:,:,s) = mean(z(:,:,3,:),4);   % yearly mean
    netcdf.close(ncid);  %关闭文件    
end
z200adf = caladf(z200)/9.8; % gpm
z200adf(:,:,1:4) = []; z200adf(:,:,end-3:end) = [];
%% EOF
clear curlz0 u10 v10 slhf slp0 sshf ssr z
lats = 1:90;
dweit = depthData(28:41)-depthData(27:40); % depth weight
Tsub_long = permute(nansum(Tempadlf_long(:,lats,41,:)+Tempadlf_long(:,lats,27:40,:).*(permute(dweit,[3 2 1])),3)/1300,[1 2 4 3]); % 700-2000m, 1940-2020
Tsub = permute(nansum(Tempadlf(:,lats,41,:)+Tempadlf(:,lats,27:40,:).*(permute(dweit,[3 2 1])),3)/1300,[1 2 4 3]); % 700-2000m
[eof_maps,pc,expvar]=eofnan(Tsub);  
EOF1z = eof_maps(:,:,1)*std(pc(1,:)); 
PC1z = pc(1,:)/std(pc(1,:)); 
%%
close all;
ftsz = 12; ticks = 0.1;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(EOF1z,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
% m_line([0:1:360],-46,'linewidth',2,'color','k'); 
% m_line([0:1:360],-61,'linewidth',2,'color','k'); 
% m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
% m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['G:\figures\IAP\Yearly\20231008_AMO_SOdeep\Tem700-2000m_EOFmap_0_90S_zscore.png'],'-dpng','-r300')
%% DI index
lats = [29:44];
Tsub_r = cat(1,Tsub(160:300,:,:),Tsub(301:360,:,:),Tsub(1:159,:,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
[spacz spac] = areamean(Tsub_r,1:141,lats,latData); 
[siaz sia] = areamean(Tsub_r,142:360,lats,latData); 
DI = zscore(spac-sia);
% long 1940-2020
Tsub_longr = cat(1,Tsub_long(160:300,:,:),Tsub_long(301:360,:,:),Tsub_long(1:159,:,:));
[spacz_long spac_long] = areamean(Tsub_longr,1:141,lats,latData); 
[siaz_long sia_long] = areamean(Tsub_longr,142:360,lats,latData); 
DI_rawlong = spac_long-sia_long;
DI_long = zscore(spac_long-sia_long);
%%  time series PC1 & AMO
close all;
Fig = figure('position',[2700 100 800 400])
plot([5:77],PC1z,'color','K','LineWidth',1.5)
hold on
plot([5:77],AMOfz,'color',[.8 .05 .05],'LineWidth',1.5)
plot(zeros(1,81),'k','LineWidth',1)
set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
legend('PC1','AMO','Location','north','Orientation','horizontal')
legend('boxoff');
print(Fig,['G:\figures\IAP\Yearly\20231008_AMO_SOdeep\TimeSeries_PC1_AMO.png'],'-dpng','-r300')
%%  time series DI & TPI & AMO
TPIfz_long = zscore(TPIf_long);
AMOfz_long = zscore(AMOf_long);
close all;
Fig = figure('position',[2700 100 800 400])
plot([5+20:77+20],DI,'-','color','k','LineWidth',1.5)
hold on
plot([1:97],TPIfz_long(1:97),'color',[.8 .05 .05],'LineWidth',1.5)
plot([1:97],AMOfz_long(1:97),'color',[.13 .27 .86],'LineWidth',1.5)
plot(zeros(1,101),'k','LineWidth',1)
plot([1+20:4+20],DI_long(1:4),'--','color','k','LineWidth',1.5)
plot([78+20:81+20],DI_long(78:81),'--','color','k','LineWidth',1.5)
plot([97:101],TPIfz_long(97:101),'--','color',[.8 .05 .05],'LineWidth',1.5)
plot([97:101],AMOfz_long(97:101),'--','color',[.13 .27 .86],'LineWidth',1.5)
[r1,p,n_eff2] = corr_eff(DI,TPIfz,.1)
text(5,-2.5,['R_P_m_A_I_ _&_ _I_P_O = ',num2str(roundn(r1,-2)),'**'],'fontsize',14,'color',[.8 .05 .05])
[r2,p,n_eff2] = corr_eff(DI,AMOfz,.1)
text(60,-2.5,['R_P_m_A_I_ _&_ _A_M_O = ',num2str(roundn(r2,-2))],'fontsize',14,'color',[.13 .27 .86])
% scatter(num_pc1p,PC1z(num_pc1p),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
% scatter(num_pc1n,PC1z(num_pc1n),'*','MarkerEdgeColor',[0,46,166]/255,'LineWidth',1.5);
set(gca,'XLim',[1,101]);
set(gca,'XTick',[1:10:101]);
set(gca,'XTickLabel',[1920:10:2020],'FontSize',14);
set(gca,'YLim',[-3,3]);
set(gca,'YTick',[-3:1:3]);
xlabel('YEAR','FontSize',12),ylabel('Normalized index','FontSize',14);
legend('PmAI','IPO','AMO','Location','north','Orientation','horizontal')
legend('boxoff');
print(Fig,['G:\figures\IAP\Yearly\20230911_IPO_SouthernOcean_46_61S\TimeSeries_DI_IPO_AMO.png'],'-dpng','-r300')
%% time series -- 
t1 = DI_long; t2 = PC1z; str1 = 'PmAI'; str2 = 'PC1'; pngname = 'DI_PC1';
[r,p,n_eff2] = corr_eff(t1(5:end-4),t2,0.05)  %  90% confidence

% t1 = spac_long; t2 = sia_long; str1 = 'Pacific'; str2 = 'Atlantic & Indian Ocean'; pngname = 'Pac_IO+Atl';
% [r,p,n_eff2] = corr_eff(t1(5:77),t2(5:77),0.1)  %  90% confidence
% corr_eff(spacd_0,siad_0,0.1)

close all;
Fig = figure('position',[700 100 800 400]);
plot([5:77],t1(5:end-4),'-','color','r','LineWidth',1.5)
hold on
plot(DIzd_0,'r--','LineWidth',1)
plot([5:77],t2(1:end),'color','b','LineWidth',1.5)
% plot([5:77],t2(5:end-4),'color','b','LineWidth',1.5)

plot([1:4],t1(1:4),'--','color','r','LineWidth',1.5)
plot([78:81],t1(end-3:end),'--','color','r','LineWidth',1.5)
plot(zeros(1,81),'k','LineWidth',1)

% plot([1:4],t2(1:4),'--','color','b','LineWidth',1.5)
% plot([78:81],t2(end-3:end),'--','color','b','LineWidth',1.5)
% plot(spacd_0,'r--','LineWidth',1)
% plot(siad_0,'b--','LineWidth',1)
% text(5,-0.17,['Corr_a_n_n_u_a_l = ','-0.17'],'fontsize',14)
% text(55,-0.17,['Corr_d_e_c_a_d_a_l = ',num2str(roundn(r,-2)),'*'],'fontsize',14)

set(gca,'XLim',[1,81]);
set(gca,'XTick',[1:10:81]);
set(gca,'XTickLabel',[1940:10:2020],'FontSize',14);
set(gca,'YLim',[-3,4]);
set(gca,'YTick',[-3:1:4]);
xlabel('YEAR','FontSize',12),ylabel('T_7_0_0  (K)','FontSize',14);
legend(str1,'PmAI_a_n_n_u_a_l',str2,'Location','north','Orientation','horizontal')
% legend(str1,str2,'Location','north','Orientation','horizontal')
legend('boxoff');
print(Fig,['G:\figures\IAP\Yearly\20230911_IPO_SouthernOcean_46_61S\TimeSeries_',pngname,'.png'],'-dpng','-r300')
%% DI&IPO scatter plot
close all;
num_pc1p = find(DI>0);
num_pc1n = find(DI<0);
figure('position',[2700 100 500 500])
pc1p1 = find(DI >= 0.5 & DI < 0.9); pc1p2 = find(DI >= 0.9 & DI < 1.3); 
pc1p3 = find(DI >= 1.3 & DI < 1.7); pc1p4 = find(DI >= 1.7); 
numhb = [num_pc1n;num_pc1p];
ch = scatter(DI(numhb),TPIfz(numhb),[],DI(numhb),'filled')
hold on
yb = polyfit(DI,TPIfz,1);
plot(DI,polyval(yb,DI),'k','linewidth',2)
% scatter(AMOfz(find(DI>0 & DI<0.5)),TPIfz(find(DI>0 & DI<0.5)),'ro')
% scatter(AMOfz(find(DI>-0.5 & DI<0)),TPIfz(find(DI>-0.5 & PC1z<0)),'bo')
set(gca,'XLim',[-2,2]);
set(gca,'XTick',[-2:1:2]);
set(gca,'YLim',[-2,2]);
set(gca,'YTick',[-2:1:2]);
text(1,-2,['Slope = ',num2str(roundn(yb(1),-2)),'**'],'fontsize',12)
text(2.15,-0.2,'PmAI','FontSize',12)
text(-0.4,2.25,'IPO','FontSize',12)
%
new_fig_handle = convert_to_std_coordinate_system(gca,0);
caxis([-2,2]);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
cb = colorbar;
set(cb,'Ticks',[-2:0.4:2],'TickLabels',[-2:0.4:2],'Position',[0.88,0.1,0.03,0.8],'FontSize',12);
set(cb.Label,'String','PmAI','Rotation',270,'position',[4 0 0],'FontSize',12)
saveas(new_fig_handle,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\scatter_DI_IPO.png'],'png')
%% DI time series with CESM PAC pacemaker
itIAP = DI_rawlong; itLE = DIf_LE; itensm = DIensm; itpac = DIpac; str1 = 'PmAI'; str2 ='IAP'; yval = .3;
close all;
Fig = figure('position',[2700 100 800 400])
maxval = max(itpac,[],2);
minval = min(itpac,[],2);
xval = [1:94];
hf = fill([xval,fliplr(xval)],[minval' fliplr(maxval')],[0.87,0.68,0.68]);
set(hf,'edgecolor',[0.87,0.68,0.68])
hold on
plot(zeros(1,101),'k','LineWidth',1)
% for s = 1:20
%     plot([5:90],itpac(5:90,s),'-','LineWidth',0.5)
%     plot([1:4],itpac(1:4,s),'-','LineWidth',0.5)
%     plot([91:94],itpac(91:94,s),'-','LineWidth',0.5)
% end
p1 = plot([25:97],itIAP(5:end-4),'k-','LineWidth',3);
p2 = plot([5:90],itensm(5:90),'-','color',[.8 .05 .05],'LineWidth',3);
p3 = plot([5:90],itLE(5:90),'-','color','b','LineWidth',1.5);
plot([21:25],itIAP(1:5),'k:','LineWidth',3)
plot([97:101],itIAP(end-4:end),'k:','LineWidth',3)
plot([1:5],itensm(1:5),':','color',[.8 .05 .05],'LineWidth',3)
plot([90:94],itensm(end-4:end),':','color',[.8 .05 .05],'LineWidth',3)
plot([1:5],itLE(1:5),'--','color','b','LineWidth',1.5)
plot([90:94],itLE(end-4:end),'--','color','b','LineWidth',1.5)
plot(44*ones(20,1),[-.3:.6/19:.3],'color',[0.18,0.60,0.49],'LineWidth',1.5); % 1963, Agung volcanic
plot(63*ones(20,1),[-.3:.6/19:.3],'color',[0.18,0.60,0.49],'LineWidth',1.5); % 1982, El Chichon volcanic
plot(72*ones(20,1),[-.3:.6/19:.3],'color',[0.18,0.60,0.49],'LineWidth',1.5); % 1991, Pinatubo volcanic
text(46,-0.28,'Agung','color',[0.18,0.60,0.49],'fontsize',14,'rotation',90)
text(65,-0.28,'El Chichon','color',[0.18,0.60,0.49],'fontsize',14,'rotation',90)
text(74,-0.28,'Pinatubo','color',[0.18,0.60,0.49],'fontsize',14,'rotation',90)
set(gca,'XLim',[1,101]);
set(gca,'XTick',[1:10:101]); 
set(gca,'XTickLabel',[1920:10:2020],'FontSize',14);
set(gca,'YLim',[-.3,.3]);
set(gca,'YTick',[-.3:.1:.3]);
xlabel('YEAR','FontSize',12),ylabel([str1, ' (K)'],'FontSize',14);
legend([p1,p2,p3],str2,['POGA-LE'],'LE','Location','northeast','Orientation','vertical')
legend('boxoff');
print(Fig,['G:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\Times_IAP_Pacemaker_LE_',str1,'.png'],'-dpng','-r300')
%% IPO time series  with CESM PAC pacemaker
itIAP = TPIf_long; itLE = TPIf_LE; itensm = TPIensm; itpac = TPIpac; str1 = 'IPO'; str2 = 'Hadisst'; yval = 1;
close all;
Fig = figure('position',[2700 100 800 400])
maxval = max(itpac,[],2);
minval = min(itpac,[],2);
xval = [1:94];
hf = fill([xval,fliplr(xval)],[minval' fliplr(maxval')],[0.87,0.68,0.68]);
set(hf,'edgecolor',[0.87,0.68,0.68])
hold on
plot(zeros(1,101),'k','LineWidth',1)
% for s = 1:20
%     plot([5:90],itpac(5:90,s),'-','LineWidth',0.5)
%     plot([1:4],itpac(1:4,s),'-','LineWidth',0.5)
%     plot([91:94],itpac(91:94,s),'-','LineWidth',0.5)
% end
p1 = plot([1:97],itIAP(1:end-4),'k-','LineWidth',3);
p2 = plot([5:90],itensm(5:90),'-','color',[.8 .05 .05],'LineWidth',3);
p3 = plot([5:90],itLE(5:90),'-','color','b','LineWidth',1.5);
plot([97:101],itIAP(end-4:end),'k:','LineWidth',3)
plot([1:5],itensm(1:5),':','color',[.8 .05 .05],'LineWidth',3)
plot([90:94],itensm(end-4:end),':','color',[.8 .05 .05],'LineWidth',3)
plot([1:5],itLE(1:5),'--','color','b','LineWidth',1.5)
plot([90:94],itLE(end-4:end),'--','color','b','LineWidth',1.5)
plot(44*ones(20,1),[-1:2/19:1],'color',[0.18,0.60,0.49],'LineWidth',1.5); % 1963, Agung volcanic
plot(63*ones(20,1),[-1:2/19:1],'color',[0.18,0.60,0.49],'LineWidth',1.5); % 1982, El Chichon volcanic
plot(72*ones(20,1),[-1:2/19:1],'color',[0.18,0.60,0.49],'LineWidth',1.5); % 1991, Pinatubo volcanic
text(46,-0.9,'Agung','color',[0.18,0.60,0.49],'fontsize',14,'rotation',90)
text(65,-0.9,'El Chichon','color',[0.18,0.60,0.49],'fontsize',14,'rotation',90)
text(74,-0.9,'Pinatubo','color',[0.18,0.60,0.49],'fontsize',14,'rotation',90)
set(gca,'XLim',[1,101]);
set(gca,'XTick',[1:10:101]); 
set(gca,'XTickLabel',[1920:10:2020],'FontSize',14);
set(gca,'YLim',[-yval,yval]);
set(gca,'YTick',[-yval:yval/2:yval]);
xlabel('YEAR','FontSize',12),ylabel([str1, ' (K)'],'FontSize',14);
legend([p1,p2,p3],str2,['POGA-LE'],'LE','Location','northeast','Orientation','vertical')
legend('boxoff');
print(Fig,['G:\figures\CESM\Yearly\PAC-pacemaker\20230606_IPO_SouthernOcean_46_61S\ensmean\Times_IAP_Pacemaker_LE_',str1,'.png'],'-dpng','-r300')
%% PC1 .reg. zonal mean 160E-60W
lonm1 = permute(nanmean(Tempadlf(160:300,lats,1:27,:),1),[2 3 4 1]);
% lonm1m = permute(nanmean(Temp(160:300,lats,1:27,:),1),[2 3 4 1]); % climotology
% corr with PC1
var = permute(lonm1,[3 1 2]);
clear par11 h05 
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = reg1_ttest(pc(1,:),var(:,i,j),0.05,1); % filtered
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
title('Corr.PC1.zonal mean temperature (South Pac)')
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\Tem0-700m_EOFmap_0_90S_Pac.png'],'-dpng','-r300')
%% PC1/DI/IPO .reg. meridional mean 
lats = [29:44]; names = '46S-61S'; index = DI; pngname = 'DI';
latm = latmean(Tempadlf,lats,latData);
% latmc = latmean(Temp,lats,latData); % climotology
% corr with 
var = permute(latm,[3 1 2]);
clear par11 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i,j),0.1,1);
    end
end
max(par11,[],'all')
min(par11,[],'all')
%%
close all;
ticks = 0.1; ftsz=12;
par11_r = cat(1,par11(160:300,:),par11(301:360,:),par11(1:159,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
h05_r = cat(1,h05(160:300,:),h05(301:360,:),h05(1:159,:));
h05_r(1:2:end) = 0;

Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,par11_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([300 300],[0 -700],'linewidth',1.5,'color','k')
[xx,yy] = find(h05_r==1);
plot(lonData_r(xx),-depthData(yy),'.','color',[.5 .5 .5],'markersize',3)
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','K','Units','normalized','position',[4 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[160,360+160]);
set(gca,'XTick',[160:60:360+160]);
set(gca,'XTickLabel',{'160^oE','140^oW','80^oW','20^oW','40^oE','100^oE','160^oE'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
% title(['PC1.reg.meridional mean temperature (',names,')'])
print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.Tem0-700m.png'],'-dpng','-r300')
%% AMO reg 700-2000m Tem
var = permute(Tsub(:,1:60,:),[3 1 2]);
index = PC1z; pngname = 'PC1';
clear par11 h01 par12 h02
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h01(i,j),t] = reg1_ttest(index,var(:,i,j),0.1,1);
    end
end
%%
lats = 1:60;
close all;
ftsz = 12; ticks = 0.15;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par11,[-0.4:0.01:0.4],ticks,lonData,latData(lats),12)
testdots(h01,[.4 .4 .4],lonData,latData)
% m_line([0:1:360],-46,'linewidth',2,'color','k'); 
% m_line([0:1:360],-61,'linewidth',2,'color','k'); 
% m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
% m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['G:\figures\IAP\Yearly\20231008_AMO_SOdeep\',pngname,'_reg_Tem700-2000m.png'],'-dpng','-r300')

%% PC1/DI .reg. taux tauy slp curlz
index = PC1z; pngname = 'PC1';
d1 = size(taux,1); d2 = size(taux,2); d3 = size(taux,3);
varu = (reshape(taux,d1*d2,d3))';
varv = (reshape(tauy,d1*d2,d3))';
varslp = (reshape(slpadf,d1*d2,d3))';
varcurlz = (reshape(curlz,d1*d2,d3))';
clear par1 h1 par2 h2 par3 h3 par4 h4
parfor i = 1:d1*d2;
    [par1(i),h1(i),t] = reg1_ttest(index,varu(:,i),0.1,1);
    [par2(i),h2(i),t] = reg1_ttest(index,varv(:,i),0.1,1);
    [par3(i),h3(i),t] = reg1_ttest(index,varslp(:,i),0.1,1);
    [par4(i),h4(i),t] = reg1_ttest(index,varcurlz(:,i),0.1,1);
end
par11 = reshape(par1,d1,d2);  par22 = reshape(par2,d1,d2);
par33 = reshape(par3,d1,d2);  par44 = reshape(par4,d1,d2);
%% slp + tau
% close all;
ftsz = 12; ticks1 = 60; ticks2 = 0.08;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
par33(361,:) = par33(1,:);
% contourVARra(var,val,ctr,lon1,lon2,lat1,lat2,lonData,latData)
contourVARra(par33,[-ticks1*5:ticks1/5:ticks1*5],ticks1,0,360,-90,90,lonEmap,latDataE)
%     m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
%     hold on
%     m_contourf(lonEmap,latDataE,par33',[-ticks1*5:ticks1/5:ticks1*5],'linestyle','none');
%     caxis([-ticks1,ticks1]);
%     load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
%     colormap(bl_re5);
%     m_coast('linewidth',1,'color','k');
%     m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
% contourfS(par33,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonEmap,latDataE,12)
% testdots(h3,[.4 .4 .4],lonEmap,latDataE)
% m_line([0:1:360],-46,'linewidth',2,'color','k'); 
% m_line([0:1:360],-61,'linewidth',2,'color','k'); 
% m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
% m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
% set_colorbar([0.83 0.08 0.03 0.88],'Sea Level Pressure  (Pa)',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
% d = 5; dd = 0.0003;
% [x,y] = meshgrid(lonEmap,latDataE);
% umap = par11.*cosd(y'); % zonal wind weight
% m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par22(1:d:end,1:d:end)'./dd,0,...
%     'color','k','linewidth',1,'MaxHeadSize',5) 
% text(1.2,0.1,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
%     'margin',5,'position',[-3.073,1.38,-1])
% m_quiver(50,-20,0.5./dd*cosd(87),0,0,'color','r','linewidth',1,'MaxHeadSize',5) 
% print(Fig,['G:\figures\IAP\Yearly\20231008_AMO_SOdeep\',pngname,'.reg.windstress&slp.png'],'-dpng','-r300')
%% curlz & wind
close all;
ftsz = 12; ticks1 = 8; ticks2 = 0.08;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
par44(361,:) = par44(1,:);
contourfS(par44*10^9,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonEmap,latDataE,12)
testdots(h4,[.4 .4 .4],lonDataE,latDataE)
d = 5; dd = 0.0003;
[x,y] = meshgrid(lonDataE,latDataE);
umap = par11.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par22(1:d:end,1:d:end)'./dd,0,...
    'color','k','linewidth',1,'MaxHeadSize',5) 
% text(1.2,0.1,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
%     'margin',5,'position',[-3.073,1.38,-1])
% m_quiver(10,-60,0.5./dd*cosd(87),0,0,'color','r','linewidth',1,'MaxHeadSize',5)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'Wind Stress Curl  (10^-^9 Pa m^-^1)',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
% print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.windstress&curlz.png'],'-dpng','-r300')
%% taux
close all;
ftsz = 12; ticks1 = 5;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
par11(361,:) = par11(1,:);
contourfS(par11*10^3,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonEmap,latDataE,12)
testdots(h1,[.4 .4 .4],lonDataE,latDataE)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'Zonal Wind Stress  (10^-^3 Pa)',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
% print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.Taux.png'],'-dpng','-r300')
%% tauy
ftsz = 12; ticks1 = 5;
close all;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
par22(361,:) = par22(1,:);
contourfS(par22*10^3,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonEmap,latDataE,12)
testdots(h2,[.4 .4 .4],lonDataE,latDataE)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'TAUY  (10^-^3 Pa)',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
% print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.Tauy.png'],'-dpng','-r300')
%% PC1 reg meridional mean (v wind stress) 46-61S 
lats = [29:44]; names = '46S-61S'; index = PC1z;
latm = latmean(tauy,lats,latData);
% latmc = latmean(Temp,lats,latData); % climotology
% corr with TPI
var = latm';
clear par11 h05
for i = 1:size(var,2);
    [par11(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i),0.05,1);
end
max(par11,[],'all')
min(par11,[],'all')
%
close all;
ticks = 0.1;
par11_r = cat(1,par11(160:300,:),par11(301:360,:),par11(1:159,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
Fig = figure('position',[100 100 800 400]);
plot(lonData_r,par11_r,'g','linewidth',1.5)
hold on
plot(lonData_r,zeros(360),'k','linewidth',1.5)
set(gca,'XLim',[150,360+149+1]);
set(gca,'XTick',[150:60:360+149+1]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-1.8,1.8]*10^-4);
set(gca,'YTick',[-1.8:0.6:1.8]*10^-4,'FontSize',14);
ylabel('Meridional wind stress (Pa)')
% print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\PC1.reg.meridionalmean.Vwindstress.50_60S.png'],'-dpng','-r300')
%% PC1/DI/IPO reg zonal mean (SST)
close all;
reg_zonalmean(sstadf,PC1z,[-0.1,0.15],[-0.1:0.05:0.15],'SST (K)','PC1.reg.zonalmean_SST',latData)
reg_zonalmean(sstadf,DI,[-0.1,0.15],[-0.1:0.05:0.15],'SST (K)','DI.reg.zonalmean_SST',latData)
reg_zonalmean(sstadf,TPIfz,[-0.1,0.15],[-0.1:0.05:0.15],'SST (K)','IPO.reg.zonalmean_SST',latData)
%% PC1/DI/IPO reg zonal mean (net heat flux)
close all;
reg_zonalmean(nhflxadf,PC1z,[-1.6,2.5],[-1.5:0.5:2.5],'Net Heat Flux (W m^-^2)','PC1.reg.zonalmean_NetHeatFlux',latDataE)
reg_zonalmean(nhflxadf,DI,[-1.5,2.5],[-1.5:0.5:2.5],'Net Heat Flux (W m^-^2)','DI.reg.zonalmean_NetHeatFlux',latDataE)
reg_zonalmean(nhflxadf,TPIfz,[-1.5,1.5],[-1.5:0.5:1.5],'Net Heat Flux (W m^-^2)','IPO.reg.zonalmean_NetHeatFlux',latDataE)
%% PC1/DI/IPO reg zonal mean (u wind)
close all;
reg_zonalmean(taux,PC1z,[-2.1,2]*10^-3,[-2:.5:2]*10^-3,'Zonal Wind Stress (Pa)','PC1.reg.zonalmean_Taux',latDataE)
reg_zonalmean(taux,DI,[-2.1,2]*10^-3,[-2:.5:2]*10^-3,'Zonal Wind Stress (Pa)','DI.reg.zonalmean_Taux',latDataE)
reg_zonalmean(taux,TPIfz,[-2.3,2]*10^-3,[-2:.5:2]*10^-3,'Zonal Wind Stress (Pa)','IPO.reg.zonalmean_Taux',latDataE)
%% PC1/DI/IPO reg zonal mean (curlz)
close all;
reg_zonalmean(curlz,PC1z,[-4,4]*10^-9,[-4:1:4]*10^-9,'Wind Stress Curl (Pa/m)','PC1.reg.zonalmean_curlz',latDataE)
reg_zonalmean(curlz,DI,[-4,4]*10^-9,[-4:1:4]*10^-9,'Wind Stress Curl (Pa/m)','DI.reg.zonalmean_curlz',latDataE)
reg_zonalmean(curlz,TPIfz,[-8,8]*10^-9,[-8:2:8]*10^-9,'Wind Stress Curl (Pa/m)','IPO.reg.zonalmean_curlz',latDataE)
%% PC1/DI/IPO reg zonal mean (slp)
close all;
reg_zonalmean(slpadf,PC1z,[-60,50],[-50:10:50],'Sea Level Pressure (Pa)','PC1.reg.zonalmean_SLP',latDataE)
reg_zonalmean(slpadf,DI,[-50,50],[-50:10:50],'Sea Level Pressure (Pa)','DI.reg.zonalmean_SLP',latDataE)
reg_zonalmean(slpadf,TPIfz,[-60,50],[-50:10:50],'Sea Level Pressure (Pa)','IPO.reg.zonalmean_SLP',latDataE)
%% PC1/DI/IPO reg zonal mean (v wind)
close all;
reg_zonalmean(tauy,PC1z,[-2.1,2]*10^-3,[-2:.5:2]*10^-3,'Zonal Wind Stress (Pa)','PC1.reg.zonalmean_Tauy',latDataE)
reg_zonalmean(tauy,DI,[-2.1,2]*10^-3,[-2:.5:2]*10^-3,'Zonal Wind Stress (Pa)','DI.reg.zonalmean_Tauy',latDataE)
reg_zonalmean(tauy,TPIfz,[-2.1,2]*10^-3,[-2:.5:2]*10^-3,'Zonal Wind Stress (Pa)','IPO.reg.zonalmean_Tauy',latDataE)
%% legend label
close all;
Fig = figure('position',[100 100 800 400]);
y = zeros(length(latData),1);
plot(latData,y,'r','linewidth',1.5)
hold on
plot(latData,y,'b','linewidth',1.5)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:10:-35]);
set(gca,'XTickLabel',{'70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',14);
set(gca,'YLim',[-2,2]);
set(gca,'FontSize',14);
legend('Pacific','Atlantic & Indian Ocean','Location','north','Orientation','vertical','position',[0.47,0.69,0.28,0.17])
legend('boxoff')
saveas(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\legend1.png'])


%% u v wind
close all;
ticks = 0.35; ftsz = 12;
Fig = figure('position',[100 100 800 400]);
% m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
hold on
[c h] = m_contourf(lonData,latData,par1(:,1:180)',[-ticks*5:ticks/5:ticks*5],'linestyle','none');
caxis([-ticks,ticks])
hold on
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-ticks:ticks/5:ticks],'position',[0.78 0.09 0.03 0.85],'fontsize',12);
set(ch.Label,'String','m s^-^1 STD ^-^1','Units','normalized','position',[4 0.5 0],'Rotation',-90,'fontsize',ftsz);
m_coast('linewidth',1,'color','k');
m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
d = 6; dd = 0.1;
[x,y] = meshgrid(lonData,latData);
umap = par1.*cosd(y'); % zonal wind weight
m_quiver(x(1:d:end,1:d:end),y(1:d:end,1:d:end),umap(1:d:end,1:d:end)'./dd,par2(1:d:end,1:d:end)'./dd,0,...
    'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5) 
m_text(10,85,'0.5 m/s','fontsize',12,'edgecolor','k','backgroundcolor','w',...
    'margin',5,'position',[-3.073,1.38,-1])
m_quiver(10,87,0.5./dd*cosd(87),0,0,'color',[0.29 0.47 0.05],'linewidth',1,'MaxHeadSize',5)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
% print(Fig,['G:\figures\IAP\Yearly\20221121_IPO_SouthernOcean\PC1_corr_UV.png'],'-dpng','-r300')
%% PC1 .reg.SST
index = PC1z; pngname = 'PC1';
var = permute(sstadf,[3 1 2]);
clear par5 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par5(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i,j),0.1,1);
    end
end
max(par5,[],'all')
min(par5,[],'all')
%
close all;
ftsz = 12; ticks1 = 0.2;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');-45
contourfS(par5,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
testdots(h05,[.4 .4 .4],lonData,latData)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'Sea Surface Temperature (K)',5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
% print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.sst.png'],'-dpng','-r300')
%% PC1/DI/IPO.reg. net heat flux
index = PC1z; pngname = 'PC1';
var = permute(nhflxadf,[3 1 2]);
clear par6 h06
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par6(i,j),h06(i,j),t] = reg1_ttest(index,var(:,i,j),0.1,1);
    end
end
max(par6,[],'all')
min(par6,[],'all')
%
close all;
ftsz = 12; ticks1 = 8;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
par6(361,:) = par6(1,:);
contourfS(par6,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonEmap,latDataE,12)
testdots(h06,[.4 .4 .4],lonDataE,latDataE)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'Net Heat Flux  (W m^-^2)',5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
% print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.nhflx.png'],'-dpng','-r300')
%% PC1 1std composite 
num_PC1p = find(PC1z > 1);
num_PC1n = find(PC1z < -1);
clear evs_PC1p evs_PC1n
[evs_PC1p] = event_index(num_PC1p,5) % consecutive 5 days
[evs_PC1n] = event_index(num_PC1n,5)
clear map1
for s = 1:length(evs_PC1p)
    map1(:,:,s) = mean(Tsub(:,:,evs_PC1p{s}),3);
end
close all;
[Fig] = com_map(mean(map1,3),0.3,lonData,latData(1:90))
% print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\PC1p_0.5std.com.Tem700m.png'],'-dpng','-r300')
clear map2
for s = 1:length(evs_PC1n)
    map2(:,:,s) = mean(Tsub(:,:,evs_PC1n{s}),3);
end
[Fig] = com_map(mean(map2,3),0.3,lonData,latData(1:90))
% print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\PC1n_0.5std.com.Tem700m.png'],'-dpng','-r300')
%% PC1 1 std composite (meridional mean Tem)
clear map1
for s = 1:length(evs_PC1n)
    map1(:,:,s) = mean(latm(:,:,evs_PC1n{s}),3);
end
close all;
ticks = 0.2;
map = mean(map1,3);
map_r = cat(1,map(160:300,:),map(301:360,:),map(1:159,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
Fig = figure('position',[100 100 800 400]);
contourf(lonData_r,-depthData,map_r',[-ticks*10:ticks/10:ticks*10],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([300 300],[0 -700],'linewidth',1.5,'color','k')
line([150 509],[-700 -700],'linewidth',1.5,'color','k')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar('location','eastoutside');
set(ch,'Ticks',[-ticks:ticks/5:ticks]);
% set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
set(ch.Label,'String','K','Units','normalized','position',[4 0.5 0],'Rotation',-90,'fontsize',ftsz);
set(gca,'XLim',[150,360+149]);
set(gca,'XTick',[150:60:360+149]);
set(gca,'XTickLabel',{'150^oE','150^oW','90^oW','30^oW','30^oE','90^oE','150^oE'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')
print(Fig,['G:\figures\IAP\Yearly\20230206_IPO_SouthernOcean\PC1n_1std.com.meridionalmeanT.png'],'-dpng','-r300')
%% DI/TPI corr other time series (scatter)
lats = [137:152];  % ERA5 lat
var = nhflxadf; varname = 'NHFLX'; titlename = 'Net Heat Flux';
var_r = cat(1,var(160:300,:,:),var(301:360,:,:),var(1:159,:,:)); 
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
[tz1 t1] = areamean(var_r,1:141,lats,latDataE); % Pac
[tz2 t2] = areamean(var_r,142:360,lats,latDataE); % Atl & IO
%
% tx1 = spacz(1:72); tx2 = tz1; pocname = 'Pac';
% tx1 = siaz(1:72); tx2 = tz2; pocname = 'IOAtl';
% tx1 = PC1z(1:72); tx2 = zscore(t1-t2); pocname = 'PC1';
% tx1 = DI; tx2 = zscore(t1-t2); xname = 'PmAI'; pocname = 'DI';
tx1 = TPIfz; tx2 = zscore(t1-t2); xname = 'IPO'; pocname = 'IPO';
%
xlimv = [-3,3];
ylimv = [-4,3];
[r,p,n_eff2] = corr_eff(tx1,tx2,0.05)  %  90% confidence  
close all;
Fig = figure('position',[700 100 400 400]);
scatter(tx1,tx2,'*','Color',[0,46,166]/255,'LineWidth',2);
hold on
yb = polyfit(tx1,tx2,1);
plot(tx1,polyval(yb,tx1),'k','linewidth',3)
text(-2.5,-2.8,['Corr = ',num2str(roundn(r,-2))],'fontsize',12)
text(-2.5,-3.4,['Slope = ',num2str(roundn(yb(1),-2))],'fontsize',12)
set(gca,'XLim',xlimv);
set(gca,'FontSize',12);
set(gca,'YLim',ylimv);
xlabel(xname,'FontSize',12),ylabel(titlename,'FontSize',12);
print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pocname,'.regwith.',varname,'.scatter.png'],'-dpng','-r300')

%% MOHT calculate
% ro*cp*OHC*Vy
% Vy = My/ro = -taux/f/ro 
datadir='G:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
f = sw_f(latData);
tauxraw = cat(1,fliplr(taux_raw(2:360,1:end-1,:)),fliplr(taux_raw(1,1:end-1,:)));
My = -tauxraw./f'; % Ekman 质量输运
ro = 1020; % sea water density kg/m3 (Antonov et al. 2004)
cp = 4187; % specific heat of water J/kg/K (Antonov et al. 2004)
Vy = My/ro;
Vy(:,85:96,:) = nan; % 赤道不适用ekman输送
londist = 2*pi*6.371393*10^6*cos(latData*pi/180)/360; % distance between 2 longitudes 
levdist(1) = 1; levdist(2:27) = dweit; % z distance
clear OHC OHCT MOHT
OHC = permute(ro*cp*sum((Temp(:,:,1:27,:)).*permute(levdist,[3 1 2]),3),[1 2 4 3])/700;
OHCT = OHC.*Vy;
OHCTw = OHCT.*londist';
OHCTw_r = cat(1,OHCTw(160:300,:,:),OHCTw(301:360,:,:),OHCTw(1:159,:,:));
MOHT1 = permute(nansum(OHCTw_r(1:141,:,:),1),[2 3 1]); % Pac unit:W
MOHT2 = permute(nansum(OHCTw_r(142:360,:,:),1),[2 3 1]); % IO+Atl unit:W
%%
close all;
MOHT = permute(nansum(OHCTw_r,1),[2 3 1]); % unit:W
Fig = figure('position',[100 100 600 300])
map = mean(MOHT,2)/10^15;
% map(85:96,:) = nan;
plot(latData,map,'linewidth',1.5)
set(gca,'XLim',[-90,90],'YLim',[-2,2],'XGrid','on','YGrid','on');
set(gca,'XTick',[-80:20:80],'YTick',[-2:0.5:2]);
ylabel('MOHT (PW)');xlabel('Latitude')
% print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\MOHT_climotology.png'],'-dpng','-r300')
%% MOHT anomaly,detrend, and 8-yr filter
MOHT1a = MOHT1-mean(MOHT1,2);
MOHT2a = MOHT2-mean(MOHT2,2);
dT = 1; cf = 1/8;
clear MOHT1ad MOHT1adf MOHT2ad MOHT2adf
for i = 1:180
    MOHT1ad(:,i) = detrend(MOHT1a(i,:)'); % linear detrend
    MOHT1adf(:,i) = lanczosfilter(MOHT1ad(:,i),dT,cf,[],'low'); % 8 year filtered
    MOHT2ad(:,i) = detrend(MOHT2a(i,:)'); % linear detrend
    MOHT2adf(:,i) = lanczosfilter(MOHT2ad(:,i),dT,cf,[],'low'); % 8 year filtered
end
MOHT1adf(1:4,:) = []; MOHT1adf(end-3:end,:) = [];
MOHT2adf(1:4,:) = []; MOHT2adf(end-3:end,:) = [];
%%
taux_r = cat(1,taux(160:300,:,:),taux(301:360,:,:),taux(1:159,:,:));
lonmh1 = permute(nanmean(taux_r(1:141,:,:),1),[2 3 1]);
lonmh2 = permute(nanmean(taux_r(142:360,:,:),1),[2 3 1]);
index = DI; pngname = 'DI';
clear par7 h07 par8 h08 par7t h07t par8t h08t
for i = 1:size(MOHT1adf,2); 
    [par7(i),h07(i),t] = reg1_ttest(index,MOHT1adf(:,i),0.05,0);
    [par8(i),h08(i),t] = reg1_ttest(index,MOHT2adf(:,i),0.05,0);
    [par7t(i),h07t(i),t] = reg1_ttest(index,lonmh1(i,:)',0.05,0);
    [par8t(i),h08t(i),t] = reg1_ttest(index,lonmh2(i,:)',0.05,0);
end
%%
close all;
Fig = figure('position',[100 100 800 400]);
p1 = plot(latData,par7/10^15,'r-','linewidth',1.5);
hold on
p2 = plot(latData,par8/10^15,'b-','linewidth',1.5);
plot(-61*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
plot(-46*ones(20),[-.025:.057/19:.032],'k--','linewidth',1.5);
h07(1:30) = 0; h08(1:30) = 0;
scatter(latData(find(h07 == 1)),par7(find(h07 == 1))/10^15,'o','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h08 == 1)),par8(find(h08 == 1))/10^15,'o','MarkerEdgeColor','b','linewidth',1)
plot(latData,zeros(1,length(latData)),'k','linewidth',1.5)
% set(gca,'YLim',[-.025,.032],'YTick',[-.03:.01:.03]); % IPO
set(gca,'YLim',[-.02,.02],'YTick',[-.03:.01:.03]); % DI
set(gca,'YTick',yticks,'FontSize',14);
ylabel('MOHT (PW)')
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
% print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.zonalmean_MOHT_Ekman.png'],'-dpng','-r300')

yyaxis right
p3 = plot(latDataE(1:180),par7t,'-','color',[0.90,0.51,0.31],'linewidth',1.5); % ERA lat
p4 = plot(latDataE(1:180),par8t,'-','color',[0.17,0.73,0.69],'linewidth',1.5); % ERA lat
scatter(latDataE(find(h07t == 1)),par7t(find(h07t == 1)),'o','MarkerEdgeColor',[0.90,0.51,0.31],'linewidth',1)
scatter(latDataE(find(h08t == 1)),par8t(find(h08t == 1)),'o','MarkerEdgeColor',[0.17,0.73,0.69],'linewidth',1)
set(gca,'YLim',[-.0023,.0023],'YTick',[-.002:.001:.002],'ycolor','k');
set(gca,'YTick',yticks,'FontSize',14);
ylabel('Zonal Wind Stress (Pa)','color','k','rotation',270,'position',[-31.8 0 -1]);
print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.reg.zonalmean_MOHT&Taux.png'],'-dpng','-r300')
%%
Fig = figure('position',[100 100 800 400]);
p1 = plot(latData,par7/10^15,'r-','linewidth',1.5);
hold on
p2 = plot(latData,par8/10^15,'b-','linewidth',1.5);
set(gca,'YLim',[-0.4,0.4],'FontSize',14);
legend([p1,p2],'Pacific MOHT','Atlantic & Indian Ocean MOHT','Location','north','Orientation','vertical','position',[0.49,0.75,0.35,0.12])
legend('boxoff')
print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\legend2.png'],'-dpng','-r300')

Fig = figure('position',[100 100 800 400]);
p3 = plot(latDataE(1:180),par7t,'-','color',[0.90,0.51,0.31],'linewidth',1.5); % ERA lat
hold on
p4 = plot(latDataE(1:180),par8t,'-','color',[0.17,0.73,0.69],'linewidth',1.5); % ERA lat
set(gca,'YLim',[-.003 .003],'FontSize',14);
legend([p3,p4],'Pacific Zonal Wind Stress','Atlantic & Indian Ocean Zonal Wind Stress','Location','north','Orientation','vertical','position',[0.49,0.75,0.35,0.12])
legend('boxoff')
print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\legend3.png'],'-dpng','-r300')


%% (teleconnection) IPO reg Z200 zonal deviation
varz = permute(z200adf,[3 1 2]);
index = TPIfz; pngname = 'IPO';
clear par1 h1 par2 h2 par3 h3 par4 h4
for i = 1:size(varz,2);
    for j = 1:size(varz,3);
        [par1(i,j),h1(i,j),t] = reg1_ttest(index,varz(:,i,j),0.05,1);
    end
end
%
close all;
ticks1 = 10; ftsz = 12;
Fig = figure('position',[100 100 650 300])
contourVARra(par1,[-ticks1*5:ticks1/5:ticks1*5],ticks1,0,360,-90,90,lonDataE,latDataE)
testdots(h1,[.4 .4 .4],lonDataE,latDataE)
set(gca,'position',[0.08,0.11,0.775,0.815]);
set_colorbar([0.87 0.095 0.026 0.832],'Z200 (gpm)',4,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
print(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\IPO.reg.Z200za.png'],'-dpng','-r300')














function [adf] = caladf(var)
% var lon*lat*time
vara = double(var - mean(var(:,:,42:71),3)); % remove climatology from 1981 to 2010
varr = permute(vara,[3 1 2]);
dT = 1; cf = 1/8;
for i = 1:size(varr,2);
    for j = 1:size(varr,3);
        ad(:,i,j) = detrend(varr(:,i,j)); % linear detrend
        adf(:,i,j) = lanczosfilter(ad(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
adf = permute(adf,[2 3 1]);
end
function [ad] = detrend3d(var)
    varr = permute(var,[3 1 2]);
    for i = 1:size(varr,2);
        for j = 1:size(varr,3);
            ad(:,i,j) = detrend(varr(:,i,j)); % linear detrend
        end
    end
    ad = permute(ad,[2 3 1]);
end
function [af] = filter3d(var,dT,cf)
% var lon*lat*time
varr = permute(var,[3 1 2]);
for i = 1:size(varr,2);
    for j = 1:size(varr,3);
        af(:,i,j) = lanczosfilter(varr(:,i,j),dT,cf,[],'low'); % 8 year filtered
    end
end
af = permute(af,[2 3 1]);
end
function [varar] = LoadNOAA(filename,varname)
% filename is where the file locates
% varname is name of the variable in raw file
% filename = 'G:\data\NOAA\20CRv3\20thC_ReanV3\Monthlies\miscSI-MO\prmsl.mon.mean.nc'; % 1836-2015
    ncid=netcdf.open(filename,'NOWRITE');
    ncdisp(filename);
    var = ncread(filename,varname);
    var = permute(mean(reshape(var(:,:,1249:end),[360 181 12 912/12]),3),[1 2 4 3]); % yearly 1940-2015
    vara = double(var - mean(var(:,:,:),3)); % remove climatology from 1981 to 2010
%42:71
    varar = permute(vara,[3 1 2]);
end
function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',90,'radius',90,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re5.mat');
    colormap(bl_re5);
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function contourfSPolar(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',55,'rotangle',0);
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
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',12);
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
function [Fig] = diffdepth(Temcom,lonData,latData,lats,val,ticks)
    Fig = figure('position',[100 100 600 300]);
%     contourVARra(Temcom(:,:,nlevel),[-val:ticks/10:val],1,0,360,-90,90,lonData,latData);
    contourfSPolar(Temcom,[-val:ticks/10:val],ticks,lonData,latData(lats),12)
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
function [res] = gradient_lat(lat,var)
% 此函数用于计算沿经向的梯度
% Input:
% var is lat*lon 
% lat is 单调递增，自南向北
% Output:
% res 即为梯度结果。
% degree to radian
rad=pi/180;
% computation of curl
[lt, ln]=size(var);
a=diff(lat);
aa=NaN*ones(length(a)-1,1);
for ii=1:length(a)-1
    if (a(ii) == a(ii+1))
        aa(ii)=a(ii);
    else
        error('Latitude difference is not consistance')
    end % endif
    dlat=mean(aa);
end % endfor
clear ii
deltay=dlat*111176;
res=NaN(lt, ln);
% Centeral difference method in x and y 
for ii=2:lt-1
    for jj=1:ln
        res(ii, jj)= (var(ii+1, jj)-var(ii-1, jj))/(2*deltay) ;
    end % endfor
end % endfor
clear ii jj
% Forward difference method in x and y 
for jj=1:ln
    res(1, jj)= (var(2, jj)-var(1, jj))/deltay ; 
end 
clear ii jj
% Backward difference method in x and y
for jj=1:ln
    res(lt, jj)= (var(lt, jj)-var(lt-1, jj))/deltay ;
end
clear ii jj
end
function reg_zonalmean(nhflx,index,ylimval,yticks,titlename,pngname,latData)
lonmh1 = permute(nanmean(nhflx(160:300,:,:),1),[2 3 1]);
nhflx_r = cat(1,nhflx(300:360,:,:),nhflx(1:299,:,:));
lonmh2 = permute(nanmean(nhflx_r(1:220,:,:),1),[2 3 1]);
% reg with PC1
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
Fig = figure('position',[100 100 800 400]);
plot(latData,par11,'r','linewidth',1.5)
hold on
plot(latData,par12,'b','linewidth',1.5)
plot(latData,zeros(1,length(latData)),'k','linewidth',1.5)
plot(-61*ones(20),[ylimval(1):(ylimval(2)-ylimval(1))/19:ylimval(2)],'k--','linewidth',1.5);
plot(-46*ones(20),[ylimval(1):(ylimval(2)-ylimval(1))/19:ylimval(2)],'k--','linewidth',1.5);
scatter(latData(find(h01 == 1)),par11(find(h01 == 1)),'o','MarkerEdgeColor','r','linewidth',1)
scatter(latData(find(h02 == 1)),par12(find(h02 == 1)),'o','MarkerEdgeColor','b','linewidth',1)
set(gca,'XLim',[-70,-35]);
set(gca,'XTick',[-70:5:-35]);
set(gca,'XTickLabel',{'70^oS','65^oS','60^oS','55^oS','50^oS','45^oS','40^oS','35^oS'},'FontSize',14);
set(gca,'YLim',ylimval);
set(gca,'YTick',yticks,'FontSize',14);
ylabel(titlename)
% legend('Pacific','Atlantic & Indian Ocean','Location','north','Orientation','vertical','position',[0.61,0.75,0.28,0.17],'edgecolor','w')
saveas(Fig,['G:\figures\IAP\Yearly\20230606_IPO_SouthernOcean_46_61S\',pngname,'.png'])
end

function [eyr] = event_index(var,days)
% This function is used to find events which is defined by consecutive 
% days/years.
% Input var is all days/years that satisfied the condition. e.g. 1std
% Input days is the number of defined consecutive days/years
% Output eyr is the choosed events
enum = 1; k = 0;
for s = 2:length(var);
    diff = var(s)-var(s-1);
    if diff == 1 
        k = k + 1;
        if k >= days-1
            eyr{enum} = var(s-k:s);
        end
    else 
        if k >= days-1
            eyr{enum} = var(s-1-k:s-1);
            enum = enum + 1;
        end
        k = 0;
    end
end
end

function [Fig] = com_map(var,ticks1,lonData,latData)
ftsz = 12;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(var,[-ticks1*5:ticks1/5:ticks1*5],ticks1,lonData,latData,12)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'K',4.5,ftsz,[-ticks1:ticks1/5:ticks1],[-ticks1:ticks1/5:ticks1])
end
function testdots(h,clor,lonData,latData)
    dots = mod(find(h == 1),length(lonData));
    dots(find(dots == 0)) = length(lonData);
    m_plot(lonData(dots),latData(ceil(find(h == 1)/length(lonData))),'.','markersize',2,'color',clor); % F-test dots
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