clc,clear,close all;
addpath(genpath('D:\1_matlab\help'));
%%
datadir='D:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% anomaly 
k=length(filelist);
clear Temp
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = double(ncread(filename1,'/temp in degC')); %读入变量
end
%Tempa = double(Temp - nanmean(Temp(:,:,:,42:71),4)); % remove climatology from 1981-2010 
datadir='D:\data\IAP\salinity\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
k=length(filelist);
clear salt
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    salt(:,:,:,s) = double(ncread(filename1,'/absolute salinity')); %读入变量
end
%% Monthly
datadir='D:\data\IAP\temperature\monthly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'CZ16*_09.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = ncread(filetemp,'lon');
latData = ncread(filetemp,'lat');
depthData = ncread(filetemp,'depth_std');
% anomaly 
k=length(filelist);
clear Temp
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = double(ncread(filename1,'temp')); %读入变量
end
Temp = permute(Temp,[2 3 1 4]);
datadir='D:\data\IAP\salinity\Monthly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'CZ16*_09.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
k=length(filelist);
clear salt
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    salt(:,:,:,s) = double(ncread(filename1,'salinity')); %读入变量
end
salt = permute(salt,[2 3 1 4]);
%% Season July-September
datadir='D:\data\IAP\temperature\monthly\'; %指定批量数据所在的文件夹
filelist9=dir([datadir,'CZ16*_09.nc']); %指定批量数据的类型
filelist8=dir([datadir,'CZ16*_08.nc']); %指定批量数据的类型
% filelist7=dir([datadir,'CZ16*_07.nc']); %指定批量数据的类型
ncdisp([datadir,filelist9(1).name]);
filetemp=[datadir,filelist9(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = ncread(filetemp,'lon');
latData = ncread(filetemp,'lat');
depthData = ncread(filetemp,'depth_std');
% anomaly 
k=length(filelist9);
clear Tempall
for s=1:k;   
    s
%     filename7=[datadir,filelist7(s).name];
    filename8=[datadir,filelist8(s).name];
    filename9=[datadir,filelist9(s).name];
%     Tempall(:,:,:,s,1) = double(ncread(filename7,'temp')); %读入变量
    Tempall(:,:,:,s,1) = double(ncread(filename8,'temp')); %读入变量
    Tempall(:,:,:,s,2) = double(ncread(filename9,'temp')); %读入变量
end
Tempsm = mean(Tempall,5);
Temp = permute(Tempsm,[2 3 1 4]);
clear Tempall Tempsm
datadir='D:\data\IAP\salinity\Monthly\'; %指定批量数据所在的文件夹
filelist9=dir([datadir,'CZ16*_09.nc']); %指定批量数据的类型
filelist8=dir([datadir,'CZ16*_08.nc']); %指定批量数据的类型
% filelist7=dir([datadir,'CZ16*_07.nc']); %指定批量数据的类型
ncdisp([datadir,filelist9(1).name]);
filetemp=[datadir,filelist9(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
k=length(filelist9);
clear saltall
for s=1:k;   
    s
%     filename7=[datadir,filelist7(s).name];
    filename8=[datadir,filelist8(s).name];
    filename9=[datadir,filelist9(s).name];
%     saltall(:,:,:,s,1) = double(ncread(filename7,'salinity')); %读入变量
    saltall(:,:,:,s,1) = double(ncread(filename8,'salinity')); %读入变量
    saltall(:,:,:,s,2) = double(ncread(filename9,'salinity')); %读入变量
end
saltsm = mean(saltall,5);
salt = permute(saltsm,[2 3 1 4]);
clear saltall saltsm
%%
% % pressure
% [dep,lat] = meshgrid(depthData,latData);
% pres = sw_pres(dep,lat);
% size(pres)
% % potential temperature
% Pres = double(permute(repmat(permute(pres,[3 4 1 2]),360,1),[1 3 4 2])); % convert pres as the same dimension of salt
% pr = Pres(:,:,41); % reference pressure 2000m
% PR = double(permute(repmat(permute(pr,[3 4 1 2]),41,1),[3 4 1 2])); % convert pr as the same dimension of salt
% for s = 1:k
%     PT(:,:,:,s) = sw_ptmp(salt(:,:,:,s),Temp(:,:,:,s),Pres,PR); % potential temperature
%     Prho(:,:,:,s) = sw_dens(salt(:,:,:,s),PT(:,:,:,s),PR); % potential density
% end
% %
% nlevels = length(depthData);
% PTr = reshape(permute(PT,[3 1 2 4]),nlevels,360*180*81);
% Prhor = reshape(permute(Prho(:,:,:,53:72),[3 1 2 4]),nlevels,360*180*nyr);
%% 
nlevels = length(depthData);
startyr = 1940; endyr = 2020;
nyr = endyr-startyr+1;
Tempr = reshape(permute(Temp(:,:,:,startyr-1939:endyr-1939),[3 1 2 4]),nlevels,360*180*nyr);
saltr = reshape(permute(salt(:,:,:,startyr-1939:endyr-1939),[3 1 2 4]),nlevels,360*180*nyr);
clear rho Temp salt
for s = 1:360*180*nyr;
    rho(:,s) = density(Tempr(:,s), saltr(:,s), depthData);
end
clear saltr
%%
dro1(1:2,:) = zeros(2,360*180*nyr);
dro2 = rho(3:nlevels,:) - rho(3,:); % ref depth 10m
dro = cat(1,dro1,dro2);
mld = zeros(nlevels,360*180*nyr);
mld(find(dro > 0.03)) = 1; % delta rho > 0.03
%
clear mldnum mldlev dro* rho
mldnum = ones(1,360*180*nyr)*nan;
mldlev = ones(1,360*180*nyr)*nan;
parfor s = 1:360*180*nyr;
    if isempty(find(mld(:,s),1)) == 0
    mldnum(s) = find(mld(:,s),1);
    mldlev(s) = ( depthData(mldnum(s))+depthData(mldnum(s)-1) )/ 2;
    elseif isempty(find(isnan(Tempr(:,s)),1)) == 0
        mldlev(s) = depthData(find(isnan(Tempr(:,s)),1));
    else
        mldlev(s) = depthData(end);
    end
end
mldlev(find(isnan(Tempr(1,:)) == 1)) =nan;
clear Tempr
%
mldr = reshape(mldlev,[360 180 nyr]);
mldrcli = nanmean(mldr,3);
mldlatmcli = latmean(mldrcli,lats,latData); % cli
mldlatmclir = -cat(1,mldlatmcli(160:360,:),mldlatmcli(1:159,:)); 
%%
close all
collev = [0,20,40,60,80,100,150,200,250,300,400,500,700,900,1100,1500];
ncol = length(collev);
intc = floor(256/ncol);
jet = jet;
jetsub = jet(1:intc:ncol*intc,:);
[~,ZI2]=histc(mldrcli,collev);
% ZI2=ZI2-1;
Fig = figure('position',[100 100 650 300])
m_proj('equidistant Cylindrical','lat',[-90 90],'long',[0 360]);
hold on
m_contourf(lonData,latData,ZI2','linestyle','none');
caxis([0,16]);
colormap(jetsub);
m_coast('patch',[1,1,1],'linewidth',1,'edgecolor','k');
m_grid('xtick',7,'ytick',7,'box','on','linestyle','none','fontsize',10);
cbar = colorbar;
set(cbar,'Ticks',0:length(collev)-1,'TickLabels',collev) ;
set(cbar.Label,'String','Mixed Layer Depth (m)','Units','normalized','position',[3.5 0.5 0],'Rotation',-90,'fontsize',10);
% print(Fig,['D:\figures\IAP\climo\MixedLayerDepth_JAS.png'],'-dpng','-r300')
%% 
clear dro* mld mldlev mldnum rho salt Temp
mlda = double(mldr - nanmean(mldr(:,:,42:71),3)); % remove climatology from 1981-2010 
mlda1 = double(permute(mlda,[3 1 2]));
dT = 1; % interval
cf = 1/8;
mlda1 = reshape(mlda1,[81,360*180]);
clear mldadl mldadlf
parfor i = 1:size(mlda1,2);
    mldadl(:,i) = detrend(mlda1(:,i)); % linear detrend
    mldadlf(:,i) = lanczosfilter(mldadl(:,i),dT,cf,[],'low'); % 8 year filtered
end
mldadlf = permute(reshape(mldadlf,[81,360,180]),[2 3 1]); 
mldadlf(:,:,1:4) = []; mldadlf(:,:,end-3:end) = [];
%% load 8-yr filtered temperature
lats = [29:44]; names = '46S-61S'; index = DI; pngname = 'DI';
latm = latmean(Tempadlf,lats,latData);
% corr with 
var = permute(latm,[3 1 2]);
clear par11 h05
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par11(i,j),h05(i,j),t] = reg1_ttest(index,var(:,i,j),0.1,1);
    end
end
% mld
% mldlatm = latmean(mldadlf,lats,latData); 
% var = mldlatm';
% clear pm hm
% for i = 1:size(var,2);
%     [pm(i),hm(i),t] = reg1_ttest(index,var(:,i),0.1,1);
% end
% max(pm,[],'all')
% min(pm,[],'all')
%%
close all;
ticks = 0.1; ftsz=12;
par11_r = cat(1,par11(160:300,:),par11(301:360,:),par11(1:159,:));
lonData_r = [lonData(160:300);lonData(301:360);lonData(1:159)+360];
h05_r = cat(1,h05(160:300,:),h05(301:360,:),h05(1:159,:));
h05_r(1:2:end) = 0;
mldlatmclir = -cat(1,mldlatmcli(160:360,:),mldlatmcli(1:159,:)); 

Fig = figure('position',[100 100 800 400]);
ax = axes('Position',[0.08 0.1 0.72 0.87],'fontsize',14,'box','on');
caxval = 0.1;
[c h] = contourf(lonData_r,-depthData(1:nlev),par11_r',[-1:caxval/10:1],'linestyle','none');
caxis([-caxval,caxval])
hold on
line([300 300],[0 -700],'linewidth',1.5,'color','k');
plot(lonData_r,mldlatmclir,'linewidth',1.5,'color','k','LineStyle','--')
[xx,yy] = find(h05_r==1);
plot(lonData_r(xx),-depthData(yy),'.','color',[.5 .5 .5],'markersize',3);
load('D:/1_matlab/help/colorbar_mat/bl_re5.mat');
colormap(bl_re5);
ch = colorbar;
set(ch,'Ticks',[-caxval:caxval/5:caxval],'position',[0.83 0.1 0.02 0.87]);
set(ch.Label,'String','K','Units','normalized','position',[4.5 0.5 0],'Rotation',-90,'fontsize',12);
set(gca,'XLim',[160,360+160]);
set(gca,'XTick',[160:60:360+160]);
set(gca,'XTickLabel',{'160^oE','140^oW','80^oW','20^oW','40^oE','100^oE','160^oE'},'FontSize',14);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',14);
set(gca,'YTickLabel',[700:-100:0],'FontSize',14);
ylabel('Depth (m)')

% yyaxis right
% maph = pm; 
% maphr = cat(2,maph(160:360),maph(1:159));
% plot(lonData_r,-maphr,'linewidth',2.5,'color',[.5 .5 .5],'LineStyle','--')
% set(gca,'YLim',[-11.5,4.6],'YTick',[-11.5:2.3:4.6],'YColor',[.5 .5 .5]); % yearly
% % set(gca,'YLim',[-60,24],'YTick',[-60:12:24],'YColor',[.5 .5 .5]);
% yh = ylabel('Mixed Layer Depth (m)')
% set(yh,'Rotation',270,'position',[562,-1.73,0]); % yearly
% set(yh,'Rotation',270,'position',[562,-12,0]); 
print(Fig,['D:\figures\IAP\Yearly\20240130_IPO_SO\',pngname,'.reg.MeridionalMeanTem_YEARLYmixedlayer.png'],'-dpng','-r300')
% print(Fig,['D:\figures\IAP\Yearly\20240130_IPO_SO\',pngname,'.reg.MeridionalMeanTem_AugSepmixedlayer.png'],'-dpng','-r300')
%%
var = permute(mldadlf(:,1:55,:),[3 1 2]);
index = TPIfz; pngname = 'IPO';
clear par1 h01 par12 h02
for i = 1:size(var,2);
    for j = 1:size(var,3);
        [par1(i,j),h01(i,j),t] = reg1_ttest(index,var(:,i,j),0.1,1);
    end
end
max(par1,[],'all')
min(par1,[],'all')
%
lats = 1:55;
close all;
ftsz = 12; ticks = 10;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(par1,[-ticks*5:ticks/10:ticks*5],ticks,lonData,latData(lats),12)
testdots(h01,[.4 .4 .4],lonData,latData)
m_line([0:1:360],-46,'linewidth',2,'color','k'); 
m_line([0:1:360],-61,'linewidth',2,'color','k'); 
m_line(-60,[-61:1:-46],'linewidth',2,'color','k'); 
m_line(160,[-61:1:-46],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'Mixed Layer Depth (m)',4.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['D:\figures\IAP\Yearly\20240130_IPO_SO\',pngname,'_reg_MixedLayerDepth_yearly.png'],'-dpng','-r300')
% print(Fig,['D:\figures\IAP\Yearly\20240130_IPO_SO\',pngname,'_reg_MixedLayerDepth_AugSep.png'],'-dpng','-r300')
% 

function [ts] = latmean(var,lats,latData)
% average along longitude (all latitudes)
% var is lon*lat*depth*time
% ts is lon*depth*time
    var1 = var(:,lats,:,:); 
    var2 = var(:,lats,:,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = permute(nansum((cos(latData(lats)'/180*pi)).*var1,2)./nansum(cos(latData(lats)'/180*pi).*var2,2),[1 3 4 2]);
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