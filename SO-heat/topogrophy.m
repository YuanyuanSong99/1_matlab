clc,clear,close all;
filename1 = '/Users/yysong/Downloads/ETOPO_2022_v1_60s_N90W180_bed.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
lat = ncread(filename1,'lat');
lon = ncread(filename1,'lon');
z = ncread(filename1,'z');
%%
close all
ftsz = 20;
Fig = figure('position',[10 50 600 400]);
ax = axes('Position',[0.0005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
zmap = z;
zmap(find(z>0)) = nan;
map = cat(1,zmap(1:end,:),zmap(1,:));
lonex = cat(1,lon(1:end),lon(1));
vlim = 6000;
m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
hold on
m_contourf(lonex(1:50:end),lat(1:50:end),map(1:50:end,1:50:end)',[-vlim:vlim/99:0],'linestyle','none');
caxis([-vlim,0])
%
h = demcmap([-100,100],100);
colormap(h(1:50,:))
% m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','none','xaxisloc','top');
set_colorbar([0.81 0.08 0.03 0.88],'Sea floor Topography (m)',6,ftsz,[-6000:1000:0],[-6000:1000:0])
print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240812/SeaFloorTopography.png'],'-dpng','-r300')
%% wind climatology
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
datadir='/Volumes/Togo4T/data/ERA5_1x1/windstress/'; %指定批量数据所在的文件夹
filelist=dir([datadir,'*.nc']); %指定批量数据的类型
ncdisp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
ncid=netcdf.open(filetemp,'NC_NOWRITE');
lonDataE = ncread(filetemp,'longitude'); % 0-360
lonEmap = lonDataE; lonEmap(361) = lonEmap(1);
latDataE = ncread(filetemp,'latitude'); % 90 - -90
%
k=length(filelist);
for s=1:k-3;  % 1940-2020
    s
    filename1=[datadir,filelist(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    metss = ncread(filename1,'metss');     uwst(:,:,s) = mean(metss,3);   % yearly
    mntss = ncread(filename1,'mntss');     vwst(:,:,s) = mean(mntss,3);   % yearly 
    curlzE(:,:,s) = ra_windstrcurl(latDataE,lonDataE,uwst(:,:,s)',vwst(:,:,s)',1); % lat*lon
    netcdf.close(ncid);  %关闭文件    
end
curlz_Er = permute(curlzE,[2 1 3]);
ucli = nanmean(uwst,3);
vcli = nanmean(vwst,3);
uvcli = sqrt(ucli.^2+vcli.^2);
curlz_cli = nanmean(curlz_Er,3);
%% zonal wind stress
close all
ftsz = 20;
Fig = figure('position',[10 50 600 400]);
ax = axes('Position',[0.0005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
map = cat(1,ucli(1:end,:),ucli(1,:));
lonData_r = cat(1,lonDataE(1:end),lonDataE(1));
ticks = 0.2;
contourfS(map,[-ticks*5:ticks/20:ticks*5],ticks,lonData_r,latDataE,20)
% m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
set_colorbar([0.81 0.08 0.03 0.88],'Zonal wind stress (Pa)',5.5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240812/ZonalWindStressClimo.png'],'-dpng','-r300')

%% wind stress curl
close all
ftsz = 20;
Fig = figure('position',[10 50 600 400]);
ax = axes('Position',[0.0005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
map = cat(1,curlz_cli(1:end,:),curlz_cli(1,:));
lonData_r = cat(1,lonDataE(1:end),lonDataE(1));
ticks = 2;
contourfS(map,[-ticks*5:ticks/20:ticks*5]*10^-7,ticks*10^-7,lonData_r,latDataE,20)
% m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
set_colorbar([0.81 0.08 0.03 0.88],'Wind stress curl (10^-^7 Pa m^-^1)',5.5,ftsz,[-ticks:ticks/5:ticks]*10^-7,[-ticks:ticks/5:ticks])
print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240812/WindStressCurl_Climo.png'],'-dpng','-r300')



function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','none','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end