clc,clear,close all;
filename1 = '/Users/yysong/Desktop/data/aviso.u.1993to2022.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
lat = ncread(filename1,'latitude');
lon = ncread(filename1,'longitude');
u = ncread(filename1,'ugos');
v = ncread('/Users/yysong/Desktop/data/aviso.v.1993to2022.nc','vgos');
uvcli = sqrt(u.^2+v.^2);
%% 
close all
ftsz = 20;
Fig = figure('position',[10 50 600 400]);
ax = axes('Position',[0.0005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
map = cat(1,uvcli(1:end,:),uvcli(1,:));
%map = cat(1,u(1:end,:),u(1,:));
lonData_r = cat(1,lon(1:end),lon(1));
ticks = 0.2;
contourfS(map,[0:ticks/20:ticks*5],ticks,lonData_r,lat,20)
% m_line([0:1:360],-35,'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line([0:1:360],-55,'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line(-70,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
% m_line(150,[-55:1:-35],'linewidth',2.5,'color',[.4 .4 .4]); 
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
set_colorbar([0.81 0.08 0.03 0.88],'Surface geostrophic velocity (m/s)',5.5,ftsz,[0:ticks/5:ticks],[0:ticks/5:ticks])
print(Fig,['/Users/yysong/Desktop/figures/SOtrend_all_Data/20240812/SurfaceVelocityClimo.png'],'-dpng','-r300')


function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([0,ctr]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4(1:10,:)));
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','none','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end