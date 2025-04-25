clc,clear,close all;
addpath G:\1_matlab\help;
addpath G:\1_matlab\help\seawater\;
load("MatFile\lonData.mat");
load("MatFile\latData.mat");
load("MatFile\depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
dweit = depthData(2:end)-depthData(1:end-1);
ncid=netcdf.open('H:\CESM-post\LE\Temperature\his_ensmean.nc','NOWRITE');
ncdisp('H:\CESM-post\LE\Temperature\his_ensmean.nc');
TemMMM = ncread('H:\CESM-post\LE\Temperature\his_ensmean.nc','TEMP'); 
Temsub = permute(nansum(TemMMM(:,:,1,:)*5+TemMMM(:,:,2:37,:).*permute(dweit,[3 2 1 4]),3)/707,[1 2 4 3]);
%%
Tsubr1 = cat(1,Temsub(181:301,:,:,:),Temsub(302:360,:,:,:),Temsub(1:180,:,:,:));
[spacz spac] = areamean(Tsubr1,1:121,27:42,latData);
[siaz sia] = areamean(Tsubr1,122:360,27:42,latData);
[r1,p1,n_eff1] = corr_eff(spacz,siaz,0.1)  %  90% confidence
DIex = spac-sia;
close all;
plot(spac,'r');
hold on
plot(sia,'b');
% plot(DIex,'k')
clear par11 h01
index = DIex; 
Tsubr11 = permute(Tsubr1,[3 1 2]);
for i = 1:size(Tsubr1,1);
    for j = 1:size(Tsubr1,2);
        [par11(i,j),h01(i,j),t] = reg1_ttest(index,Tsubr11(:,i,j),0.05,1);
    end
end

par11r = reshape(par11,[360,90]);
%%
close all;
    map = -par11r;
    Fig = figure('position',[100 100 600 300]);
    ticks = 2;
    contourfSPolar(map,[-ticks*10:ticks/10:ticks*10],ticks,lonData,latData(1:90),12)
    set_colorbar([0.83 0.08 0.03 0.88],[],4.5,12,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
    m_line([0:1:360],-45,'linewidth',2,'color','k');
    m_line([0:1:360],-60,'linewidth',2,'color','k');
    m_line(-60,[-60:1:-45],'linewidth',2,'color','k');
    m_line(150,[-60:1:-45],'linewidth',2,'color','k');




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
