inlats = 35:55;
adall = squeeze(nansum(nansum(nansum(advTV(300:359,inlats,:,:),1),2),3)); 
advall = 60*60*24*365*cumsum(adall(41:end))/10^21;
close all
Fig = figure('position',[700 100 800 400]);
% plot(tdall,'r')
% hold on
% plot(focall,'b')
plot(advall,'g')
% plot(focall(1:46)+advall','k')
set(gca,'XLim',[1,61],'ylim',[-90,90],'ytick',[-90:30:90]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
ylabel('ZJ')
% legend('tendency','forcing','transport','forcing+transport','location','north')
% legend('boxoff')
%% adv
test1 = squeeze(nanmean(nanmean(advTV(:,:,:,:).*permute(dweit,[3 2 1]),3),4)); 
% close all;
ftsz = 12; ticks = .001;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(test1/10^15,[-ticks*5:ticks/10:ticks*5],ticks,lonData(2:360),latData(2:90),12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-55,'linewidth',2,'color','k'); 
m_line(-60,[-55:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-55:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'^oC s^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
%% dTdx
test1 = squeeze(nanmean(nanmean(URES(2:end,2:splat,2:inlev,:).*dTdx(:,:,:,:).*permute(dweit,[3 2 1]),3),4)); 
max(test1,[],'all')
min(test1,[],'all')
% close all;
ftsz = 12; ticks = .01/10^4;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(test1,[-ticks*5:ticks/10:ticks*5],ticks,lonData(2:360),latData(2:90),12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-55,'linewidth',2,'color','k'); 
m_line(-60,[-55:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-55:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'^oC m^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
%%
test1 = squeeze(nanmean(nanmean(VRES(2:end,2:splat,2:inlev,:).*dTdy(:,:,:,:).*permute(dweit,[3 2 1]),3),4)); 
max(test1,[],'all')
min(test1,[],'all')
% close all;
ftsz = 12; ticks = .01/10^4;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(test1,[-ticks*5:ticks/10:ticks*5],ticks,lonData(2:360),latData(2:90),12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-55,'linewidth',2,'color','k'); 
m_line(-60,[-55:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-55:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'^oC m^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
%%
test1 = squeeze(nanmean(nanmean(WRES(2:end,2:splat,2:inlev,:).*dTdz(:,:,:,:).*permute(dweit,[3 2 1]),3),4)); 
max(test1,[],'all')
min(test1,[],'all')
% close all;
ftsz = 12; ticks = .01/10^4;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(test1,[-ticks*5:ticks/10:ticks*5],ticks,lonData(2:360),latData(2:90),12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-55,'linewidth',2,'color','k'); 
m_line(-60,[-55:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-55:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'^oC m^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
%%
test1 = squeeze(nanmean(nanmean((VRES(2:end,2:splat,2:inlev,:).*dTdy(:,:,:,:)+URES(2:end,2:splat,2:inlev,:).*dTdx(:,:,:,:)).*permute(dweit,[3 2 1]),3),4)); 
max(test1,[],'all')
min(test1,[],'all')
% close all;
ftsz = 12; ticks = .01/10^4;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS(test1,[-ticks*5:ticks/10:ticks*5],ticks,lonData(2:360),latData(2:90),12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
m_line([0:1:360],-55,'linewidth',2,'color','k'); 
m_line(-60,[-55:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-55:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'^oC m^-^1',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])

function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end