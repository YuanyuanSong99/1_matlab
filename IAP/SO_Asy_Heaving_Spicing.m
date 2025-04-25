load('G:/1_matlab/IAP/MatFile/heaving_T_durack.mat','heaving_T_durack');
load('G:/1_matlab/IAP/MatFile/spicing_T_durack.mat','spicing_T_durack');
load('G:/1_matlab/IAP/MatFile/lonData.mat','lonData');
load('G:/1_matlab/IAP/MatFile/latData.mat','latData');
load('G:/1_matlab/IAP/MatFile/depthData.mat','depthData');
%% OHC
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
%------------------------- 0-700m -----------------------------------------
% spiT = permute(spicing_T_durack,[2 1 3 4]);
% hevT = permute(heaving_T_durack,[2 1 3 4]);
% depthstr = '0-700';
% dweit = depthData(2:27)-depthData(1:26); % depth weight
%------------------------- 0-2000m -----------------------------------------
spiT = permute(spicing_T_durack,[2 1 3 4]);
hevT = permute(heaving_T_durack,[2 1 3 4]);
depthstr = '0-2000';
dweit = depthData(2:41)-depthData(1:40); % depth weight

spiTsub = cp*ro*permute(nansum(cat(3,spiT(:,:,1,:),spiT(:,:,2:41,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]);
hevTsub = cp*ro*permute(nansum(cat(3,hevT(:,:,1,:),hevT(:,:,2:41,:).*(permute(dweit,[3 2 1]))),3),[1 2 4 3]);

yrspan = 20;
yrs1 = 1960; 
yrs2 = 2000;

spimap = nanmean(spiTsub(:,:,yrs2-1939:yrs2-1939+yrspan-1),3)-nanmean(spiTsub(:,:,yrs1-1939:yrs1-1939+yrspan-1),3);
hevmap = nanmean(hevTsub(:,:,yrs2-1939:yrs2-1939+yrspan-1),3)-nanmean(hevTsub(:,:,yrs1-1939:yrs1-1939+yrspan-1),3);
max(spimap,[],'all')
min(spimap,[],'all')
%% test with Duan et al. 2022 (2000-2014 - 1958-1972)
close all;
Fig = figure('Position',[300 300 600 400])
contourVARra((spimap)/10^9,[-4:0.1:4],2,0,160,-60,60,lonData,latData);colorbar;
title('spicing Indian Ocean')
print(Fig,['G:\figures\IAP\Yearly\20231104_SO_warming\DJ_spicing.png'],'-dpng','-r300')

%% SO asymmetry
close all;
ftsz = 12; ticks = 2;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
contourfS((spimap)/10^9,[-ticks*5:0.01:ticks*5],ticks,lonData,latData,12)
m_line([0:1:360],-35,'linewidth',2,'color','k'); 
% m_line([0:1:360],-50,'linewidth',2,'color','k'); 
m_line(-60,[-90:1:-35],'linewidth',2,'color','k'); 
m_line(180,[-90:1:-35],'linewidth',2,'color','k'); 
set_colorbar([0.83 0.08 0.03 0.88],'10^9 J',5,ftsz,[-ticks:ticks/5:ticks],[-ticks:ticks/5:ticks])
% testdots(h0,[.4 .4 .4],lonData,latData);
print(Fig,['G:\figures\IAP\Yearly\20231107_SO_warming\Delta_OHC_spi',depthstr,'m_',num2str(yrs1),'_',num2str(yrs2),'_20yrs.png'],'-dpng','-r300')



function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',4,'yticklabel',[],'box','on','linestyle','--','xaxisloc','top');
end
function contourVARra(var,val,ctr,lon1,lon2,lat1,lat2,lonData,latData)
    m_proj('equidistant Cylindrical','lat',[lat1 lat2],'long',[lon1 lon2]);
    hold on
    m_contourf(lonData,latData,var',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('G:/1_matlab/help/colorbar_mat/bl_re4.mat');
    colormap(flipud(bl_re4));
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',7,'ytick',7,'box','on','linestyle','none');
end
function set_colorbar(cbr_po,strings,str_po,ftsz,ticks_val,tickslabel)
    hb1 = colorbar('location','westoutside');
    set(hb1,'Units','normalized','position',cbr_po,'fontsize',ftsz,'Ticks',ticks_val,'TickLabels',tickslabel);
    set(hb1.Label,'String',strings,'Units','normalized','position',[str_po 0.5 0],'Rotation',-90,'fontsize',ftsz);
end