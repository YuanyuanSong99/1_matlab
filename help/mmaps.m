%% Equidistant
% From Yanan_t2m_trend_gw.m
fig11 = figure(1);
    set(fig11,'Unit','centimeters')  
    set(fig11,'Position',[0,0,33,18])  
    axes1 = axes('Parent',fig11,'Unit','centimeters','Position',[2 1 30 16]);
    datat=[p1(136:180,:,i,1);p1(1:136,:,i,1)]';
    m_proj('Equidistant','lon',rangelon,'lat',rangelat);
    m_contourf(lon1,lat1,datat*10,v1,'linestyle','none');
    m_coast('line','color',[0.5 0.5 0.5],'LineWidth',1.5);
    m_grid('xtick',13,'ytick',10,'tickdir','out','fontsize',16);
    caxis([-4 4])
    colormap(cc)
    colorbar('SouthOutside','position',[2/33 1.5/18 30/33 0.5/18],'fontsize',16)
    text(-1.8,2.4,strcat(strname(i),' surface temperature trends (1979¨C2017)'),'fontsize',24)
    research_rangee([40,60],[60,120],'-k',1.5)
    clear axes1
    set(fig11,'PaperUnits','centimeters','PaperPosition',[0 0 33 18])
    print(fig11,'-dpng','-r300',['C:\new\4EAcooling\0805_ttrend\trendtu\tad_trend19792017_',num2str(i),'.png']);
    clear fig11
    close(1)
%% stereographic 
% From Yanan_compose_ubuv.m    Including significance testing
axes5 = axes('Parent',fig11,'Unit','centimeters','Position',[2 2 10 10]);
datati=[p(:,:,1);p(1,:,1)]';
m_proj('stereographic','lon',60,'lat',90,'radius',60,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
hold on
datazi=[p(:,:,3);p(1,:,3)]';
dataz1=datazi;dataz2=datazi;
dataz1(dataz1<0)=0;dataz2(dataz2>0)=0;
m_contour(lon1,lat1,dataz1,'LevelList',v2,'LineStyle','-','color','k','LineWidth',1.2);
hold on
m_contour(lon1,lat1,dataz2,'LevelList',v3,'LineStyle','-.','color','k','LineWidth',1.2);
m_coast('line','color','k');
m_grid('xtick',12,'ytick',7,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

% From Davini2Dindex.m (set positions of title and  colorbar)
close all
fig11=figure(1);
axes5 = axes('Parent',fig11,'Units','normalized','position',[0.05 0.2 0.8 0.6]);
m_proj('stereographic','lon',0,'lat',90,'radius',60,'rotangle',0);
cc=1:0.5:25;
[C,h]=m_contourf(lon,lat,bl_freq'*100,cc,'linestyle','none');
hold on

m_coast('linewidth',1,'color','k');
m_grid('xtick',12,'ytick',4,'yticklabel',[],'tickdir','out','linest','--',...
    'tickdir','in','fontsize',10);
% set colorbar
hb = colorbar('location','southoutside');
%set(hb,'Units','normalized','position',[0.7 0.06 0.01 0.85]);
set(hb,'Units','normalized','position',[0.2 0.1 0.5 0.03]);
colormap(jet);
hb1 = title(' Large-scale Blocking Frequency (%)       ');
set(hb1,'Units','normalized','position',[0.55 1.1])
%% test dots
    dots = mod(find(h2 == 1),144);
    dots(find(dots == 0)) = 144;
    m_plot(lonData(dots),latData(ceil(find(h2 == 1)/144)),'.','markersize',4,'color','k'); % F-test dots
%% area box
function areabox(lon1,lon2,lat1,lat2)
    m_line(lon1,[lat1:1:lat2],'linewidth',2,'color','g'); 
    m_line(lon2,[lat1:1:lat2],'linewidth',2,'color','g');
    m_line([lon1:1:lon2],lat1,'linewidth',2,'color','g');
    m_line([lon1:1:lon2],lat2,'linewidth',2,'color','g');    
end
