% Data 
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
load("/Volumes/Togo4T/1_matlab/MatData/CESM/lonData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/latData.mat");
load("/Volumes/Togo4T/1_matlab/MatData/CESM/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
% LE ensmean (EX)
addpath /Volumes/Togo4T/1_matlab/help;
filename1 = '/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_1920-2005.nc';
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
TemMMM1 = ncread(filename1,'TEMP');
depthData = ncread(filename1,'z_t')/100; % cm -> m
filename2 = '/Volumes/CESM-post2/CESM-post/LE/TEMP/TEMPensmean_2006-2100.nc';
ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
TemMMM2 = ncread(filename2,'TEMP');
TemMMMall = cat(4,TemMMM1(:,:,1:47,:),TemMMM2(:,:,1:47,:));
Temp = TemMMMall;
clear TemMMM1 TemMMM2 TemMMMall
%%  Tem
lats=1:90;
Tsub_r = cat(1,Temp(152:360,lats,:,:),Temp(1:151,lats,:,:));
lonData_r = [lonData(152:360);lonData(1:151)+360];
Tpac = squeeze(nanmean(Tsub_r(1:140,:,:,:),1));
Taio = squeeze(nanmean(Tsub_r(141:360,:,:,:),1));
%% Tem trend 2d
startyr = 1960;
endyr = 2100;
var1 = permute(Tpac(:,:,startyr-1919:endyr-1919),[3 1 2]);
var2 = permute(Taio(:,:,startyr-1919:endyr-1919),[3 1 2]);
x = [1:size(var1,1)]';
clear trd1
for i = 1:size(var1,2);
    for j = 1:size(var1,3);
        par1=polyfit(x,var1(:,i,j),1); % regression parameters
        trd1(i,j) = par1(1); 
        par2=polyfit(x,var2(:,i,j),1); % regression parameters
        trd2(i,j) = par2(1); 
    end
end

max(trd1,[],'all')
min(trd1,[],'all')
%%
close all
ftsz = 20; ticks=0.03;
Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:37),trd1(:,1:37)',[-ticks*13:ticks/13:ticks*13],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([-35 -35],[0 -700],'linewidth',1.5,'color','k','linestyle','--')
line([-55 -55],[0 -700],'linewidth',1.5,'color','k','linestyle','--')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
colormap(flipud(bl_re4));
set(gca,'XLim',[-80,-30]);
set(gca,'XTick',[-80:10:-30]);
set(gca,'XTickLabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',ftsz);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',ftsz);
set(gca,'YTickLabel',[700:-100:0],'FontSize',ftsz);
ylabel('Depth (m)')
print(Fig,['/Users/yysong/Desktop/Inter-basin contrast/20241022/schematic_Pac_1960-2100.png'],'-dpng','-r300')

Fig = figure('position',[100 100 800 400]);
contourf(latData(lats),-depthData(1:37),trd2(:,1:37)',[-ticks*13:ticks/13:ticks*13],'linestyle','none')
caxis([-ticks,ticks])
hold on
line([-35 -35],[0 -700],'linewidth',1.5,'color','k','linestyle','--')
line([-55 -55],[0 -700],'linewidth',1.5,'color','k','linestyle','--')
% [ac ah] = contour(lonData,-depthData,(nanmean(latmc,3))',14,'ShowText','on','LevelList',[0:2:30],'LineColor','k','LabelSpacing',240);
load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
%bl_re(18:4:26,:) = []; bl_re(1:4:9,:) = [];
colormap(flipud(bl_re4));

set(gca,'XLim',[-80,-30]);
set(gca,'XTick',[-80:10:-30]);
set(gca,'XTickLabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',ftsz);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',ftsz);
set(gca,'YTickLabel',[700:-100:0],'FontSize',ftsz);
ylabel('Depth (m)')
print(Fig,['/Users/yysong/Desktop/Inter-basin contrast/20241022/schematic_AIO_1960-2100.png'],'-dpng','-r300')
%% Tem 2080-2100 minus 1960-1980
map_pac = nanmean(Tpac(:,:,2080-1919:2100-1919),3)-nanmean(Tpac(:,:,1960-1919:1980-1919),3);
map_aio = nanmean(Taio(:,:,2080-1919:2100-1919),3)-nanmean(Taio(:,:,1960-1919:1980-1919),3);
redyellow = [
    0.3, 0, 0;         % 深酒红色
    0.4, 0, 0.05;      % 深红色
    0.6, 0, 0.1;       % 稍亮的红色
    0.8, 0.1, 0.1;     % 更亮的红色
    1, 0.2, 0;         % 鲜红色
    1, 0.4, 0.1;       % 明显的橙红色
    1, 0.5, 0.2;       % 橙色
    1, 0.6, 0.3;       % 浅橙色
    1, 0.7, 0.4;       % 更浅的橙色
    1, 0.75, 0.45;     % 浅橙黄色
    1, 0.8, 0.5;       % 浅黄色
    1, 0.85, 0.6;      % 更浅的黄色
    1, 0.9, 0.7;       % 非常浅的黄色
    1, 0.95, 0.8;      % 最浅黄色
    1, 0.97, 0.85;     % 极浅黄色（接近黄色，但没有白色）
    1, 0.99, 0.9;      % 更极浅黄色
    1, 1, 0.95;        % 极淡黄色（没有白色）
    1, 1, 0.98;        % 额外添加的非常浅黄色（新色块）
];
close all
ftsz = 20; ticks=4.0;
% Pacific
Fig = figure('position',[100 100 800 600]);
contourf(latData(1:60),-depthData(1:37),fliplr(map_pac(1:60,1:37)'),[0:ticks/18:ticks*15],'linestyle','none')
caxis([0,ticks])
hold on
% 绘制指定等值线
% contour_levels = [3.0, 2.4, 1.8, 1.2, 0.6];
% [contour_lines, contour_handle] = contour(latData(lats), -depthData(1:37), map_pac(:,1:37)',...
%     contour_levels, 'LineColor', 'k','ShowText','on','LabelSpacing',240);

line([-65 -65],[0 -700],'linewidth',1.5,'color','k','linestyle','--')
line([-85 -85],[0 -700],'linewidth',1.5,'color','k','linestyle','--')

colormap(flipud(redyellow))

% ch = colorbar('location','eastoutside');
% set(ch,'Ticks',[-ticks:ticks/18:ticks]);
set(gca,'XLim',[-90,-30],'LineWidth',2);
set(gca,'XTick',[-90:10:-30]);
set(gca,'XTickLabel',{'30^oS','40^oS','50^oS','60^oS','70^oS','80^oS','90^oS'},'FontSize',ftsz);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',ftsz);
set(gca,'YTickLabel',[700:-100:0],'FontSize',ftsz);
ylabel('Depth (m)')
print(Fig,['/Users/yysong/Desktop/Inter-basin contrast/20241022/schematic_Pac_1980_20yrs_minus1960_20yrs_3.png'],'-dpng','-r300')

%% Atlantic-Indian
Fig = figure('position',[100 100 800 600]);
contourf(latData(lats),-depthData(1:37),map_aio(:,1:37)',[0:ticks/18:ticks*15],'linestyle','none')
caxis([0,ticks])
hold on

% 绘制指定等值线
% contour_levels = [ 3, 2.4, 1.8, 1.2];
% [contour_lines, contour_handle] = contour(latData(lats), -depthData(1:37), map_aio(:,1:37)',...
%     contour_levels, 'LineColor', 'k','ShowText','on','LabelSpacing',240);

line([-35 -35],[0 -700],'linewidth',1.5,'color','k','linestyle','--')
line([-55 -55],[0 -700],'linewidth',1.5,'color','k','linestyle','--')
colormap(flipud(redyellow))

 % ch = colorbar('location','eastoutside');
 % set(ch,'Ticks',[0:ticks/18:ticks]);

set(gca,'XLim',[-90,-30],'LineWidth',2);
set(gca,'XTick',[-90:10:-30]);
set(gca,'XTickLabel',{'90^oS','80^oS','70^oS','60^oS','50^oS','40^oS','30^oS'},'FontSize',ftsz);
set(gca,'YLim',[-700,0]);
set(gca,'YTick',[-700:100:0],'FontSize',ftsz);
set(gca,'YTickLabel',[700:-100:0],'FontSize',ftsz);
ylabel('Depth (m)')
print(Fig,['/Users/yysong/Desktop/Inter-basin contrast/20241022/schematic_AIO_1980_20yrs_minus1960_20yrs_3.png'],'-dpng','-r300')
%%  Tem 
% ------------------------- 0-700m ----------------------------------------
depthstr = '0-700';
dweit = depthData(2:37)-depthData(1:36);
lats = 1:60;
TemMMMsub = permute(nansum(cat(3,Temp(:,:,1,:)*5,Temp(:,:,2:37,:).*permute(dweit,[3 2 1 4])),3)/707,[1 2 4 3]);

%%
TemMMMsub=squeeze(Temp(:,:,1,:));
Tsub_r = cat(1,TemMMMsub(182:360,lats,:),TemMMMsub(1:181,lats,:));
lonData_r = [lonData(182:360);lonData(1:181)+360];
% Tem plot
map=nanmean(Tsub_r(:,:,2080-1919:2100-1919),3)-nanmean(Tsub_r(:,:,1960-1919:1980-1919),3);
close all;
ftsz = 12; ticks = 3;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on','LineWidth',2);
contourfS(map,[0:ticks/18:ticks*5],ticks,lonData_r,latData(lats),12)
%print(Fig,['/Users/yysong/Desktop/Inter-basin contrast/20241022/Tem.png'],'-dpng','-r300')

function contourfS(sic,val,ctr,lonData,latData,ftsz)
    m_proj('stereographic','lon',-30,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData,latData,sic',val,'linestyle','none');
    caxis([0,ctr]);
    redyellow = [
    0.3, 0, 0;         % 深酒红色
    0.4, 0, 0.05;      % 深红色
    0.6, 0, 0.1;       % 稍亮的红色
    0.8, 0.1, 0.1;     % 更亮的红色
    1, 0.2, 0;         % 鲜红色
    1, 0.4, 0.1;       % 明显的橙红色
    1, 0.5, 0.2;       % 橙色
    1, 0.6, 0.3;       % 浅橙色
    1, 0.7, 0.4;       % 更浅的橙色
    1, 0.75, 0.45;     % 浅橙黄色
    1, 0.8, 0.5;       % 浅黄色
    1, 0.85, 0.6;      % 更浅的黄色
    1, 0.9, 0.7;       % 非常浅的黄色
    1, 0.95, 0.8;      % 最浅黄色
    1, 0.97, 0.85;     % 极浅黄色（接近黄色，但没有白色）
    1, 0.99, 0.9;      % 更极浅黄色
    1, 1, 0.95;        % 极淡黄色（没有白色）
    1, 1, 0.98;        % 额外添加的非常浅黄色（新色块）
];

    colormap(flipud(redyellow))
    m_coast('patch','w');
    m_coast('linewidth',2,'color','k');
    m_grid('fontsize',ftsz,'xtick',[],'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','none','xaxisloc','top');
end
