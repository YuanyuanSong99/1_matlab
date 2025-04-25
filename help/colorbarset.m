%% set colorbar position
% after contour
hb = colorbar('location','southoutside');
set(hb,'Units','normalized','position',[0.29 0.02 0.38 0.02]); % left bottom width height
%% 看到喜欢的colorbar就可以截图，读进来用。
close all;
imread('G:/blre4.png');
colortemp = ans(1:38:end,:,:);
color=colortemp(:,1,:);
colorfinal=squeeze(color);
colorbar
colormap(double(colorfinal)/255) %需要转化成双精度，0-1之间的数值

bl_re4=double(unique(colorfinal,'rows','stable'))/255; % 偶数个色块
colorbar
colormap(double(bl_re4))
save('G:\1_matlab\help\colorbar_mat\bl_re4.mat','bl_re4');
%% 将colorbar具体分成10，12，14，16，18，20个色块，使数值精确与色块对应
close all;
% 如有新colorbar，(1)修改bl_re名称; (2)修改数值范围确保输出的colorbar的length为偶数
% bl_re4.mat
load('G:/1_matlab/help/colorbar_mat/bl_re4.mat');
load('G:/1_matlab/help/par0.mat');
length(bl_re4)
% 10-01
bl_re4_10_01 = colorbar10(bl_re4);
length(bl_re4_10_01)
colorbar
colormap(double(bl_re4_10_01))
% save('E:\1_matRlab\help\colorbar_mat\bl_re4_10_01.mat','bl_re4_10_01');
% showmap(map,0.5,bl_re4_10_01,[10 50 550 400],lonData,latData)
%% 10-02
rmv = ceil(length(bl_re4)/6);
bl_re4_10_02 = colorbar10(bl_re4(rmv:length(bl_re4)-rmv,:));
length(bl_re4_10_02)
save('E:\1_matlab\help\colorbar_mat\bl_re4_10_02.mat','bl_re4_10_02');
showmap(-par0,0.5,bl_re4_10_02,[600 50 550 400])
% 12-01
bl_re4_12_01 = colorbar12(bl_re4);
length(bl_re4_12_01)
save('E:\1_matlab\help\colorbar_mat\bl_re4_12_01.mat','bl_re4_12_01');
showmap(-par0,0.6,bl_re4_12_01,[10 50 550 400])
% 12-02
rmv = ceil(length(bl_re4)/12);
bl_re4_12_02 = colorbar12(bl_re4(rmv:length(bl_re4)-rmv,:));
length(bl_re4_12_02)
save('E:\1_matlab\help\colorbar_mat\bl_re4_12_02.mat','bl_re4_12_02');
showmap(-par0,0.6,bl_re4_12_02,[600 50 550 400])
% 14-01
bl_re4_14_01 = colorbar14(bl_re4);
length(bl_re4_14_01)
save('E:\1_matlab\help\colorbar_mat\bl_re4_14_01.mat','bl_re4_14_01');
showmap(-par0,0.7,bl_re4_14_01,[10 50 550 400])
% 14-02
rmv = ceil(length(bl_re4)/15);
bl_re4_14_02 = colorbar14(bl_re4(rmv:length(bl_re4)-rmv,:));
length(bl_re4_14_02)
save('E:\1_matlab\help\colorbar_mat\bl_re4_14_02.mat','bl_re4_14_02');
showmap(-par0,0.7,bl_re4_14_02,[600 50 550 400])
% 16-01
bl_re4_16_01 = colorbar16(bl_re4);
length(bl_re4_16_01)
save('E:\1_matlab\help\colorbar_mat\bl_re4_16_01.mat','bl_re4_16_01');
showmap(-par0,0.8,bl_re4_16_01,[10 50 550 400])
% 16-02
rmv = ceil(length(bl_re4)/15);
bl_re4_16_02 = colorbar16(bl_re4(rmv:length(bl_re4)-rmv,:));
length(bl_re4_16_02)
save('E:\1_matlab\help\colorbar_mat\bl_re4_16_02.mat','bl_re4_16_02');
showmap(-par0,0.8,bl_re4_16_02,[600 50 550 400])
% 18-01
bl_re4_18_01 = colorbar18(bl_re4);
length(bl_re4_18_01)
save('E:\1_matlab\help\colorbar_mat\bl_re4_18_01.mat','bl_re4_18_01');
showmap(-par0,0.9,bl_re4_18_01,[10 50 550 400])
% 18-02
rmv = ceil(length(bl_re4)/31);
bl_re4_18_02 = colorbar18(bl_re4(rmv:length(bl_re4)-rmv,:));
length(bl_re4_18_02)
save('E:\1_matlab\help\colorbar_mat\bl_re4_18_02.mat','bl_re4_18_02');
showmap(-par0,0.9,bl_re4_18_02,[600 50 550 400])






function [color_new] = colorbar10(color_org)
    color_new = double(color_org([1:ceil(length(color_org)/9):end],:));
    length(color_new)
    color_new(5:6,:) = ones(2,3);
end
function [color_new] = colorbar12(color_org)
    color_new = double(color_org([1:floor(length(color_org)/11):end],:));
    color_new(6:7,:) = ones(2,3);
end
function [color_new] = colorbar14(color_org)
    color_new = double(color_org([1:floor(length(color_org)/13):end],:));
    color_new(7:8,:) = ones(2,3);
end
function [color_new] = colorbar16(color_org)
    color_new = double(color_org([1:floor(length(color_org)/15):end],:));
    color_new(8:9,:) = ones(2,3);
end
function [color_new] = colorbar18(color_org)
    color_new = double(color_org([1:floor(length(color_org)/17):end],:));
    color_new(9:10,:) = ones(2,3);
end
function showmap(map,ctl,xxx,pos,lonData,latData)
ftsz = 14;
Fig = figure('position',pos);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
m_proj('stereographic','lon',90,'lat',90,'radius',70,'rotangle',0);
hold on
m_contourf(lonData,latData,map',[-0.8:0.1:0.8],'linestyle','none');
colormap(xxx)
caxis([-ctl,ctl]);
colorbar
m_coast('linewidth',1,'color','k');
m_grid('xtick',13,'ytick',0,'yticklabel',[],'box','on','linestyle','--');
end


