close all;
ftsz = 12; ticks = 0.25;
Fig = figure('position',[10 50 550 400]);
ax = axes('Position',[0.005 0.07 0.9 0.87],'fontsize',ftsz,'box','on');
sic=trd*10; val=[-ticks*5:ticks/10:ticks*5]; ctr=ticks; 
    m_proj('stereographic','lon',0,'lat',-90,'radius',60,'rotangle',0);
    hold on
    m_contourf(lonData_r,latData(lats),sic',val,'linestyle','none');
    caxis([-ctr,ctr]);
    load('/Volumes/Togo4T/1_matlab/help/colorbar_mat/bl_re4.mat');
    bl_re4(18:4:26,:) = []; bl_re4(1:4:9,:) = []; 
    colormap(flipud(bl_re4));
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('fontsize',ftsz,'xtick',13,'tickdir','out','ytick',2,'yticklabel',[],'box','on','linestyle','none','xaxisloc','top');

view(0, 30);
%%
% 清理工作区并关闭所有图形窗口
clear; close all; clc;

% 读取并显示二维图像
img = imread('/Users/yysong/Desktop/Inter-basin contrast/20241022/Tem.png'); % 替换为你想使用的图像文件
figure;
imshow(img);
title('原始二维图像');

% 在三维空间中创建一个平面并将图像映射到平面上
figure;
[X, Y] = meshgrid(1:size(img, 2), 1:size(img, 1)); % 使用图像的像素尺寸创建网格
Z = zeros(size(X)); % 创建一个平面Z坐标为0

% 将图像纹理贴到三维平面
h = surface(X, Y, Z, flipud(img), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
axis equal;
hold on;

% 绕 X 轴旋转图像，使其倾斜
rotate(h, [1, 0, 0], 90); % [1, 0, 0] 表示绕 X 轴旋转，60 表示旋转角度（可调）

% 设置视角，仰角和方位角调整
view(0, 60); % 可以根据需求调整

% 设置图形属性
xlabel('X');
ylabel('Y');
zlabel('Z');
title('静态视角：图像从屏幕外向屏幕内翻转');
axis tight;
grid on;
print(figure(2),['/Users/yysong/Desktop/Inter-basin contrast/20241022/Tem_1.png'],'-dpng','-r300')


%%
% 获取图像尺寸
[img_height, img_width, ~] = size(img);

% 创建三维平面，将图像映射到平面上
figure;
[X, Y] = meshgrid(1:img_width, 1:img_height); % 使用图像的像素尺寸创建网格
Z = zeros(size(X)); % 创建一个平面Z坐标为0

% 将图像纹理贴到三维平面
surface('XData', X, 'YData', Y, 'ZData', Z, 'CData', flipud(img), ...
        'FaceColor', 'texturemap', 'EdgeColor', 'none');
axis equal;
hold on;

% 绕 X 轴旋转平面，实现仰角效果
% 创建新的 Z 坐标数据，实现绕 X 轴旋转
angle = 60; % 旋转角度（可以根据需要调整）
Z_rotated = Y * sind(angle); % 将 Y 坐标按角度转换成 Z 坐标，实现 X 轴旋转效果
Y_rotated = Y * cosd(angle); % 更新 Y 坐标

% 更新图像所在平面的 Y 和 Z 数据，实现旋转后的效果
set(gca, 'YDir', 'reverse'); % 反转 Y 轴方向，使图像正常显示
surface('XData', X, 'YData', Y_rotated, 'ZData', Z_rotated, 'CData', flipud(img), ...
        'FaceColor', 'texturemap', 'EdgeColor', 'none');

% 设置视角，仰角和方位角调整
view(0, 30); % 设置观察角度，仰角 30 度，可根据需求调整

% 设置图形属性
xlabel('X');
ylabel('Y');
zlabel('Z');
title('静态视角：图像从屏幕外向屏幕内翻转');
axis tight;
grid on;


