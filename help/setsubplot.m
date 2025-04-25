close all;
figure(1)
subplot('position',[0.1 0.74 0.35 0.16]); % left bottom width height
contourf(duramd2(49:133,:)',[-250:20:250]);
hold on
contour(duramd2(49:133,:)',[80 80],'k','LineWidth',2); %�ر���ע80gpm
colormap(othercolor('BuOr_12'));
%colorbar;
caxis([-250,250]); %����ɫ�귶Χ���ǳ���������
set(gca,'XTick',[1:12:90]);
set(gca,'XTickLabel',{'60^oW','30^oW','0^o','30^oE','60^oE','90^oE','120^oE','150^oE'},'FontSize',8);
set(gca,'YTick',[1:2:dzn]);
set(gca,'YTickLabel',[-(dzn-1)/2:2:(dzn-1)/2],'FontSize',8);
xlabel('Longitude','FontSize',8),ylabel('Lag(day)','FontSize',8);
title('PVy(U) > 0.5','FontSize',12);

subplot('position',[0.5 0.74 0.35 0.16]);
contourf(duramx2(49:133,:)',[-200:20:250]);
hold on
contour(duramx2(49:133,:)',[80 80],'k','LineWidth',2);
colormap(othercolor('BuOr_12'));
%colorbar;
caxis([-200,250]);
set(gca,'XTick',[1:12:90]);
set(gca,'XTickLabel',{'60^oW','30^oW','0^o','30^oE','60^oE','90^oE','120^oE','150^oE'},'FontSize',8);
set(gca,'YTick',[1:2:dzn]);
set(gca,'YTickLabel',[-(dzn-1)/2:2:(dzn-1)/2],'FontSize',8);
xlabel('Longitude','FontSize',8),ylabel('Lag(day)','FontSize',8);
title('PVy(U) < -0.5','FontSize',12);

subplot('position',[0.1 0.52 0.35 0.16]);
contourf(duramd3(49:133,:)',[-200:20:250]);
hold on
contour(duramd3(49:133,:)',[80 80],'k','LineWidth',2);
colormap(othercolor('BuOr_12'));
%colorbar;
caxis([-200,250]);
set(gca,'XTick',[1:12:90]);
set(gca,'XTickLabel',{'60^oW','30^oW','0^o','30^oE','60^oE','90^oE','120^oE','150^oE'},'FontSize',8);
set(gca,'YTick',[1:2:dzn]);
set(gca,'YTickLabel',[-(dzn-1)/2:2:(dzn-1)/2],'FontSize',8);
xlabel('Longitude','FontSize',8),ylabel('Lag(day)','FontSize',8);
title('U > 0.5','FontSize',12);
%
subplot('position',[0.5 0.52 0.35 0.16]);
contourf(duramx3(49:133,:)',[-200:20:250]);
hold on
contour(duramx3(49:133,:)',[80 80],'k','LineWidth',2);
colormap(othercolor('BuOr_12'));
%colorbar;
caxis([-200,250]);
set(gca,'XTick',[1:12:90]);
set(gca,'XTickLabel',{'60^oW','30^oW','0^o','30^oE','60^oE','90^oE','120^oE','150^oE'},'FontSize',8);
set(gca,'YTick',[1:2:dzn]);
set(gca,'YTickLabel',[-(dzn-1)/2:2:(dzn-1)/2],'FontSize',8);
xlabel('Longitude','FontSize',8),ylabel('Lag(day)','FontSize',8);
title('U < -0.5','FontSize',12);

subplot('position',[0.5 0.08 0.35 0.16]);
contourf(duramx4(49:133,:)',[-200:20:250]);
hold on
contour(duramx4(49:133,:)',[80 80],'k','LineWidth',2);
colormap(othercolor('BuOr_12'));
%colorbar;
caxis([-200,250]);
set(gca,'XTick',[1:12:90]);
set(gca,'XTickLabel',{'60^oW','30^oW','0^o','30^oE','60^oE','90^oE','120^oE','150^oE'},'FontSize',8);
set(gca,'YTick',[1:2:dzn]);
set(gca,'YTickLabel',[-(dzn-1)/2:2:(dzn-1)/2],'FontSize',8);
xlabel('Longitude','FontSize',8),ylabel('Lag(day)','FontSize',8);
title('-Ty < -0.5','FontSize',12);
%
subplot('position',[0.1 0.08 0.35 0.16]);
contourf(duramd4(49:133,:)',[-200:20:250]);
hold on
contour(duramd4(49:133,:)',[80 80],'k','LineWidth',2);
colormap(othercolor('BuOr_12'));
%colorbar;
caxis([-200,250]);
set(gca,'XTick',[1:12:90]);
set(gca,'XTickLabel',{'60^oW','30^oW','0^o','30^oE','60^oE','90^oE','120^oE','150^oE'},'FontSize',8);
set(gca,'YTick',[1:2:dzn]);
set(gca,'YTickLabel',[-(dzn-1)/2:2:(dzn-1)/2],'FontSize',8);
xlabel('Longitude','FontSize',8),ylabel('Lag(day)','FontSize',8);
title('-Ty > 0.5','FontSize',12);

subplot('position',[0.5 0.30 0.35 0.16]);
contourf(duramx5(49:133,:)',[-200:20:250]);
hold on
contour(duramx5(49:133,:)',[80 80],'k','LineWidth',2);
colormap(othercolor('BuOr_12'));
%colorbar;
caxis([-200,250]);
set(gca,'XTick',[1:12:90]);
set(gca,'XTickLabel',{'60^oW','30^oW','0^o','30^oE','60^oE','90^oE','120^oE','150^oE'},'FontSize',8);
set(gca,'YTick',[1:2:dzn]);
set(gca,'YTickLabel',[-(dzn-1)/2:2:(dzn-1)/2],'FontSize',8);
xlabel('Longitude','FontSize',8),ylabel('Lag(day)','FontSize',8);
title('-Uyy < -0.5','FontSize',12);
%
subplot('position',[0.1 0.30 0.35 0.16]);
contourf(duramd5(49:133,:)',[-200:20:250]);
hold on
contour(duramd5(49:133,:)',[80 80],'k','LineWidth',2);
colormap(othercolor('BuOr_12'));
%colorbar;
caxis([-250,250]);
set(gca,'XTick',[1:12:90]);
set(gca,'XTickLabel',{'60^oW','30^oW','0^o','30^oE','60^oE','90^oE','120^oE','150^oE'},'FontSize',8);
set(gca,'YTick',[1:2:dzn]);
set(gca,'YTickLabel',[-(dzn-1)/2:2:(dzn-1)/2],'FontSize',8);
xlabel('Longitude','FontSize',8),ylabel('Lag(day)','FontSize',8);
title('-Uyy > 0.5','FontSize',12);
% set colorbar
hb = colorbar('location','southoutside');
%set(hb,'Units','normalized','position',[0.7 0.06 0.01 0.85]);
set(hb,'Units','normalized','position',[0.1 0.02 0.75 0.02]); % left bottom width height
sfig=figure(1);
set(sfig,'PaperPositionMode','manual');
set(sfig,'PaperUnits','inches');
set(sfig,'PaperPosition',[0 0 15 20]);
set(sfig,'PaperSize',[15,20]);
saveas(gcf,['D:\sea ice\figures\whole_t\PVy(U)\cos\duration\rawduration',num2str(latData(lats(end))),'_',num2str(latData(lats(1)))...
    '&',num2str(lonData(lons(1))),'_',num2str(lonData(lons(end))),'.png']);