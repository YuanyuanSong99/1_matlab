function new_fig_handle = convert_to_std_coordinate_system( fig_handle,gridflag)
% 本函数目的是把 matlab 做的图的坐标轴在原点十字带箭头显示(与数学的做图习惯一致)
% 2019.2.2 in ucas
% author:LuSong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    gridflag = 0;
end
figure('Name','std_coordinate_system','NumberTitle','off') % Create a new figure
new_fig_handle = copyobj( fig_handle , gcf );%拷贝图形到一个新的窗口
xL=xlim ;%取原来的x和y的范围
yL=ylim ;
xt=get(gca,'xtick') ;%取原来的tick
yt=get(gca,'ytick') ;
%set(gca,'YColor','w') ;
%set(gca,'XColor','w') ;
%set(gca,'xticklabel','');
%set(gca,'yticklabel','');
if gridflag == 0
set(gca,'xtick',[]);
set(gca,'ytick',[]);
else
grid on;
set(gca,'xticklabel','');
set(gca,'yticklabel','');
end

% 把 x 和 y 坐标轴的两个方向各延长 10% （为了视觉上好看）
extend_x = ( xL(2)-xL(1) ) * 0.1 ;
extend_y = ( yL(2)-yL(1) ) * 0.1 ;
xxL = xL + [ -extend_x extend_x] ;
yyL = yL + [ -extend_y extend_y] ;
set(gca,'xlim', xxL) ;%设置Range增加10%
set(gca,'ylim', yyL) ;
set(gca,'Position',[0.05 0.1 0.8 0.8])
pos = get(gca,'Position') ;%获取绘图区域，pos(1)指左下角点的x坐标，pos(2)指左下角点的y坐标,pos(3)指宽度，pos(4)指高度
box off;%不展示轴的边界
x_shift = abs( yyL(1)/(yyL(2)-yyL(1)) ) ;
y_shift = abs( xxL(1)/(xxL(2)-xxL(1)) ) ;
temp_1 = axes('Position', pos + [ 0 , pos(4) * x_shift , 0 , - pos(4)* x_shift ] ) ;%新建一个axes
xlim(xxL) ;
box off ;
xt(find(xt==0))=[];
set(temp_1,'XTick',xt,'fontsize',12,'Color','None','YTick',[]) ;%显示我想要的范围刻度
set(temp_1,'YColor','w') ;
temp_2 = axes( 'Position', pos + [ pos(3) * y_shift , 0 , -pos(3)*y_shift*0.99999 , 0 ] ) ;
ylim(yyL) ;
box off ;
yt(find(yt==0))=[];
set(temp_2,'YTick',yt,'fontsize',12,'Color','None','XTick',[]) ;
set(temp_2,'XColor','w') ;
text(-0.07,-0.2,'0','FontSize',12)
Base_pos = get(new_fig_handle,'Position') ;
arrow_pos_in_x_dircetion = Base_pos(2) - Base_pos(4) * yyL(1)/(yyL(2)-yyL(1)) ;
arrow_pos_in_y_dircetion = Base_pos(1) - Base_pos(3) * xxL(1)/(xxL(2)-xxL(1)) ;
annotation('arrow',[Base_pos(1), Base_pos(1)+Base_pos(3)] , [arrow_pos_in_x_dircetion , arrow_pos_in_x_dircetion ] , 'Color','k');
annotation('arrow',[arrow_pos_in_y_dircetion, arrow_pos_in_y_dircetion ] , [Base_pos(2) , Base_pos(2)+Base_pos(4)] , 'Color','k');
%下面的处理是为了处理边框
% temp_addition = axes( 'Position', pos) ;
% box on ;
% set(temp_addition,'Color','None') ;
% set(temp_addition,'xtick',[],'ytick',[]) ;
temp_addition = axes( 'Position', pos) 
box off ;
set(temp_addition,'Color','None') ;
set(temp_addition,'xtick',[],'ytick',[]) ;
set(temp_addition,'XColor','w') ;
set(temp_addition,'YColor','w')
box off
end