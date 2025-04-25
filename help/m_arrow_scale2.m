%% M_ARROW_SCALE2使用说明
%                   给矢量场图添加箭头比例尺
%                        编程人：律成林（中国科学院大气物理研究所 博士研究生）
%==========================================================================
%    M_ARROW_SCALE2(HP, u, v, u0)
%    M_ARROW_SCALE2(HP, u, v, u0, Name, Value)
%    [arrow, txt, ax, layer] = M_ARROW_SCALE2(___)
%==========================================================================
%                                一、用法
%在使用m_vec函数绘制了矢量场图后（[HP, HT] = m_vec(s, lon, lat)），把输出
%的HP代入到m_arrow_scale2函数中绘制箭头比例尺：
%    M_ARROW_SCALE2(HP, u, v, u0)
%其中，u、v分别为为原始矢量场的横向分量和纵向分量，u0为箭头比例尺的指示值。
%在默认情况下，箭头比例尺位于该图右上角的新绘制的白框内，箭头在上，标注文本
%在下。用户也可以通过设定Name-Value值，来进行调整（详见第二节）：
%    M_ARROW_SCALE2(HP, u, v, u0, Name, Value)
%也可以通过直接设置对象的属性来进行调整（详见第三节）：
%    [arrow, txt, ax, layer] = M_ARROW_SCALE2(___)
%    set(arrow, Name, Value)
%==========================================================================
%                             二、Name-Value值
%例如设置标注文本，Name是'label'，Value是'5 m/s'，则输入的语句为：
%    M_ARROW_SCALE2(HP, u, v, u0, 'label', '5m/s')
%执行的效果是标注文本为 5m/s。
%----------------------- Name可以不区分大小写 -----------------------------
%Name            含义　　　　　Value
%'Label'         标注内容　　　如'5 m/s'等，char或cell类型（默认为u0对
%                　　　　　　　　　应的字符）
%'FontName'      标注字体　　　char或cell类型（默认为'Helvetica'）
%'FontSize'      标注字号　　　数字。单位：榜（默认根据图窗大小确定）
%'FontColor'     标注字颜色　　可为'r', 'b'等，也可为[R, G, B]向量（默
%                　　　　　　　　　认为黑色，'k' / [0, 0, 0]）
%'FontWeight'    标注是否加粗　'normal'不加粗（默认），'bold'加粗
%'FontAngle'     标注是否倾斜　'normal'不倾斜（默认），'italic'倾斜
%'Position'      框的位置　　　可用的Value有：
%                　'lu'/'left_up'/'lt'/'left_top'/'左上'：       图窗左上方
%                　'ru'/'right_up'/'rt'/'right_top'/'右上'：     图窗右上方
%                　'ld'/'left_down'/'lb'/'left_bottom'/'左下'：  图窗左下方
%                　'rd'/'right_down'/'rb'/'right_bottom'/'右下'：图窗右下方
%                　[x0, y0, xlength, ylength]：    框左下顶点坐标为[x0, y0]
%                　　　　　　　且长和高分别为xlength和ylength，在figure内的相
%                　　　　　　　对位置，范围为0～1
%'LabelPosition' 标注位置　　　默认标注在箭头下方。可用的Value有：
%                　'u'/'up'/'t'/'top'/'上'：            标注在箭头上方
%                　'd'/'down'/'b'/'bottom'/'下'：       标注在箭头下方
%                　[x, y]：   x, y是figure内的相对位置，分别从左、下开始计算，
%                　　　       范围为0～1。标注文本以该位置为中心坐标对齐
%'ArrowPosition' 箭头位置　　　箭头默认位于标注的上方。可用的Value有：
%                　'u'/'up'/'t'/'top'/'上'：            箭头在标注上方
%                　'd'/'down'/'b'/'bottom'/'下'：       箭头在标注下方
%                　[x, y]：   x, y是figure内的相对位置，分别从左、下开始计算，
%                　　　       范围为0～1。箭头以该位置为中心坐标对齐
%'BoxStyle'      框线线型　　　'-'实线（默认），'--'虚线，'-.'点划线，':'点线，
%                　　　　　　　　　'none'线条不可见
%'BoxColor'      框线颜色　　　同上面的'FontColor'的Value，此外还可为'none'来
%                　　　　　　　　　不显示框
%'BoxWidth'　　　框线宽度　　　实数。单位：榜（默认为0.5）
%'FaceColor'     填充颜色　　　同上面的'BoxColor'的Value（默认为白色，'w'）
%'FaceAlpha'     填充透明度　　0～1，0为完全透明，1为完全不透明（默认为1）
%==========================================================================
%                            三、输出对象
%    [arrow, txt, ax, layer] = M_ARROW_SCALE2(___)
%arrow：    在图窗ax内画出的箭头，类型为 1×1 Patch
%txt：      文本框，类型为 1×1 Text
%ax：       图窗，默认为透明无边框，坐标皆对应figure的相对位置，类型为 Axes
%layer：    装饰图层，最先绘制在ax内，类型为 1×1 Patch
%用户可通过 a = get(arrow) 查看arrow中的所有属性值，并用
%    set(arrow, Name, Value)
%来进行对对象属性值的设定。例如，设定箭头超出图窗的部分不显示：
%    set(ax, 'Clipping', 'on')
%对于txt、ax、layer同样如此。用户也可以在
%https://ww2.mathworks.cn/help/matlab/index.html
%网站上输入Patch并回车，点击检索到的第一项，然后在网页的最下方点击'Patch属性'
%超链接，查阅对象Patch所有的属性的含义和用法。(其他对象属性含义的查阅同上)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%完成日期：2022年1月22日
%see also: m_vec
 
function varargout = m_arrow_scale2(var_vec, u, v, u0, varargin)
%% 判断用户输入的参数格式是否正确
if nargin < 4
    warning('您输入的参数少于4个。m_arrow_scale的用法如下：')
    help m_arrow_scale2
    return
end
if nargin > 4
    if mod(length(varargin), 2) > 0
        error('您输入的属性中有值未输入。')
    end
end
%% 切换到矢量图所在的图窗
axes(var_vec.Parent);
%% 获取窗口和图窗的一些属性值
% 窗口位置
gcf_position = get(gcf, 'Position');
gca_position = get(gca, 'Position');
% 图窗长和高
gca_x_length = gcf_position(3) * gca_position(3);
gca_y_length = gcf_position(4) * gca_position(4);
% 图窗x和y的范围
gca_xlim = get(gca, 'XLim');
gca_ylim = get(gca, 'YLim');
%% 判断横纵轴单位长度的长是否相等
if strcmp(get(gca, 'DataAspectRatioMode'), 'manual')
    dx_dy_dz = get(gca, 'DataAspectRatio');
else
    dx_dy_dz = [gca_x_length / (gca_xlim(2) - gca_xlim(1)), ...
        gca_y_length / (gca_ylim(2) - gca_ylim(1)), 1];
end
if dx_dy_dz(1) ~= dx_dy_dz(2)
    warning(['请注意：图窗横纵轴单位长度的长并不相等，比例为 ', ...
        num2str(dx_dy_dz(1)), ' : ', num2str(dx_dy_dz(2)), ' 。'])
end
%% 获取图窗单位长度对应的像素数
% 图窗纵坐标和横坐标的单位长度的长之比
dy_by_dx = dx_dy_dz(2) / dx_dy_dz(1);
% 图窗的未缩放维度和缩放维度的判定（小于0，x未缩放；大于0，y未缩放）
zoom_dimension = gca_x_length / (gca_xlim(2) - gca_xlim(1)) * ...
    (gca_ylim(2) - gca_ylim(1)) * dy_by_dx - gca_y_length;
% 求出x单位长度对应的像素数
if zoom_dimension <= 0
    gca_x_length_real = gca_x_length;
    gca_y_length_real = gca_x_length / (gca_xlim(2) - gca_xlim(1)) * ...
        (gca_ylim(2) - gca_ylim(1)) * dy_by_dx;
    dx_real = gca_x_length / (gca_xlim(2) - gca_xlim(1));
else
    gca_x_length_real = gca_y_length / dy_by_dx / ...
        (gca_ylim(2) - gca_ylim(1)) * (gca_xlim(2) - gca_xlim(1));
    gca_y_length_real = gca_y_length;
    dx_real = gca_x_length_real / (gca_xlim(2) - gca_xlim(1));
end
% 窗口位置
gca_position = [gca_position(1) + (gca_position(3) - gca_x_length_real / gcf_position(3)) / 2, ...
    gca_position(2) + (gca_position(4) - gca_y_length_real / gcf_position(4)) / 2, ...
    gca_x_length_real / gcf_position(3), gca_y_length_real / gcf_position(4)];
% 图窗高
gca_y_length = gcf_position(4) * gca_position(4);
%% 获取箭头参数
error_nan = true;
for i = 1 : size(var_vec.XData, 2)
    if ~ isnan(var_vec.XData(1, i)) && ~ isnan(var_vec.YData(1, i))
        if u(i) ^ 2 + v(i) ^ 2 > 0
            error_nan = false;
            break
        end
    end
end
if error_nan
    error('你输入的向量阵全是0或缺省值')
end
x = var_vec.XData(:, i);
y = var_vec.YData(:, i);
% 第一个箭头总长度对应的像素、对应的向量大小、箭头比例尺总长度对应的像素
arrow_length = sqrt((x(4, :) - (x(1, :) + x(7, :)) / 2) .^ 2 + ...
    (y(4, :) - (y(1, :) + y(7, :)) / 2) .^ 2) * dx_real;
arrow_value = sqrt(u(i) ^ 2 + v(i) ^ 2);
arrow_length = arrow_length / arrow_value * u0;
% 箭头尾线宽、头部宽、头部长对应的像素
arrow_width = sqrt((x(7, :) - x(1, :)) .^ 2 + ...
    (y(7, :) - y(1, :)) .^ 2) * dx_real;
head_width = sqrt((x(3, :) - x(5, :)) .^ 2 + ...
    (y(3, :) - y(5, :)) .^ 2) * dx_real;
head_length = sqrt((x(4, :) - (x(3, :) + x(5, :)) / 2) .^ 2 + ...
    (y(4, :) - (y(3, :) + y(5, :)) / 2) .^ 2) * dx_real;
% 箭头边线宽度，边线颜色，填充颜色，和不按照坐标区域裁剪对象
arrow_linewidth = var_vec.LineWidth;
arrow_edgecolor = var_vec.EdgeColor;
arrow_facecolor = var_vec.FaceColor;
%% 确定矢量图的箭头比例尺的各种参数并绘制箭头比例尺
% 标注文本内容、字体、字号、字颜色、字加粗、字倾斜、框线条类型、框颜色、
% 框宽、填充颜色、填充透明度、框位置、文本位置的默认值
new_label = num2str(u0);
new_fontname = 'Helvetica';
new_fontsize = gca_y_length * 0.06 * 0.6;
new_fontcolor = [0, 0, 0];
new_fontweight = 'normal';
new_fontangle = 'normal';
new_boxstyle = '-';
new_boxcolor = [0, 0, 0];
new_boxwidth = 0.5;
new_facecolor = [1, 1, 1];
new_facealpha = 1;
new_position_unsure = 'ru';
new_labelposition_unsure = 'd';
% 箭头比例尺的原始坐标
new_XData = [0; arrow_length - head_length; arrow_length - head_length; ...
    arrow_length; arrow_length - head_length; ...
    arrow_length - head_length; 0] / gcf_position(3);
new_YData = [(head_width + arrow_width) / 2; ...
    (head_width + arrow_width) / 2; head_width; head_width / 2; ...
    0; (head_width - arrow_width) / 2; ...
    (head_width - arrow_width) / 2] / gcf_position(4);
% 用户给定的参数
if nargin > 4
    i = 1;
    while i <= length(varargin)
        if ~ isa(varargin{i}, 'char')
            error("属性名称应为char类型。如文本'lable'，字号'fontsize'等。")
        else
            switch lower(varargin{i})
                case 'label'
                    new_label = varargin{i + 1};
                case 'fontname'
                    new_fontname = varargin{i + 1};
                case 'fontsize'
                    new_fontsize = varargin{i + 1};
                case 'fontcolor'
                    new_fontcolor = varargin{i + 1};
                case 'fontweight'
                    new_fontweight = varargin{i + 1};
                case 'fontangle'
                    new_fontangle = varargin{i + 1};
                case 'position'
                    new_position_unsure = varargin{i + 1};
                case 'labelposition'
                    new_labelposition_unsure = varargin{i + 1};
                case 'arrowposition'
                    new_arrowposition_unsure = varargin{i + 1};
                    new_arrowposition_manual = true;
                case 'boxstyle'
                    new_boxstyle = varargin{i + 1};
                case 'boxcolor'
                    new_boxcolor = varargin{i + 1};
                case 'boxwidth'
                    new_boxwidth = varargin{i + 1};
                case 'facecolor'
                    new_facecolor = varargin{i + 1};
                case 'facealpha'
                    new_facealpha = varargin{i + 1};
                otherwise
                    error('未找到您所要的属性名。')
            end
            i = i + 2;
        end
    end
end
% 新图窗的长和高
new_xlim = arrow_length / gcf_position(3) * 1.2;
new_ylim = (new_fontsize / 0.6 + head_width * 1.5) / gcf_position(4);
% 新图窗位置
if isa(new_position_unsure, 'char')
    switch lower(new_position_unsure)
        case {'lu', 'left_up', 'lt', 'left_top', '左上'}
            x_begin = gca_position(1);
            y_begin = gca_position(2) + gca_position(4) - new_ylim;
            ax_position = [x_begin, y_begin, new_xlim, new_ylim];
        case{'ru', 'right_up', 'rt', 'right_top', '右上'}
            x_begin = gca_position(1) + gca_position(3) - new_xlim;
            y_begin = gca_position(2) + gca_position(4) - new_ylim;
            ax_position = [x_begin, y_begin, new_xlim, new_ylim];
        case{'ld', 'left_down', 'lb', 'left_bottom', '左下'}
            x_begin = gca_position(1);
            y_begin = gca_position(2);
            ax_position = [x_begin, y_begin, new_xlim, new_ylim];
        case{'rd', 'right_down', 'rb', 'right_bottom', '右下'}
            x_begin = gca_position(1) + gca_position(3) - new_xlim;
            y_begin = gca_position(2);
            ax_position = [x_begin, y_begin, new_xlim, new_ylim];
        otherwise
            error("图窗位置参数输入不正确，有'右上'，'右下'等。")
    end
else
    if ~ isa(new_position_unsure, 'numeric')
        error(["属性'position'后面对应的值应为4个元素的向量，[最左的", ...
            "位置，最下的位置，长，高]，是标准化的数，原则上介于0～1之间。"])
    else
        if length(new_position_unsure) ~= 4
            error(["属性'position'后面对应的值应为4个元素的向量，[最左的", ...
            "位置，最下的位置，长，高]，是标准化的数，原则上介于0～1之间。"])
        else
            x_begin = new_position_unsure(1) + ...
                new_position_unsure(3) / 2 - new_xlim / 2;
            y_begin = new_position_unsure(2) + ...
                new_position_unsure(4) / 2 - new_ylim / 2;
            ax_position = new_position_unsure;
        end
    end
end
ax = axes('position', ax_position, 'XLim', ...
    [ax_position(1), ax_position(1) + ax_position(3)], ...
    'YLim', [ax_position(2), ax_position(2) + ax_position(4)]);
an = patch(ax, [ax_position(1), ax_position(1) + ax_position(3), ...
    ax_position(1) + ax_position(3), ax_position(1)], [ax_position(2), ...
    ax_position(2), ax_position(2) + ax_position(4), ax_position(2) + ...
    ax_position(4)], new_boxcolor, 'LineStyle', new_boxstyle, ...
    'LineWidth', new_boxwidth, 'FaceColor', ...
    new_facecolor, 'FaceAlpha', new_facealpha);
% 文本位置和箭头坐标位置
if isa(new_labelposition_unsure, 'char')
    switch lower(new_labelposition_unsure)
        case {'u', 'up', 't', 'top', '上'}
            txt_position = [x_begin + new_xlim / 2, y_begin + ...
                (head_width * 1.5 + new_fontsize / 1.2) / gcf_position(4)];
            arr_position = [x_begin + new_xlim / 12, ...
                y_begin + head_width * 0.5 / gcf_position(4)];
        case{'d', 'down', 'b', 'bottom', '下'}
            txt_position = [x_begin + new_xlim / 2, ...
                y_begin + new_fontsize / 1.2 / gcf_position(4)];
            arr_position = [x_begin + new_xlim / 12, ...
                y_begin + new_fontsize / 0.6 / gcf_position(4)];
        otherwise
            error("文本位置参数输入不正确，请输入'上'或'下'。")
    end
else
    if ~ isa(new_labelposition_unsure, 'numeric')
        error("属性'labelposition'后面的值应为2个元素的向量，[x, y]。")
    else
        if length(new_labelposition_unsure) ~= 2
            error("属性'labelposition'后面的值应为2个元素的向量，[x, y]。")
        else
            txt_position = new_labelposition_unsure;
            if new_labelposition_unsure(2) > y_begin + new_ylim / 2
                arr_position = [x_begin + new_xlim / 12, ...
                    y_begin + head_width * 0.5 / gcf_position(4)];
            else
                arr_position = [x_begin + new_xlim / 12, ...
                    y_begin + new_fontsize / 0.6 / gcf_position(4)];
            end
        end
    end
end
if exist('new_arrowposition_unsure', 'var')
    if isa(new_labelposition_unsure, 'char')
        switch lower(new_labelposition_unsure)
            case {'u', 'up', 't', 'top', '上'}
                if exist('new_arrowposition_manual', 'var')
                    arr_position = [x_begin + new_xlim / 12, ...
                        y_begin + head_width * 0.5 / gcf_position(4)];
                else
                    txt_position = [x_begin + new_xlim / 2, ...
                        y_begin + new_fontsize / 1.2 / gcf_position(4)];
                    arr_position = [x_begin + new_xlim / 12, ...
                        y_begin + head_width * 0.5 / gcf_position(4)];
                end
            case{'d', 'down', 'b', 'bottom', '下'}
                if exist('new_arrowposition_manual', 'var')
                    arr_position = [x_begin + new_xlim / 12, ...
                        y_begin + new_fontsize / 0.6 / gcf_position(4)];
                else
                    txt_position = [x_begin + new_xlim / 2, y_begin + ...
                        (head_width * 1.5 + new_fontsize / 1.2) / ...
                        gcf_position(4)];
                    arr_position = [x_begin + new_xlim / 12, ...
                        y_begin + new_fontsize / 0.6 / gcf_position(4)];
                end
            otherwise
                error("箭头位置参数输入不正确，请输入'上'或'下'。")
        end
    else
        if ~ isa(new_arrowposition_unsure, 'numeric')
            error("属性'arrowposition'后面的值应为[x, y]。")
        else
            if length(new_arrowposition_unsure) ~= 2
                error("属性'arrowposition'后面的值应为[x, y]。")
            else
                arr_position = [new_arrowposition_unsure(1) - ...
                    arrow_length / gcf_position(3) / 2, ...
                    new_arrowposition_unsure(2) - head_width / ...
                    gcf_position(4) / 2];
                if ~ exist('new_arrowposition_manual', 'var')
                    if new_arrowposition_unsure(2) >= y_begin + ...
                            new_ylim / 2
                        txt_position = [x_begin + new_xlim / 2, ...
                            y_begin + new_fontsize / 1.2 / ...
                            gcf_position(4)];
                    elseif new_arrowposition_manual
                        txt_position = [x_begin + new_xlim / 2, ...
                            y_begin + (head_width * 1.5 + ...
                            new_fontsize / 1.2) / gcf_position(4)];
                    end
                end
            end
        end
    end
end
txt = text(ax, txt_position(1), txt_position(2), new_label, ...
    'HorizontalAlignment', 'center', 'FontName', new_fontname, ...
    'FontSize', new_fontsize, 'Color', new_fontcolor, 'FontWeight', ...
    new_fontweight, 'FontAngle', new_fontangle);
arr = patch(ax, new_XData + arr_position(1), new_YData + ...
    arr_position(2), arrow_facecolor, 'EdgeColor', arrow_edgecolor, ...
    'LineWidth', arrow_linewidth, ...
    'Clipping', 'off');
set(ax,'Clipping', 'off', 'XTick', [], 'YTick', [], ...
    'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
switch nargout
    case 1
        varargout{1} = arr;
    case 2
        varargout{1} = arr;
        varargout{2} = txt;
    case 3
        varargout{1} = arr;
        varargout{2} = txt;
        varargout{3} = ax;
    case 4
        varargout{1} = arr;
        varargout{2} = txt;
        varargout{3} = ax;
        varargout{4} = an;
end
end