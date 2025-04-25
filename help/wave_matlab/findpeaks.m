
function n = findpeaks(x)
% Find peaks.找到极值 ,n为极值点所在位置
% n = findpeaks(x)
n    = find(diff(diff(x) > 0) < 0);
u    = find(x(n+1) > x(n));
n(u) = n(u)+1;