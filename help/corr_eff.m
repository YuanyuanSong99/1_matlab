%   [r,p,n_eff] = corr_eff(x1,x2,a)
%   This function is used to calculate correlation coefficient between x1 
% and x2, under ttest by effective freedom degree.
%   这个函数用于计算两序列的相关系数，并使用有效自由度进行t检验。
%   x1 and x2 are 1D data.
%   a is alpha for ttest.
%   r is the correlation coefficient between x1 and x2.
%   p is the result of ttest. If p == 1, they pass the test; if p == 0, they 
% fail to pass the test.
%   n_eff is effective freedom degree. If n_eff <= 2, we think the
% correlation cofficient is invalid.

% Author: Syy, IAP, CHINA.
function [r,p,n_eff2] = corr_eff(x1,x2,a);
si1 = size(x1,1);
if si1 ~= 1
    x1 = x1';
else
    si1 = size(x1,2);
end
si2 = size(x2,1);
if si2 ~= 1
    x2 = x2';
else
    si2 = size(x2,2);
end
p = 0;
r0 = corrcoef(x1,x2);
r = r0(1,2);
zi1 = corrcoef(x1(1:end-1),x1(2:end));
zi2 = corrcoef(x2(1:end-1),x2(2:end));
n_eff1 = si1*(1-zi1(1,2)*zi2(1,2))/(1+zi1(1,2)*zi2(1,2)); % statistical effective freedom degree (too strict)
for t = 1:si1-2;
    temp1 = corrcoef(x1(1:end-t),x1(t+1:end));
    x1x1(t) = temp1(1,2); 
    temp2 = corrcoef(x2(1:end-t),x2(t+1:end));
    x2x2(t) = temp2(1,2); 
    C(t) = (si1-t)/si1;
end
n_eff2 = si1/(1+2*sum(C.*x1x1.*x2x2)); % effective freedom degree from Li Jianping et al.
if n_eff2 > 2
    t_val = sqrt(n_eff2-2)*r/sqrt(1-r^2);
    if abs(t_val) > tinv(1-a,n_eff2-2); % 单侧 t test
        p = 1;
    end
end
    
