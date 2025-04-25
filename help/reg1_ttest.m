% 一元一次回归方程 & 系数的t检验
% x and y should be 1D data.
% alp is alpha(test level) for t-test
% if opt = 0, freedom is n0-2; if opt = 1, freedom is the effective freedom
% for decadal variation calculated from Li Jianping et al. 
% par are coefficients.
% t_val is t statistic.
% If the coefficient passes the test, h = 1; else, h = 0.
function [par,h,t_val] = reg1_ttest(x,y,alp,opt)
    n0 = length(x);
    x = reshape(x,1,n0);
    y = reshape(y,1,n0);
    par1 = polyfit(x,y,1);
    par = par1(1);
    yp = polyval(par,x);
    h = 0;
    if opt == 0
        t_val = par(1)/sqrt((sum((y-yp).^2)/(n0-2))/sum((x-mean(x)).^2));
        if abs(t_val) > tinv(1-alp/2,n0-2); % 双侧 t test
            h = 1;
        end
    else
        for t = 1:n0-2;
            temp1 = corrcoef(x(1:end-t),x(t+1:end));
            xx(t) = temp1(1,2);
            temp2 = corrcoef(y(1:end-t),y(t+1:end));
            yy(t) = temp2(1,2);
            C(t) = (n0-t)/n0;
        end
        n_eff2 = n0/(1+2*sum(C.*xx.*yy)); % effective freedom degree from Li Jianping et al.
        t_val = par(1)/sqrt((sum((y-yp).^2)/(n_eff2))/sum((x-mean(x)).^2));
        if abs(t_val) > tinv(1-alp,n_eff2); % 单侧 t test
            h = 1;
        end
    end
end