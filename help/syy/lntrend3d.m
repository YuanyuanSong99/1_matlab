function [trd h0] = lntrend3d(varin,alp)
% This function is to test linear trend of a 3D variable
% varin is input variable with lon*lat*t
% alp is the significance level of student's t test. For example, alpha
% is 0.05 meaning 95% significance level.
var = permute(varin,[3 1 2]);
d1 = size(var,1); d2 = size(var,2); d3 = size(var,3);
h0 = zeros(d2,d3);
x = [1:d1]';
for i=1:d2;
    for j=1:d3;
        m3=var(:,i,j);
        par(i,j,:)=polyfit(x,m3,1); % two regression parameters 
        trd(i,j) = par(i,j,1); % linear trend
        y1=polyval(permute(par(i,j,:),[3 2 1]),x);
        Q = sum((m3-y1).^2); % t test
        c = 1/sum((x-mean(x)).^2);
        T = par(i,j,1)/sqrt(c)/sqrt(Q/(d1-2));
        if abs(T) > tinv(1-alp/2,d1-2); % 双侧 t test
            h0(i,j) = 1; % significant grid with 1; insignificant 0
        end
    end
end
end