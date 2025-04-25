 %   [llcoef,p,n_eff] = leadlagcor(x1,x2,lyrs,a)
%   x1 and x2 should be 1D data; 
%   lyrs should be 11,21,31,41...;
%   a is alpha for ttest;
%   llcoef will be "x1 lead x2"
%   p is the result of ttest. If p == 1, they pass the test; if p == 0, they 
% fail to pass the test.
%   n_eff is effective freedom degree. If n_eff <= 2, we think the
% correlation cofficient is invalid.

% Author: Syy, IAP, CHINA.
function [llcoef,p,n_eff] = leadlagcor(x1,x2,lyrs,a);
addpath /Volumes/Togo4T/1_matlab/help     %  call function corr_eff.m
si1 = size(x1,1);
if si1 ~= 1
    x1 = x1';
end
si2 = size(x2,1);
if si2 ~= 1
    x2 = x2';
end
leng=size(x1,2);
for yr=1:lyrs
    [c2d(yr),p2d(yr),n_eff2(yr)] = corr_eff(x1(yr:leng),x2(1:leng-yr+1),a); % x2 lead
    [c4d(yr),p4d(yr),n_eff4(yr)] = corr_eff(x1(1:leng-lyrs+yr),x2(lyrs+1-yr:leng),a); % x2 lag
end
c2d(1)=[]; 
p2d(1)=[];
n_eff2(1)=[];
llcoef=[c4d,c2d];
llcoef=double(llcoef)';
p=[p4d,p2d];
p=double(p');
n_eff=[n_eff4,n_eff2];
n_eff=double(n_eff');
end
%  plot
% plot(llcoef,'k','linewidth',2);
% hold on
% c11d=find(p90==1);
% scatter(c11d,llcoef(c11d),'o','g','linewidth',2);
% c9d=find(p95==1);
% scatter(c9d,llcoef(c9d),'*','b','linewidth',2);
% ht1d=text(14,0.9,'>95%  *','fontsize',12);
% set(ht1d,'color',[0 0 1])
% ht2d=text(22,0.9,'>90%  o','fontsize',12);
% set(ht2d,'color',[0 1 0])
% c55d=abs(llcoef(1:lyrs));
% c6d=find(c55d==max(c55d));
% cn=size(c6d,2);
% if cn > 1
%     c66d=c6d(1);
% else
%     c66d=c6d;
% end
% c7d=roundn(llcoef(c66d),-2);
% scatter(c66d,c7d,'+','r','linewidth',2); % maximum abs
% if c7d > 0
%     text(c66d,c7d+0.1,num2str(c7d));
% else
%     text(c66d,c7d-0.1,num2str(c7d));
% end
% set(gca,'XLim',[1,lyrs*2-1]);
% set(gca,'XTick',[1:5:lyrs*2-1]);
% set(gca,'XTickLabel',[-(lyrs-1):5:lyrs-1],'FontSize',12);
% set(gca,'YLim',[-1,1]);
% set(gca,'YTick',[-1:0.5:1]);
% xlabel(['x1 leads              Lag year             ','x1 lags'],'FontSize',12);
% ylabel('R','FontSize',12);
