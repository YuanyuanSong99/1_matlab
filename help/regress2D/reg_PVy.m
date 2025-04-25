%   This function is for regression onto DJF PVy(U500),NCEP.
%   x is 1D series.
%   lagyr means you can change the setup year for regression. Remember to
% maitain the same length with x.
%   alp is test level.
%   par is coefficients.
%   h is test results.
function [par,h] = reg_PVy(x,lagyr,alp)
addpath E:\1_matlab\help
load('E:\1_matlab\NCEP\PVy_U500_djf.mat'); % PVy
load('E:\1_matlab\NCEP\lon.mat');
load('E:\1_matlab\NCEP\lat.mat');
m0=PVy;
yrsnum=size(PVy,2);
% calculate anomaly
l=0;
k=0;
for s=1:yrsnum
    fe1=size(m0{1,s},3);
    if fe1 == 91 % notice feb 28or29
        l=l+1;
        m1(:,:,l)=m0{1,s}(:,:,91); % 29th
        m2(:,:,1:90,l)=m0{1,s}(:,:,1:90); % before 29th
    else
        k=k+1;
        m3(:,:,:,k)=m0{1,s};
    end
end
hm1=mean(m1,3); %29th
hm2=cat(4,m2,m3(:,:,1:90,:));
hm3=mean(hm2,4);
for s=1:yrsnum
    fe1=size(m0{1,s},3);
    if fe1 == 91
        m4{1,s}(:,:,1:90)=m0{1,s}(:,:,1:90)-hm3; % anomaly, nonseasonalized
        m4{1,s}(:,:,91)=m0{1,s}(:,:,91)-hm1;
    else
        m4{1,s}(:,:,1:90)=m0{1,s}(:,:,1:90)-hm3;
    end
    m5(:,:,s)=mean(m4{1,s},3);
end
for i = 1:144
    for j = 1:37
        m6(i,j,:) = detrend(permute(m5(i,j,:),[3 2 1])); % detrend
    end
end
% regress
x1(:,1)=zscore(x); % 标准化
for i=1:144
    for j=1:37
        yvar=permute(m6(i,j,lagyr:end),[3 2 1]);    
        [par(i,j,:),h(i,j),t_val] = reg1_ttest(x1,yvar,alp);  
    end
end
m_proj('lambert conformal conic','lon',[40 100],'lat',[40 70]);
m_coast('linewidth',2,'color','k');%画粗海岸线
m_grid('xtick',4,'ytick',4,'linestyle','none');
maxval = max(par(17:41,9:21,1),[],'all');
minval = min(par(17:41,9:21,1),[],'all');
if abs(maxval) > abs(minval);
    bj = ceil(abs(maxval));
else
    bj = ceil(abs(minval));
end
[C,hh]=m_contourf(lonData(17:41),latData(9:21),par(17:41,9:21,1)',[-bj:bj/20:bj],'linestyle','none');
caxis([-bj,bj]);
hold on
[dots(:,1) dots(:,2)] = find(h == 1);
m_plot(lonData(dots(:,1)),latData(dots(:,2)),'.','markersize',4,'color','k'); 
m_coast('linewidth',2,'color','k');%画粗海岸线
m_grid('xtick',4,'ytick',4,'linestyle','none');
load('E:/1_matlab/help/colorbar_mat/bl_re.mat');
xg=cat(1,bl_re(1:112,:),bl_re(162:end,:));
colormap(double(xg)/255);
% set colorbar
hb = colorbar('location','southoutside');
set(hb,'fontsize',16,'Units','normalized','position',[0.2 0.08 0.6 0.02]);
