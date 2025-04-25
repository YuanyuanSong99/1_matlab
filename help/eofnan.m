function [eof_maps,pc,expvar,eof_maps_detrend,pc_detrend,expvar_detrend] = eofnan(X)
% X is anomaly
[nx,ny,nt] = size(X);
qsh_num1 = ismissing(X); %
qsh_num3 = sum(qsh_num1,3); % 缺省map
qsh_num4 = find(qsh_num3); % 缺省空间点标号(annual)
qsh_map = zeros(nx*ny,nt);
qsh_map(qsh_num4,:) = nan;
fqsh_num = find(qsh_map == 0); % 正常点标号
m5 = reshape(X,nx*ny,nt); % 原序列
m6 = m5;
m6(qsh_num4,:) = []; % after remove nan
fqsh_gsh = size(m6,1); % 正常空间点个数
x = [1:1:nt];
for i = 1:fqsh_gsh
par = polyfit(x,m6(i,:),1);
y = polyval(par,x);
m7(i,:) = m6(i,:) - y; % detrended
end
% detrend
kmod = nt;
[meo neo]=size(m7); % eof
[V1,E]=eig(m7'*m7);
E=fliplr(flipud(E));
V2=fliplr(m7*V1);
lamd=diag(E)/neo;  % 取主对角线上的元素
for i=1:kmod;
V(:,i)=V2(:,i) / sqrt(E(i,i));
pc_detrend(i,:)=(V(:,i)'*m7)';
expvar_detrend(i)=E(i,i)/neo/ sum(lamd);
end
meof1 = NaN(nx*ny*nt,1);
meof1(fqsh_num,1) = reshape(V,numel(V),1);
eof_maps_detrend = reshape(meof1,nx,ny,nt);
clear kmod meo neo V1 E V2 lamd V meof1 
% non-detrend
kmod = nt;
[meo neo]=size(m6); % eof
[V1,E]=eig(m6'*m6);
E=fliplr(flipud(E));
V2=fliplr(m6*V1);
lamd=diag(E)/neo;  % 取主对角线上的元素
for i=1:kmod;
V(:,i)=V2(:,i) / sqrt(E(i,i));
pc(i,:)=(V(:,i)'*m6)';
expvar(i)=E(i,i)/neo/ sum(lamd);
end
meof1 = NaN(nx*ny*nt,1);
meof1(fqsh_num,1) = reshape(V,numel(V),1);
eof_maps = reshape(meof1,nx,ny,nt);
end