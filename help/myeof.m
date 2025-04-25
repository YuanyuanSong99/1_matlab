function [eof_maps,pc,expvar]=myeof(X,kmod)
%   EOF analysis 
%   EOF gives eigenmode maps of variability and corresponding principal component
% time series for spatiotemporal data analysis.  It is designed specifically for 3D matricies
% of data such as surface atmosphere temperatures where dimensions 1 and 2 are spatial dimensions 
% (e.g., lat and lon; lon and lat; x and y, etc.), and the third dimension represents different 
% slices or snapshots of data in time. X should be anomaly or normalized data.
X=X-mean(mean(mean(X,3),2),1);
[onum anum tnum]=size(X);
X1=reshape(X,[onum*anum,tnum]);
[m nt]=size(X1);
[V1,E]=eig(X1'*X1); 
E=fliplr(flipud(E));
V2=fliplr(X1*V1);
lamd=diag(E)/nt;  % 取主对角线上的元素
for i=1:kmod; 
      V(:,i)=V2(:,i) / sqrt(E(i,i));
      pc(i,:)=(V(:,i)'*X1)';
      expvar(i)=E(i,i)/nt/ sum(lamd);
end
eof_maps=reshape(V,[onum,anum,kmod]);

end


