function [U_svd,A_svd,V_svd,B_svd,scf_svd] = SVD(t,h,kmod)
% 
% ---svd
% 
P_net = t;
Q_net = h;
clear t h;
[nx1,ny1,nt] = size(P_net);
[nx2,ny2,nt] = size(Q_net);
np = nx1*ny1;
nq = nx2*ny2;
p = reshape(P_net,np,nt);
q = reshape(Q_net,nq,nt);
% 
% ---calculate the month anomaly of var
nn = mod(nt,12);
for m = 1:12;
    mon_mean1(:,m) = mean(p(:,m:12:nt),2);
    mon_mean2(:,m) = mean(q(:,m:12:nt),2);
    if(m<= nn)
      p(:,m:12:nt) = p(:,m:12:nt)-mon_mean1(:,m)*ones(1,floor(nt/12)+1);
      q(:,m:12:nt) = q(:,m:12:nt)-mon_mean2(:,m)*ones(1,floor(nt/12)+1);
    else
      p(:,m:12:nt) = p(:,m:12:nt)-mon_mean1(:,m)*ones(1,floor(nt/12)); 
      q(:,m:12:nt) = q(:,m:12:nt)-mon_mean2(:,m)*ones(1,floor(nt/12)); 
    end
end
% 
% remove the point has NaN value from the data, and land is an index which is 1 if there is NaN means land
% ---p_svd
p = p';
p_svd = p;
land1(1:np) = 0;
for i = 1:np
    index = find(isnan(p(:,i)));
    if(length(index)~=0);
        land1(i)=1;       % this means here is land
    end;
end;
if(isempty(find(land1==1))==1)
else
  p_svd(:,find(land1==1))=[];  % delete the point where there is NaN
end;
% --q_svd
q = q';
q_svd = q;
land2(1:nq) = 0;
for i = 1:nq
    index = find(isnan(q(:,i)));
    if(length(index)~=0);
        land2(i)=1;       
    end;
end;
if(isempty(find(land2==1))==1)
else
  q_svd(:,find(land2==1))=[];  
end;
% 
% ---the data to svd
C = p_svd'*q_svd;
[U,L,V] = svd(C);
A = p_svd*U;
B = q_svd*V;
l = diag(L);
scf = l.^2 / sum(l.^2);
A_svd = A(:,1:kmod); 
B_svd = B(:,1:kmod); 
scf_svd = scf(1:kmod);
%
% ---put U back as land  to Ul (Uland) ,U_svd(nx,ny,kmod)
Ul(1:np,1:kmod) = NaN; 
k = 1;
for i = 1:np;  
    if(land1(i)==0) 
        Ul(i,:) = U(k,1:kmod); 
        k = k+1; 
    end;
end;
U_svd = reshape(Ul,nx1,ny1,kmod);
% ---put V back as land  to Vl (Vland) ,V_svd(nx2,ny2,kmod)
Vl(1:nq,1:kmod) = NaN; 
k = 1;
for i = 1:nq;  
    if(land2(i)==0) 
        Vl(i,:) = V(k,1:kmod); 
        k = k+1; 
    end;
end;
V_svd = reshape(Vl,nx2,ny2,kmod);
%