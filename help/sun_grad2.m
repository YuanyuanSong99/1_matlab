function  [ux,uy]=sun_grad2(u,lon,lat)
%-----------------------------------------------------
% [u,v]=sun_grad2(u,lon,lat)
% horizontal gradient [m,n] --> [m,n]
% assuming [1:m] is northward and positive
% Central differential
%-----------------------------------------------------

[m,n,k]=size(u); dx=ones(m,n)*NaN; dy=dx; ux=u*NaN; uy=ux;
ux(:,2:n-1,:)=u(:,3:n,:)-u(:,1:n-2,:); 
uy(2:m-1,:,:)=u(3:m,:,:)-u(1:m-2,:,:);
dlon(2:n-1)=lon(3:n)-lon(1:n-2);
dlat(2:m-1)=lat(3:m)-lat(1:m-2);

for j=2:m-1,
 for i=2:n-1, 
  dx(j,i)=dlon(i)/abs(dlon(i))*sw_dist([lat(j) lat(j)],[0,dlon(i)],'km')*1e3;
 end
dy(j,:)=dlat(j)/abs(dlat(j))*sw_dist([0, dlat(j)],[0,0],'km')*1e3;
end

for i=1:k,
ux(:,:,i)=ux(:,:,i)./dx; uy(:,:,i)=uy(:,:,i)./dy;
end
ux(:,1,:)=ux(:,2,:);ux(:,n,:)=ux(:,n-1,:);ux(1,:,:)=ux(2,:,:);ux(m,:,:)=ux(m-1,:,:);
uy(:,1,:)=uy(:,2,:);uy(:,n,:)=uy(:,n-1,:);uy(1,:,:)=uy(2,:,:);uy(m,:,:)=uy(m-1,:,:);
