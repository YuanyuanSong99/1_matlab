function [ Block_day_lon ] = Select( Zz,Nlon,nx_a,nx_b,Nday )
%Select
%   此处显示详细说明
%%LatData(5)=80N LatData(13)=60N LatData(21)=40N
n=size(Zz);
latN=5;lat0=13;latS=21;dlat=2;
%  use Zz(144,37,151,67)
Block(1:144,1:151,1:67)=0;
for iyear=1:n(4)
   for iday=1:151 
       for ilon=1:144
           if (((Zz(ilon,latN,iday,iyear)-Zz(ilon,lat0,iday,iyear))<-200)&&(Zz(ilon,lat0,iday,iyear)-Zz(ilon,latS,iday,iyear)>0))...
                        ||(((Zz(ilon,latN+dlat,iday,iyear)-Zz(ilon,lat0+dlat,iday,iyear))<-200)&&(Zz(ilon,lat0+dlat,iday,iyear)-Zz(ilon,latS+dlat,iday,iyear)>0))...
                        ||(((Zz(ilon,latN-dlat,iday,iyear)-Zz(ilon,lat0-dlat,iday,iyear))<-200)&&(Zz(ilon,lat0-dlat,iday,iyear)-Zz(ilon,latS-dlat,iday,iyear)>0))
               Block(ilon,iday,iyear)=1;
           end
        end
   end
end



Block_Nlon=Block;
for iyear=1:n(4)
    for iday=1:151
        num_lon=0;
        for ilon=nx_a:nx_b     %  2   1:144
                 if Block(ilon,iday,iyear)==1
                      num_lon=num_lon+1; 
                 else
                      if num_lon>=Nlon
                           Block_Nlon(ilon-num_lon:ilon-1,iday,iyear)=1;
                      else
                           Block_Nlon(ilon-num_lon:ilon,iday,iyear)=0;
                      end
                      num_lon=0;
                 end
        end
    end
end


B_Nlon=Block_Nlon; 
B_Nlon(1:nx_a-1,:,:)=0;
B_Nlon(nx_b+1:144,:,:)=0;
Block_day_lon=B_Nlon;

for iyear=1:n(4)
    num_day=0;
   for iday=1:151 
       if max(B_Nlon(nx_a:nx_b,iday,iyear))==1
           num_day=num_day+1;
       else
           if num_day<Nday
               Block_day_lon(:,iday-num_day:iday,iyear)=0;
           end
           num_day=0;
       end
   end
end
% B_day_lon 为 连续N天发生阻塞的天数

return

end

