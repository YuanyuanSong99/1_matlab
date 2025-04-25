function [ Uu] = Deal( U,dim)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

Zv=reshape(U,[144,37,4,366,40]);
Zm=mean(Zv,3);
Zuse=reshape(Zm,[144 37 366 40]);
Zp=Zuse/dim;  % 无量纲化
%数据是从北向南：第一个数据是北极，最后一个数据是赤道

Uu(1:144,1:37,1:131,1:39)=0; 
for iyear=1:39
    if Zp(1,1,366,iyear)==0
        Uu(:,:,1:31+20,iyear)=Zp(:,:,365-50:365,iyear);
        if Zp(1,1,366,iyear+1)==0
            Uu(:,:,52:52+31+28+20-1,iyear)=Zp(:,:,1:31+28+20,iyear+1);
        elseif Zp(1,1,366,iyear+1)~=0
            Uu(:,:,52:52+31+29+20-1,iyear)=Zp(:,:,1:31+29+20,iyear+1); 
        end
    elseif Zp(1,1,366,iyear)~=0
        Uu(:,:,1:31+20,iyear)=Zp(:,:,366-50:366,iyear);
        Uu(:,:,52:52+31+28+20-1,iyear)=Zp(:,:,1:31+28+20,iyear+1);
    end
end
%  Uu为每年冬天的【原始场】的数据：Uu(1)为1979.12-1980.2的数据，Uu(39)为2017.12-2018.2的数据，前后各加上20天
%    Uu(:,:,131,2)=0; 第二年的131天是零


return

end

