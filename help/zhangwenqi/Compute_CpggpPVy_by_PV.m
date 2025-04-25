function [ Cp,Cg,Cgp,PVy ] =Compute_CpggpPVy_by_PV( Uu,PVu)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
alfa =-1; 
omiga=7.2722*10^(-5);
R=6.371*10^6;
pi=3.1416;
fai0=55*pi/180;
Ly=5;
L=10^6;Rd=10^6;Uo=10;
F=(L/Rd)^2;
k0=1/(6.371*cos(fai0));
k=2*k0;
beta=L^2/Uo*2*omiga/R*cos(fai0); 
m=alfa*2*pi/Ly;
dy=2*pi*2.5*6.371/360;

for iyear=1:39
    for iday=1:131
        for ilon=1:144
            for ilat=2:36
                PVy(ilon,ilat,iday,iyear)=-(PVu(ilon,ilat+1,iday,iyear)-PVu(ilon,ilat-1,iday,iyear))/(2*dy);
            end
             PVy(ilon,37,iday,iyear)= PVy(ilon,36,iday,iyear);
             PVy(ilon,1,iday,iyear)= PVy(ilon,2,iday,iyear);
        end
    end
end

Cp=Uu-PVy/(k^2+m^2+F);   %  计算Cp
Cg=Uu-PVy*(m^2+F-k^2)/((k^2+m^2+F)^2);   %   计算 Cg
Cgp=Cg-Cp;   %    计算  Cgp（144，37，131，39）

end

