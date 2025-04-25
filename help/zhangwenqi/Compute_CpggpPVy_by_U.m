function [ Cp,Cg,Cgp,PVy,Uyy ] =Compute_CpggpPVy_by_U( Uu, Vu)


alfa =-1; 
omiga=7.2722*10^(-5);
R=6.371*10^6;
pi=3.1416;
fai0=55*pi/180;
Ly=5;
L=10^6;Rd=10^6;Uo=10;
  F=(L/Rd)^2;
Ld=3.1*10^6;
% F=1/(Ld^2);
k0=1/(6.371*cos(fai0));
k=2*k0;
  beta=L^2/Uo*2*omiga/R*cos(fai0); 
%  beta=2*omiga/R*cos(fai0); 
m=alfa*2*pi/Ly;
dy=2*pi*2.5*6.371/360;
for ilat=1:37
    dx(ilat)=2*pi*2.5*6.371/360*cos((90-2.5*(ilat-1))*pi/180);
end
% Uu （144，37，131，39）
for iyear=1:39
    for iday=1:131
        for ilon=1:144
            for ilat=2:36
                Uyy(ilon,ilat,iday,iyear)=(Uu(ilon,ilat-1,iday,iyear)+Uu(ilon,ilat+1,iday,iyear)-2*Uu(ilon,ilat,iday,iyear))/(dy^2);
            end
            Uyy(ilon,1,iday,iyear)=Uyy(ilon,2,iday,iyear);
            Uyy(ilon,37,iday,iyear)=Uyy(ilon,36,iday,iyear);
        end
    end
end
Vxy(1:144,1:37,1:131,1:39)=0;
for iyear=1:39
    for iday=1:131
        for ilon=2:143
            for ilat=2:36
                Vxy(ilon,ilat,iday,iyear)=(Vu(ilon+1,ilat+1,iday,iyear)+Vu(ilon-1,ilat-1,iday,iyear)-Vu(ilon-1,ilat+1,iday,iyear)-Vu(ilon+1,ilat-1,iday,iyear))/(4*dy*dx(ilat));
            end
        end
        Vxy(:,1,iday,iyear)=Vxy(:,2,iday,iyear);
        Vxy(:,37,iday,iyear)=Vxy(:,36,iday,iyear);
        Vxy(1,:,iday,iyear)=Vxy(2,:,iday,iyear);
        Vxy(144,:,iday,iyear)=Vxy(143,:,iday,iyear);
    end
end


            
            
PVy=beta+F.*Uu+Vxy-Uyy;  %计算 PVy（144，37，131，39）

Cp=Uu-PVy/(k^2+m^2+F);   %  计算Cp
Cg=Uu-PVy*(m^2+F-k^2)/((k^2+m^2+F)^2);   %   计算 Cg
Cgp=Cg-Cp;   %    计算  Cgp（144，37，131，39）



end

