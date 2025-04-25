%% À„∑Á”¶¡¶
[a,b,c]=size(U_1);
U=sqrt(U_1.^2+V_1.^2);
utao=zeros(size(U_1));vtao=zeros(size(V_1));
for i=1:a
    for j=1:b
        for t=1:c
            if  4< U(i,j,t) < 11
                Cd=1.2*10^-3;
                utao(i,j,t)=1.2*Cd*U(i,j,t)*U_1(i,j,t);
                vtao(i,j,t)=1.2*Cd*U(i,j,t)*V_1(i,j,t);
            elseif U(i,j,t) >= 11 && U(i,j,t) < 25
                Cd=(0.49+0.065*U(i,j,t))*10^-3;
                utao(i,j,t)=1.2*Cd*U(i,j,t)*U_1(i,j,t);
                vtao(i,j,t)=1.2*Cd*U(i,j,t)*V_1(i,j,t);
            elseif U(i,j,t)>=25
                Cd=2.115*10^-3;
                utao(i,j,t)=1.2*Cd*U(i,j,t)*U_1(i,j,t);
                vtao(i,j,t)=1.2*Cd*U(i,j,t)*V_1(i,j,t);
            end
            clear Cd
        end
    end
end

  for i=1:62
   curlz(:,:,i) = zh_curl(squeeze(utao(:,:,i)'),squeeze(vtao(:,:,i)'),lon,lat);
  end
curlz1=curlz.*10000000;%wind stress curl 10^-7

%----------
function vor=zh_curl(u1,v1,lon,lat)
[u1x,u1y]=sun_grad2(u1,lon,lat);
[v1x,v1y]=sun_grad2(v1,lon,lat);
vor=v1x-u1y;

