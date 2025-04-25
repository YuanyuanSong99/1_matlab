clc;clear all;
% ncdisp('D:\data\u_daily_500_2.5\u.daily.2.5.1979_500.nc'); %查看nc文件信息


datadir='D:\data\u_daily_500_2.5\';  %指定文件夹
filelist=dir([datadir,'*.nc']);  %指定文件类型
k=length(filelist);  % k 为文件数
 U(1:144,1:37,1:1464,1:k)=0;
for year=1:k
    filename=[datadir,filelist(year).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    uData    = ncread(filename,'u'); 
    LonData  = ncread(filename,'longitude');
    LatData  = ncread(filename,'latitude'); 
%     TimeData = ncread(filename,'time');
    netcdf.close(ncid); 
    nvar=size(uData); 
    U(1:144,1:37,1:nvar(3),year)=uData;
end

%% 

%（一）处理数据   Uu （144，37，131，39） 第一个数据为90N，第37个为赤道
dim=10;
[ Uu] = Deal( U,dim );  %Uu 为you量纲的西风U/U0（U0=1）原始场

[ Vu] = Deal( V,dim );  %Uu 为you量纲的西风U/U0（U0=1）原始场

%（二）计算Cg-Cp-Cgp-PVy by U 
[ Cp,Cg,Cgp,PVy,Uyy ] =Compute_CpggpPVy_by_U( Uu,Vu);

% (三) 区域平均的Cp Cg Cgp
nx_a=85;nx_b=105;
ny_a=5;ny_b=21;
[ Cp_ave_area,Cg_ave_area,Cgp_ave_area ] = C_ave_area( Cp,Cg,Cgp,nx_a,nx_b,ny_a,ny_b );

% (四) 计算Cg，Cp，Cgp 序列 lag-10到lag10
% load UB_lag0_2.mat
i=1;
for iyear=1:39
    for iday=1:131
        if AB_lag0(iday,iyear)==999    %%  U B
            for lag=1:41
            Cp_use(lag,i)=Cp_ave_area(iday-21+lag,iyear);
            Cg_use(lag,i)=Cg_ave_area(iday-21+lag,iyear);
            Cgp_use(lag,i)=Cgp_ave_area(iday-21+lag,iyear);
            PVy_use(:,:,lag,i)=PVy(:,:,iday+lag-21,iyear);
            Uu_use(:,:,lag,i)=Uu(:,:,iday+lag-21,iyear);
            Uyy_use(:,:,lag,i)=Uyy(:,:,iday+lag-21,iyear);
            end
            i=i+1;
        end
       
    end
end
Cp_pl=mean(Cp_use,2);
Cg_pl=mean(Cg_use,2);
Cgp_pl=mean(Cgp_use,2);
PVy_lag=mean(PVy_use,4);
Uu_lag=mean(Uu_use,4);
Uyy_lag=mean(Uyy_use,4);

for j=1:64
      PVy_m(j)=mean(mean(mean(PVy_use(nx_a:nx_b,ny_a:ny_b,1:6,j))));
end
PVy_mm=mean(PVy_m);
Std_PVy=std(PVy_mm);
i=1;
for in=1:nB
   if( (PVy_m(in)<0.5*Std_PVy+PVy_mm))% &&(U_x_mean(in)<mean_x_U-0.1*std_x_U) ) 
     
    
       Z_use(:,:,:,i)=Z_ano_lag(:,:,:,in);  %    Z_ano_use 177 37 21  上游强西风的个数   距平场
       PVy_use_2(:,:,:,i)=PVy_use(:,:,:,in);
       i=i+1;
   end 
end
% （五） 画图
% figure(1); 
% t=1:21;t=t';
%   c = polyfit(t, Cgp_pl, 6); %进行拟合，c为2次拟合后的系数
%   d = polyval(c, t); %拟合后，每一个横坐标对应的值即为d
%   plot(t, d, 'k-','LineWidth',2);hold on; %拟合后的曲线
% 
%   c = polyfit(t, Cp_pl, 6); %进行拟合，c为2次拟合后的系数
%   d = polyval(c, t); %拟合后，每一个横坐标对应的值即为d
%   plot(t, d, 'b-','LineWidth',2);hold on; %拟合后的曲线
%   
%   c = polyfit(t, Cg_pl, 6); %进行拟合，c为2次拟合后的系数
%   d = polyval(c, t); %拟合后，每一个横坐标对应的值即为d
%   plot(t, d, 'r-','LineWidth',2);hold on; %拟合后的曲线
% 
% plot(t,Cgp_pl,'k');hold on;
% plot(t,Cp_pl,'b-');hold on;
% plot(t,Cg_pl,'r-');hold on;
% 
% set(gca,'XLim',[1 21]);
% set(gca,'xtick',1:1:21);
% set(gca, 'xtickLabel' ,{'-10','-9','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4',...
%     '5','6','7','8','9','10'}) ;
% xlabel('lag','Fontsize',12)
% grid on;
% 
% 
% for ilon=1:144
%     for ilat=1:37
%          PVy_mean(ilon,ilat)=mean(mean(PVy(ilon,ilat,21:110,:)));             
%     end
% end
figure(2);
[lat, lon]=meshgrid(LatData,LonData); 
% %   m_proj('Equidistant Cylindrical','lat',[0 90],'lon',[-180 180]);
% %    m_coast('patch',[1 1 1],'edgecolor','none');
%   m_proj('ortho','lat',90,'lon',60); m_grid('linestyle',':','xtick',[-180:60:180],'ytick',[0:30:90]);hold on;
% 
% % pl=reshape(PVy_mean(:,:,11),[144 37]);
% [ch,s]=m_contourf(lon,lat,PVy_mean,[-10:1:10]);hold on;
% %   clabel(ch,s);
% colormap(jet);

m_proj('stereographic','lat',90,'long',0,'radius',80);
m_coast('linewidth',1,'color','k');
m_grid('linestyle','--','xtick',[-90 0 90 180],'ytick',4,...
    'tickdir','in','fontsize',10);
hold on

%[-20*10^-10:0.1*10^-10:20*10^-10]

[C,h]=m_contourf(lon,lat,PVy_mean,[-20:0.1:20],'linestyle','none');
colorbar;
colormap(jet);




