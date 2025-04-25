clc;clear all;
% ncdisp('D:\data\pv_daily_330_2.5\pv_2.5.daily.1982_pt.nc'); %查看nc文件信息


datadir='D:\data\pv_daily_330_2.5\';  %指定文件夹
filelist=dir([datadir,'*.nc']);  %指定文件类型
k=length(filelist);  % k 为文件数
 PV(1:144,1:37,1:1464,1:k)=0;
for year=1:k
    filename=[datadir,filelist(year).name];
    ncid=netcdf.open(filename,'NC_NOWRITE');
    pvData    = ncread(filename,'pv'); 
    LonData  = ncread(filename,'longitude');
    LatData  = ncread(filename,'latitude'); 
%     TimeData = ncread(filename,'time');
    netcdf.close(ncid); 
    nvar=size(pvData); 
    PV(1:144,1:37,1:nvar(3),year)=pvData;
end


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

%（一）处理数据   PV （144，37，131，39） 第一个数据为90N，第37个为赤道
dimPV=10^(-6);
[ PVu] = Deal( PV,dimPV );  %PVu 为you量纲的西风U/U0（U0=1）原始场
dimU=10;
[Uu] = Deal(U,dimU);
%（二）计算Cg-Cp-Cgp-PVy by PV
[ Cp,Cg,Cgp,PVy ] =Compute_CpggpPVy_by_PV( Uu,PVu);


% (三) 区域平均的Cp Cg Cgp
nx_a=85;nx_b=105;
ny_a=5;ny_b=21;
[ Cp_ave_area,Cg_ave_area,Cgp_ave_area ] = C_ave_area( Cp,Cg,Cgp,nx_a,nx_b,ny_a,ny_b );

% (四) 计算Cg，Cp，Cgp 序列 lag-10到lag10
load UB_lag0_2.mat
i=1;
for iyear=1:39
    for iday=1:131
        if UB_lag0(iday,iyear)==999    %%  U B
            for lag=1:21
            Cp_use(lag,i)=Cp_ave_area(iday-11+lag,iyear);
            Cg_use(lag,i)=Cg_ave_area(iday-11+lag,iyear);
            Cgp_use(lag,i)=Cgp_ave_area(iday-11+lag,iyear);
            PVy_use(:,:,lag,i)=PVy(:,:,iday+lag-11,iyear);
            Uu_use(:,:,lag,i)=Uu(:,:,iday+lag-11,iyear);
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


% (五) 画图
t=1:21;t=t';
  c = polyfit(t, Cgp_pl, 6); %进行拟合，c为2次拟合后的系数
  d = polyval(c, t); %拟合后，每一个横坐标对应的值即为d
  plot(t, d, 'k-','LineWidth',2);hold on; %拟合后的曲线

  c = polyfit(t, Cp_pl, 6); %进行拟合，c为2次拟合后的系数
  d = polyval(c, t); %拟合后，每一个横坐标对应的值即为d
  plot(t, d, 'b-','LineWidth',2);hold on; %拟合后的曲线
  
  c = polyfit(t, Cg_pl, 6); %进行拟合，c为2次拟合后的系数
  d = polyval(c, t); %拟合后，每一个横坐标对应的值即为d
  plot(t, d, 'r-','LineWidth',2);hold on; %拟合后的曲线

plot(t,Cgp_pl,'k');hold on;
plot(t,Cp_pl,'b-');hold on;
plot(t,Cg_pl,'r-');hold on;

set(gca,'XLim',[1 21]);
set(gca,'xtick',1:1:21);
set(gca, 'xtickLabel' ,{'-10','-9','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4',...
    '5','6','7','8','9','10'}) ;
xlabel('lag','Fontsize',12)
grid on;







