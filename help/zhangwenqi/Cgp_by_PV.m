clc;clear all;
% ncdisp('D:\data\pv_daily_330_2.5\pv_2.5.daily.1982_pt.nc'); %�鿴nc�ļ���Ϣ


datadir='D:\data\pv_daily_330_2.5\';  %ָ���ļ���
filelist=dir([datadir,'*.nc']);  %ָ���ļ�����
k=length(filelist);  % k Ϊ�ļ���
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


datadir='D:\data\u_daily_500_2.5\';  %ָ���ļ���
filelist=dir([datadir,'*.nc']);  %ָ���ļ�����
k=length(filelist);  % k Ϊ�ļ���
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

%��һ����������   PV ��144��37��131��39�� ��һ������Ϊ90N����37��Ϊ���
dimPV=10^(-6);
[ PVu] = Deal( PV,dimPV );  %PVu Ϊyou���ٵ�����U/U0��U0=1��ԭʼ��
dimU=10;
[Uu] = Deal(U,dimU);
%����������Cg-Cp-Cgp-PVy by PV
[ Cp,Cg,Cgp,PVy ] =Compute_CpggpPVy_by_PV( Uu,PVu);


% (��) ����ƽ����Cp Cg Cgp
nx_a=85;nx_b=105;
ny_a=5;ny_b=21;
[ Cp_ave_area,Cg_ave_area,Cgp_ave_area ] = C_ave_area( Cp,Cg,Cgp,nx_a,nx_b,ny_a,ny_b );

% (��) ����Cg��Cp��Cgp ���� lag-10��lag10
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


% (��) ��ͼ
t=1:21;t=t';
  c = polyfit(t, Cgp_pl, 6); %������ϣ�cΪ2����Ϻ��ϵ��
  d = polyval(c, t); %��Ϻ�ÿһ���������Ӧ��ֵ��Ϊd
  plot(t, d, 'k-','LineWidth',2);hold on; %��Ϻ������

  c = polyfit(t, Cp_pl, 6); %������ϣ�cΪ2����Ϻ��ϵ��
  d = polyval(c, t); %��Ϻ�ÿһ���������Ӧ��ֵ��Ϊd
  plot(t, d, 'b-','LineWidth',2);hold on; %��Ϻ������
  
  c = polyfit(t, Cg_pl, 6); %������ϣ�cΪ2����Ϻ��ϵ��
  d = polyval(c, t); %��Ϻ�ÿһ���������Ӧ��ֵ��Ϊd
  plot(t, d, 'r-','LineWidth',2);hold on; %��Ϻ������

plot(t,Cgp_pl,'k');hold on;
plot(t,Cp_pl,'b-');hold on;
plot(t,Cg_pl,'r-');hold on;

set(gca,'XLim',[1 21]);
set(gca,'xtick',1:1:21);
set(gca, 'xtickLabel' ,{'-10','-9','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4',...
    '5','6','7','8','9','10'}) ;
xlabel('lag','Fontsize',12)
grid on;







