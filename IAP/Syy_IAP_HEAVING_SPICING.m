%% spicing heaving
clear;clc;close all;
addpath(genpath('G:\1_matlab\help'));
% data_dir='J:\IAP\data\';
% data_savedir='../../data/HEAVING_SPICING/';
%% Tem
datadir='G:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
% Salt
datadir='G:\data\IAP\salinity\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
k=length(filelist);
for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Salt(:,:,:,s) = ncread(filename1,'absolute salinity'); %读入变量
end

%%
salt_iap = Salt; temp_iap = Temp;
yrs=1958; yre=2019;
% load([data_dir,'IAP_salt_globe_upper2000_1958_2019.mat'])
% load([data_dir,'IAP_temp_globe_upper2000_1958_2019.mat'])
h=depthData;
SP=permute(salt_iap,[2 1 3 4]);tt=permute(temp_iap,[2 1 3 4]);
[ny nx nz nt]=size(SP); % [longitude latitude depth time]
lon=lonData; lat=latData;
clear temp_iap salt_iap

% %% 计算密度
for  i=1:nz
    i
    for j=1:nt;
        h_1=h(i).*ones(1,ny); 
        P=sw_pres(h_1,lat')';  % 海水压强
        CT(:,:,i,j)= gsw_CT_from_t(squeeze(SP(:,:,i,j)),squeeze(tt(:,:,i,j)),P);%Conservative Temperature
        pt(:,:,i,j) = gsw_pt_from_CT(squeeze(SP(:,:,i,j)),squeeze(CT(:,:,i,j)));%potential temperature
        sigma(:,:,i,j)=gsw_sigma0(squeeze(SP(:,:,i,j)),squeeze(CT(:,:,i,j)));
    end
end
%
sigma_new1=[22:0.1:26];
sigma_new2=[26.01:0.01:27.75];
sigma_new=[sigma_new1,sigma_new2];

% %% interp on sigma grid
for i=1:ny
    i
    for j=1:nx
        for k=1:nt;
            A=squeeze(sigma(i,j,:,k));
            B=squeeze(pt(i,j,:,k));
            if    all(~isnan(A))
                [A1,pos]=sort(A);%%排序
                B1=B(pos);
                %去掉相同元素
                [A1,X1,~]=unique(A1);
                B1=B1(X1);
                t_1(1:length(sigma_new),i,j,k)=interp1(A1,B1,sigma_new);
            else
                a=isnan(A);
                b=find(a==1); %标记是nan的元素
                if b(1)==1||  b(1)==2  %第一第二个是nan
                    t_1(1:length(sigma_new),i,j,k)=nan;
                else
                    A1=A(1:b(1)-1);
                    B1=B(1:b(1)-1);
                    [A1_1,pos]=sort(A1);%%排序
                    B1_1=B1(pos);
                    %去掉相同元素
                    [A1_1,X1,~]=unique(A1_1);
                    if length(A1_1)==1
                        t_1(1:length(sigma_new),i,j,k)=nan;
                    else
                        B1_1=B1_1(X1);
                        t_1(1:length(sigma_new),i,j,k)=interp1(A1_1,B1_1,sigma_new);
                    end
                end
            end
        end
    end
end


for i=1:ny
    i
    for j=1:nx
        for k=1:nt
            A=squeeze(sigma(i,j,:,k));
            B=h;
            if    all(~isnan(A))
                [A1,pos]=sort(A);%%排序
                B1=B(pos);
                %去掉相同元素
                [A1,X1,~]=unique(A1);
                B1=B1(X1);
                h_1(1:length(sigma_new),i,j,k)=interp1(A1,B1,sigma_new);
            else
                a=isnan(A);
                b=find(a==1); %标记是nan的元素
                if b(1)==1||  b(1)==2  %第一第二个是nan
                    h_1(1:length(sigma_new),i,j,k)=nan;
                else
                    A1=A(1:b(1)-1);
                    B1=B(1:b(1)-1);
                    [A1_1,pos]=sort(A1);%%排序
                    B1_1=B1(pos);
                    %去掉相同元素
                    [A1_1,X1,~]=unique(A1_1);
                    if length(A1_1)==1
                        h_1(1:length(sigma_new),i,j,k)=nan;
                    else
                        B1_1=B1_1(X1);
                        h_1(1:length(sigma_new),i,j,k)=interp1(A1_1,B1_1,sigma_new);
                    end
                end
            end
        end
    end
end
%% 
%============= heaving =================
dT_dt_sigma=t_1(:,:,:,2:end)-t_1(:,:,:,1:end-1);

% clear t_1
for i=1:nt-1
    i
    dT_dt_sigma_sum(:,:,:,i)=nansum(dT_dt_sigma(:,:,:,1:i),4); %
end
% clear dT_dt_sigma
for z=1:nt-1;
    z
    dT_dt_sigma_sum1(:,:,:,z)=dT_dt_sigma_sum(:,:,:,z)-nanmean(dT_dt_sigma_sum,4);
end
% clear dT_dt_sigma_sum

h_1_cli=squeeze(mean(h_1,4));
t_1_mean=nanmean(t_1,4);
for z=1:size(t_1_mean,1);
    z
    if z==1;
        dt_dz_mean_sigma(z,:,:)=(t_1_mean(2,:,:)-t_1_mean(1,:,:))/(sigma_new(2)-sigma_new(1));
    end
    if z==size(t_1_mean,1);
        dt_dz_mean_sigma(z,:,:)=(t_1_mean(end,:,:)-t_1_mean(end-1,:,:))/(sigma_new(end)-sigma_new(end-1));
    end
    if z>1 && z<size(t_1_mean,1);
        dt_dz_mean_sigma(z,:,:)=(t_1_mean(z+1,:,:)-t_1_mean(z-1,:,:))/(sigma_new(z+1)-sigma_new(z-1));
    end
end
%% ======= interp on Z grid ========
dT_dz_cli=zeros(ny,nx,nz);
for i=1:ny
    i
    for j=1:nx
        A=h_1_cli(:,i,j);
        B=dt_dz_mean_sigma(:,i,j);
        %%都不是nan
        if    all(~isnan(A))
            [A1,pos]=sort(A);%%排序
            B1=B(pos);
            %去掉相同元素
            [A1,X1,~]=unique(A1);
            B1=B1(X1);
            dT_dz_cli(i,j,1:nz)=interp1(A1,B1,h);
            % 其中有nan
        else
            a=isnan(A);
            b=find(a==0); %标记是nan的元素
            if length(b)==0||  length(b)==1  %第一第二个是nan
                dT_dz_cli(i,j,1:nz)=nan;
            else
                A1=A(b);
                B1=B(b);
                [A1_1,pos]=sort(A1);%%排序
                B1_1=B1(pos);
                %去掉相同元素
                [A1_1,X1,~]=unique(A1_1);
                if length(A1_1)==1
                    dT_dz_cli(i,j,1:nz)=nan;
                else
                    B1_1=B1_1(X1);
                    dT_dz_cli(i,j,1:nz)=interp1(A1_1,B1_1,h);
                end
            end
        end
    end
end

for t=1:size(sigma,4);
    if t==1;
        dsigma_dt(:,:,:,t)=sigma(:,:,:,2)-sigma(:,:,:,1);
    end
    if t==size(sigma,4);
        dsigma_dt(:,:,:,t)=sigma(:,:,:,end)-sigma(:,:,:,end-1);
    end
    if t>1 && t<size(sigma,4);
        dsigma_dt(:,:,:,t)=(sigma(:,:,:,t+1)-sigma(:,:,:,t-1))/2;
    end
end

for i=1:size(dsigma_dt,4);
    dsigma_dt_sum(:,:,:,i)=nansum(dsigma_dt(:,:,:,1:i),4);
end
for i=1:size(dsigma_dt,4);
    dsigma_dt_new(:,:,:,i)=dsigma_dt_sum(:,:,:,i)-nanmean(dsigma_dt_sum,4);
end
for t=1:size(dsigma_dt,4);
    heaving_T_durack(:,:,:,t)=squeeze(dsigma_dt_new(:,:,:,t)).*dT_dz_cli;
end
%%
clearvars -except heaving_T_durack t_1 ny nx nz nt h_1 h
%============= spicing =================
for t=1:size(t_1,4);
    if t==1;
        dT_dt_sigmad(:,:,:,t)=t_1(:,:,:,2)-t_1(:,:,:,1);
    end
    if t==size(t_1,4);
        dT_dt_sigmad(:,:,:,t)=t_1(:,:,:,end)-t_1(:,:,:,end-1);
    end
    if t>1 && t<size(t_1,4);
        dT_dt_sigmad(:,:,:,t)=(t_1(:,:,:,t+1)-t_1(:,:,:,t-1))/2;
    end
end
for i=1:size(dT_dt_sigmad,4)
    dT_dt_sigma_sumd(:,:,:,i)=nansum(dT_dt_sigmad(:,:,:,1:i),4); 
end
for z=1:size(dT_dt_sigmad,4)
    dT_dt_sigma_sumd1(:,:,:,z)=dT_dt_sigma_sumd(:,:,:,z)-nanmean(dT_dt_sigma_sumd,4);
end

%% ======= interp on Z grid ========
spicing_T_durack=zeros(ny,nx,nz,nt);
for i=1:ny
    i
    for j=1:nx
        for k=1:nt
            A=squeeze(h_1(:,i,j,k));
            B=squeeze(dT_dt_sigma_sumd1(:,i,j,k));
            %%都不是nan
            if    all(~isnan(A))
                [A1,pos]=sort(A);%%排序
                B1=B(pos);
                %去掉相同元素
                [A1,X1,~]=unique(A1);
                B1=B1(X1);
                spicing_T_durack(i,j,1:nz,k)=interp1(A1,B1,h);
                % 其中有nan
            else
                a=isnan(A);
                b=find(a==0); %标记是nan的元素
                if length(b)==0||  length(b)==1  %第一第二个是nan
                    spicing_T_durack(i,j,1:nz,k)=nan;
                else
                    A1=A(b);
                    B1=B(b);
                    [A1_1,pos]=sort(A1);%%排序
                    B1_1=B1(pos);
                    %去掉相同元素
                    [A1_1,X1,~]=unique(A1_1);
                    if length(A1_1)==1
                        spicing_T_durack(i,j,1:nz,k)=nan;
                    else
                        B1_1=B1_1(X1);
                        spicing_T_durack(i,j,1:nz,k)=interp1(A1_1,B1_1,h);
                    end
                end
            end
        end
    end
end
save('G:/1_matlab/IAP/MatFile/heaving_T_durack.mat','heaving_T_durack');
save('G:/1_matlab/IAP/MatFile/spicing_T_durack.mat','spicing_T_durack');

