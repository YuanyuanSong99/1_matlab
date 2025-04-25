clc;
clear all;
close all;
tem=load('H:\som\contribution\data\tem_1958_2016_detrended.mat');
tem=tem.hgt;
lon=load('H:\som\contribution\data\lon.mat','lon');
lon=lon.lon;
lat=load('H:\som\contribution\data\lat.mat','lat');
lat=lat.lat;
tem=tem(21:71,181:331,:);
count1=1;
for ii=1:5:51
    count2=1;
    for jj=1:5:151
        tem1(count1,count2,:)=tem(ii,jj,:);
        count2=count2+1;
    end
    count1=count1+1;
end
tem1=reshape(tem1,31*11,5325);

% for ii=1:31*11
%     temp=tem1(ii,:);
%     temp=zscore(temp);
%     temp(temp>mean(temp)-std(temp)/2)=0;
%     tem1(ii,:)=temp;
% end

[S2,hits,sM,sD,bmus,Qe,Te]=SOM_func(tem1,2,3,6,20);
bmus1=bmus(:,20);
save('H:\som\contribution\data\bmus1.mat','bmus1');  % 保存到其他文件夹的写法

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=1;
f2=1;
f3=1;
f4=1;
f5=1;
f6=1;
tem=reshape(tem,51,151,5325);
for ii=1:5325
    if bmus1(ii)==1
        t2m1(:,:,f1)=tem(:,:,ii);
        f1=f1+1;
    end
    if bmus1(ii)==2
        t2m2(:,:,f2)=tem(:,:,ii);
        f2=f2+1;
    end
    if bmus1(ii)==3
        t2m3(:,:,f3)=tem(:,:,ii);
        f3=f3+1;
    end
    if bmus1(ii)==4
        t2m4(:,:,f4)=tem(:,:,ii);
        f4=f4+1;
    end
    if bmus1(ii)==5
        t2m5(:,:,f5)=tem(:,:,ii);
        f5=f5+1;
    end
    if bmus1(ii)==6
        t2m6(:,:,f6)=tem(:,:,ii);
        f6=f6+1;
    end
end
t2m_m(:,:,1)=mean(t2m1,3);
t2m_m(:,:,2)=mean(t2m2,3);
t2m_m(:,:,3)=mean(t2m3,3);
t2m_m(:,:,4)=mean(t2m4,3);
t2m_m(:,:,5)=mean(t2m5,3);
t2m_m(:,:,6)=mean(t2m6,3);
save('H:\som\contribution\data\tem_som_result.mat','t2m_m');