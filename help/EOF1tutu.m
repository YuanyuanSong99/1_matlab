clc;clear
lon=ncread('D:\data\ERA\hgt_daily_500_4time_2.5\hgt.daily.500.2.5_1979.nc','longitude');
lat=ncread('D:\data\ERA\hgt_daily_500_4time_2.5\hgt.daily.500.2.5_1979.nc','latitude');
x='D:\data\ERA\hgt_daily_500_4time_2.5\hgt.daily.500.2.5_';  %x、y用来构造文件名
y='.nc';
nt=39;nx=144;ny=37;nd=132;nt_nd=nt*nd;nt1=7;nt_nd1=nt1*nd;
for a=1:nt    %循环读取68年数据
    zz=1978+a;
    z=num2str(zz);
    c=[x z y];  %c：文件名字
    hgt_tem=ncread(c,'z'); % ，4维     %hgt_tem的意思是临时的hgt数据
    if rem(zz,4)==0  %闰年直接读出来，经度取出来
        hgt_tem1=reshape(hgt_tem,nx,ny,4,366);
        hgt_tem1=nanmean(hgt_tem1,3);
        hgt_tem1=reshape(hgt_tem1,nx,ny,366);
        hgthgt(:,:,:,a)=hgt_tem1;            %hgthgt是用来画图de数据，同上，4维
    else                                 %为保证数列可建，非闰年用nan帮忙，两种处理办法
        hgt_tem2=reshape(hgt_tem,nx,ny,4,365);
        hgt_tem2=nanmean(hgt_tem2,3);
        hgt_tem2=reshape(hgt_tem2,nx,ny,365);
        hgthgt(:,:,1:59,a)=hgt_tem2(:,:,1:59);
        hgthgt(:,:,60,a)=nan;
        hgthgt(:,:,61:366,a)=hgt_tem2(:,:,60:365);
    end
end
hgthgt=hgthgt(:,:,133:264,:)./9.8;
jp=nanmean(hgthgt,4);
for i=1:nt
    hgt(:,:,:,i)=hgthgt(:,:,:,i)-jp;
end
hgt1(:,:,:,1)=hgt(:,:,:,3);hgthgt1(:,:,:,1)=hgthgt(:,:,:,3);
hgt1(:,:,:,2)=hgt(:,:,:,10);hgthgt1(:,:,:,2)=hgthgt(:,:,:,10);
hgt1(:,:,:,3)=hgt(:,:,:,17);hgthgt1(:,:,:,3)=hgthgt(:,:,:,17);
hgt1(:,:,:,4)=hgt(:,:,:,20);hgthgt1(:,:,:,4)=hgthgt(:,:,:,20);
hgt1(:,:,:,5)=hgt(:,:,:,21);hgthgt1(:,:,:,5)=hgthgt(:,:,:,21);
hgt1(:,:,:,6)=hgt(:,:,:,28);hgthgt1(:,:,:,6)=hgthgt(:,:,:,28);
hgt1(:,:,:,7)=hgt(:,:,:,32);hgthgt1(:,:,:,7)=hgthgt(:,:,:,32);

for i=1:ny
    hgt2(:,i,:,:)=hgt1(:,ny+1-i,:,:);
    hgthgt2(:,i,:,:)=hgthgt1(:,ny+1-i,:,:);
    lat1(i)=lat(ny+1-i);
end
for i=1:nx
    for j=1:ny
        for k=1:nt1
            h(i,j,:,k)=smooth(hgthgt2(i,j,:,k));  
        end
    end
end
h=reshape(h,nx,ny,nt_nd1);
p=zeros(1,nt_nd1);
g=78;
for i=1:nt_nd1      % 对时间循环 
    for j=79:103     % 对经度循环   !!!!!!!!!!!!经度选的范围是15E-75E
        GHGN(j-g,1,i)=h(j,35,i)-h(j,27,i);   % q1:dita为5时de位势高度差，用于第二个条件
        GHGN(j-g,2,i)=h(j,33,i)-h(j,25,i);   % q2:dita为0时de位势高度差，用于第二个条件
        GHGN(j-g,3,i)=h(j,31,i)-h(j,23,i); % q3:dita为-5时de位势高度差，用于第二个条件
        GHGS(j-g,1,i)=h(j,27,i)-h(j,19,i);   % q1:dita为5时de位势高度差，用于第二个条件
        GHGS(j-g,2,i)=h(j,25,i)-h(j,17,i);   % q2:dita为0时de位势高度差，用于第二个条件
        GHGS(j-g,3,i)=h(j,23,i)-h(j,15,i); % q3:dita为-5时de位势高度差，用于第二个条件
        if GHGS(j-g,1,i)>0&GHGN(j-g,1,i)<-200      %deta为5时同时满足两个条件
            p(i)=1;
        elseif GHGS(j-g,2,i)>0&GHGN(j-g,2,i)<-200  %deta为0时同时满足两个条件
            p(i)=1;
        elseif GHGS(j-g,3,i)>0&GHGN(j-g,3,i)<-200  %deta为-5时同时满足两个条件
            p(i)=1;
        end                               %任意一个if成立时，p向量中对应位置会变成1，代表这条经度上这一天有阻塞发生
    end
end
p=reshape(p,nd,nt1);
for i=1:nt1    
    for j=21:112
        if p(j,i)==1                                    %在该天发生阻塞的情况下
            if p(j-1,i)==0&p(j+1,i)==1&p(j+2,i)==1        %这天为阻塞第1天满足de条件
                p(j,i)=1;
            elseif p(j+1,i)==0&p(j-1,i)==1&p(j-2,i)==1      %这天为阻塞最后一天满足de条件
                p(j,i)=1;
            elseif p(j-1,i)==1&p(j+1,i)==1                %这天为阻塞非头非尾时满足de条件
                p(j,i)=1;
            else
                p(j,i)=0;
            end
        end
    end
end
    p(1:20,:)=0;p(113:132,:)=0;
    p(21:33,3)=0;
    p(104:112,1)=0;p(111:112,3)=0;
    %p就是最后发生阻塞的天数，其中，夏季之外的前后20天没有做判断，全部标为0
    p=reshape(p,nt_nd1,1);
    k1=0;k2=0;
    for i=2:nt_nd1-1
        if p(i-1)==0&p(i)==1&p(i+1)==1
            p(i)=2;
            k1=k1+1;
            kaishi(k1)=i;
        elseif p(i-1)==1&p(i)==1&p(i+1)==0
            p(i)=3;
            k2=k2+1;
            jieshu(k2)=i;
        end
    end
    zhouqi=jieshu-kaishi+1;
    ghgs=max(max(GHGS,[],1),[],2);
    for i=1:k1
        [m1,m2]=max(ghgs(kaishi(i):jieshu(i)));
        lag0(i)=kaishi(i)+m2-1;
    end
    hgt2=reshape(hgt2,nx,ny,nt_nd1);
        for j=1:k1
            hgt_lag(:,:,21,j)=hgt2(:,:,lag0(j));
            for k=1:20
                hgt_lag(:,:,k,j)=hgt2(:,:,lag0(j)-21+k);
                hgt_lag(:,:,k+21,j)=hgt2(:,:,lag0(j)+k);
            end
        end
    hgt_lag1=nanmean(hgt_lag,4);
    fid=fopen('f:\hgt_lag.dat','wb');
    fwrite(fid,hgt_lag1,'float')
    fclose(fid);
    
    %%%%%%%%%%%%t检验
    
    for i=1:nx
        for j=1:ny
            for k=1:41
                t_tem=reshape(hgt_lag(i,j,k,:),1,k1);
                t_hgt(i,j,k)=nanmean(t_tem)/std(t_tem)*sqrt(k1);
                %if t_hgt(i,j,k)<=2.07
                 %   hgt_lag1(i,j,k)=nan;
                %end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lon=ncread('D:\data\ERA\t2m_daily_4time_2.5\tem.daily.1979.nc','longitude');
lat=ncread('D:\data\ERA\t2m_daily_4time_2.5\tem.daily.1979.nc','latitude');
x='D:\data\ERA\t2m_daily_4time_2.5\tem.daily.';  %x、y用来构造文件名
y='.nc';
for a=1:nt    %循环读取68年数据
    zz=1978+a;
    z=num2str(zz);
    c=[x z y];  %c：文件名字
    air_tem=ncread(c,'t2m'); % ，4维     %hgt_tem的意思是临时的hgt数据
    if rem(zz,4)==0  %闰年直接读出来，经度取出来
        air_tem1=reshape(air_tem,nx,ny,4,366);
        air_tem1=air_tem1(:,:,3,:);
        air_tem1=reshape(air_tem1,nx,ny,366);
        airair(:,:,:,a)=air_tem1;            %hgthgt是用来画图de数据，同上，4维
    else                                 %为保证数列可建，非闰年用nan帮忙，两种处理办法
        air_tem2=reshape(air_tem,nx,ny,4,365);
        air_tem2=air_tem2(:,:,3,:);
        air_tem2=reshape(air_tem2,nx,ny,365);
        airair(:,:,1:59,a)=air_tem2(:,:,1:59);
        airair(:,:,60,a)=nan;
        airair(:,:,61:366,a)=air_tem2(:,:,60:365);
    end
end
airair=airair(:,:,133:264,:);
for i=1:nx
    for j=1:ny
        for k=1:nd
            air(i,j,k,:)=detrend(reshape(airair(i,j,k,:),1,nt));
        end
    end
end
clear airair
air1(:,:,:,1)=air(:,:,:,3);
air1(:,:,:,2)=air(:,:,:,10);
air1(:,:,:,3)=air(:,:,:,17);
air1(:,:,:,4)=air(:,:,:,20);
air1(:,:,:,5)=air(:,:,:,21);
air1(:,:,:,6)=air(:,:,:,28);
air1(:,:,:,7)=air(:,:,:,32);
for i=1:ny
    air2(:,i,:,:)=air1(:,ny+1-i,:,:);
end
airair=reshape(air2,nx,ny,nt_nd1);
for j=1:k1
    air_lag(:,:,21,j)=airair(:,:,lag0(j));
    for k=1:20
        air_lag(:,:,k,j)=airair(:,:,lag0(j)-21+k);
        air_lag(:,:,k+21,j)=airair(:,:,lag0(j)+k);
    end
end
air_lag1=nanmean(air_lag,4);
for i=1:nx
        for j=1:ny
            for k=1:41
                t_tem1=reshape(air_lag(i,j,k,:),1,k1);
                t_air(i,j,k)=nanmean(t_tem1)/std(t_tem1)*sqrt(k1);
                if abs(t_air(i,j,k))<=1.71
                    air_lag1(i,j,k)=nan;
                end
            end
        end
end
for i=1:21
    m_proj('equidistant Cylindrical','lat',[20 90],'long',[-150 150]);
    m_coast('color',[0.5 0.5 0.5]);
    m_grid('xtick',[-150:60:150],'ytick',[30:30:90],'fontsize',8);
    hold on
    lon=[-150:2.5:150];lat=[20:2.5:90];
    m_contourf(lon,lat,air_lag1(13:133,9:37,i)',[-3:0.01:3],'linestyle','none');
    hold on

    %%%%%%

    caxis([-5.5 5.5])
    imread('f:\yeah\colorbar\colorbar1.png');
    color=ans((end)/2,:,:);                
    colorfinal=squeeze(color);
    colorfinal=flipud(colorfinal);
    colormap(double(colorfinal)/255)

    %%%%%%
            
    [cs,h]=m_contour(lon,lat,hgt_lag1(13:133,9:37,i)',[10:10:200],'-k','linewidth',0.08)
    clabel(cs,h,'fontsize',5,'LabelSpacing',10000);
    hold on
    [cs,h]=m_contour(lon,lat,hgt_lag1(13:133,9:37,i)',[-205:10:0],'--k','linewidth',0.08)
    clabel(cs,h,'fontsize',5,'LabelSpacing',10000);
    hold on
    %set(gca,'position',[0.5 0.15 0.4 0.2])
    ti=['lag',num2str(i-21)];
    title(ti) 
    saveas(gcf,['F:\yeah\阻塞前期\EOF1\',ti,'.png']);
    close
end





hgt_mean=nanmean(hgt_lag1(:,:,1:11),3);
air_mean=reshape(nanmean(air_lag(:,:,1:11,:),3),nx,ny,k1);
air_mean1=nanmean(air_mean,3);
for i=1:nx
    for j=1:ny
            t_tem_mean=reshape(air_mean(i,j,:),1,k1);
            t_air_mean(i,j)=nanmean(t_tem_mean)/std(t_tem_mean)*sqrt(k1);
                if abs(t_air_mean(i,j))<=1.71
                    air_mean1(i,j)=nan;
                end
    end
end
hgt_meanmean=nanmean(hgt_lag1(:,:,16:26),3);
air_meanmean=reshape(nanmean(air_lag(:,:,16:26,:),3),nx,ny,k1);
air_meanmean1=nanmean(air_meanmean,3);
for i=1:nx
    for j=1:ny
            t_tem_meanmean=reshape(air_meanmean(i,j,:),1,k1);
            t_air_meanmean(i,j)=nanmean(t_tem_meanmean)/std(t_tem_meanmean)*sqrt(k1);
                if abs(t_air_meanmean(i,j))<=1.71
                    air_meanmean1(i,j)=nan;
                end
    end
end
subplot(1,2,1)
    m_proj('equidistant Cylindrical','lat',[20 90],'long',[-150 150]);
    m_coast('color',[0.5 0.5 0.5]);
    m_grid('xtick',[-150:60:150],'ytick',[30:30:90],'fontsize',8);
    hold on
    lon=[-150:2.5:150];lat=[20:2.5:90];
    m_contourf(lon,lat,air_mean1(13:133,9:37)',[-3:0.01:3],'linestyle','none');
    hold on

    %%%%%%

    caxis([-3.3 3.3])
    imread('f:\yeah\colorbar\colorbar1.png');
    color=ans((end)/2,:,:);                
    colorfinal=squeeze(color);
    colorfinal=flipud(colorfinal);
    colormap(double(colorfinal)/255)

    %%%%%%
            
    [cs,h]=m_contour(lon,lat,hgt_mean(13:133,9:37)',[5:5:60],'-k','linewidth',0.08)
    clabel(cs,h,'fontsize',5,'LabelSpacing',10000);
    hold on
    [cs,h]=m_contour(lon,lat,hgt_mean(13:133,9:37)',[-60:5:-5],'--k','linewidth',0.08)
    clabel(cs,h,'fontsize',5,'LabelSpacing',10000);
    hold on
    %set(gca,'position',[0.5 0.15 0.4 0.2])
  
    title('lag-20 to lag-10')


subplot(1,2,2)
    m_proj('equidistant Cylindrical','lat',[20 90],'long',[-150 150]);
    m_coast('color',[0.5 0.5 0.5]);
    m_grid('xtick',[-150:60:150],'ytick',[30:30:90],'fontsize',8);
    hold on
    lon=[-150:2.5:150];lat=[20:2.5:90];
    m_contourf(lon,lat,air_meanmean1(13:133,9:37)',[-3:0.01:3],'linestyle','none');
    hold on

    %%%%%%

    caxis([-3.3 3.3])
    imread('f:\yeah\colorbar\colorbar1.png');
    color=ans((end)/2,:,:);                
    colorfinal=squeeze(color);
    colorfinal=flipud(colorfinal);
    colormap(double(colorfinal)/255)

    %%%%%%
            
    [cs,h]=m_contour(lon,lat,hgt_meanmean(13:133,9:37)',[5:5:60],'-k','linewidth',0.08)
    clabel(cs,h,'fontsize',5,'LabelSpacing',10000);
    hold on
    [cs,h]=m_contour(lon,lat,hgt_meanmean(13:133,9:37)',[-60:5:-5],'--k','linewidth',0.08)
    clabel(cs,h,'fontsize',5,'LabelSpacing',10000);
    hold on
    %set(gca,'position',[0.5 0.15 0.4 0.2])
  
    title('lag-5 to lag-5')
    colorbar('ytick',[-3 -2.4 -1.8 -1.2 -0.6 0 0.6 1.2 1.8 2.4 3],'location','southoutside','position',[0.39 0.3 0.25 0.02]);

%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 画图

for i=1:10
    subplot(5,2,i);
    m_proj('equidistant Cylindrical','lat',[20 90],'long',[-150 150]);
    m_coast('color',[0.5 0.5 0.5]);
    m_grid('xtick',[-150:60:150],'ytick',[30:30:90],'fontsize',8);
    hold on
    lon=[-150:2.5:150];lat=[20:2.5:90];
    m_contourf(lon,lat,air_lag1(13:133,9:37,i+10)',[-3:0.01:3],'linestyle','none');
    hold on

   %%%%%%
    caxis([-5.5 5.5])
   imread('f:\yeah\colorbar\colorbar1.png');
    color=ans((end)/2,:,:);                
    colorfinal=squeeze(color);
    colorfinal=flipud(colorfinal);
    colormap(double(colorfinal)/255)
    %%%%%%
            
    subplot(5,2,i)
    [cs,h]=m_contour(lon,lat,hgt_lag1(13:133,9:37,i+10)',[10:10:200],'-k','linewidth',0.08)
    clabel(cs,h,'fontsize',5,'LabelSpacing',10000);
    hold on
    [cs,h]=m_contour(lon,lat,hgt_lag1(13:133,9:37,i+10)',[-200:10:-10],'--k','linewidth',0.08)
    clabel(cs,h,'fontsize',5,'LabelSpacing',10000);
    hold on
    
   ti=['lag',num2str(i-11)];
    title(ti)
end
caxis([-5.5 5.5])
imread('f:\yeah\colorbar\colorbar1.png');
color=ans((end)/2,:,:);                
colorfinal=squeeze(color);
colorfinal=flipud(colorfinal);
colormap(double(colorfinal)/255)
colorbar('ytick',[-5 -4 -3 -2 -1 0 1 2 3 4 5],'location','southoutside','position',[0.425 0.04 0.2 0.02]);
subplot(5,2,10)
set(gca,'position',[0.55 0.08 0.15 0.15])
subplot(5,2,8)
set(gca,'position',[0.55 0.25 0.15 0.15])
subplot(5,2,6)
set(gca,'position',[0.55 0.42 0.15 0.15])
subplot(5,2,4)
set(gca,'position',[0.55 0.59 0.15 0.15])
subplot(5,2,2)
set(gca,'position',[0.55 0.76 0.15 0.15])
subplot(5,2,9)
set(gca,'position',[0.35 0.08 0.15 0.15])
subplot(5,2,7)
set(gca,'position',[0.35 0.25 0.15 0.15])
subplot(5,2,5)
set(gca,'position',[0.35 0.42 0.15 0.15])
subplot(5,2,3)
set(gca,'position',[0.35 0.59 0.15 0.15])
subplot(5,2,1)
set(gca,'position',[0.35 0.76 0.15 0.15])