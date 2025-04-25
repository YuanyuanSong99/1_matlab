clc;clear;clear all memory;
cro=5;
fid=fopen('C:\new\1wind_sst_sic\daily_line\U10_daily_bks_bk.txt','r');
u10=fscanf(fid,'%g',[14000,cro]);
fclose(fid);
clear fid
datay=zeros(366,39,cro);
num=0;
for i=1979:2016
    if mod(i,4)==0
        datay(1:366,i-1978,1:cro)=u10(num+1:num+366,1:cro);
        num=num+366;
    else
        datay(1:59,i-1978,1:cro)=u10(num+1:num+59,1:cro);
        datay(60,i-1978,1:cro)=nan;
        datay(61:366,i-1978,1:cro)=u10(num+60:num+365,1:cro);
        num=num+365;
    end
end
datay(1:59,39,1:cro)=u10(num+1:num+59,1:cro);
datay(60,39,1:cro)=nan;
datay(61:121,39,1:cro)=u10(num+60:num+120,1:cro);
datay(122:366,39,1:cro)=nan;
clear u10 num
datay=datay-repmat(nanmean(datay,2),1,39,1);
datay=reshape(datay,[366*39,cro]);
datay(isnan(datay(:,1))==1,:)=[];
u10=zeros(size(datay));
for i=1:cro
    y(1:14000)=datay(:,i);
    u10(:,i)=detrend(y);
end
clear datay i y

fid=fopen('C:\new\1wind_sst_sic\daily_line\V10_daily_bks_bk.txt','r');
v10=fscanf(fid,'%g',[14000,cro]);
fclose(fid);
clear fid
datay=zeros(366,39,cro);
num=0;
for i=1979:2016
    if mod(i,4)==0
        datay(1:366,i-1978,1:cro)=v10(num+1:num+366,1:cro);
        num=num+366;
    else
        datay(1:59,i-1978,1:cro)=v10(num+1:num+59,1:cro);
        datay(60,i-1978,1:cro)=nan;
        datay(61:366,i-1978,1:cro)=v10(num+60:num+365,1:cro);
        num=num+365;
    end
end
datay(1:59,39,1:cro)=v10(num+1:num+59,1:cro);
datay(60,39,1:cro)=nan;
datay(61:121,39,1:cro)=v10(num+60:num+120,1:cro);
datay(122:366,39,1:cro)=nan;
clear v10 num
datay=datay-repmat(nanmean(datay,2),1,39,1);
datay=reshape(datay,[366*39,cro]);
datay(isnan(datay(:,1))==1,:)=[];
v10=zeros(size(datay));
for i=1:cro
    y(1:14000)=datay(:,i);
    v10(:,i)=detrend(y);
end
clear datay i y

fid=fopen('C:\new\blocking\UB_1202_1004_2.txt','r');
ub=fscanf(fid,'%g',[38,213]);
fclose(fid);
clear fid
ub=ub';
ub=[zeros(153,38);ub];
ub=reshape(ub,[366*38,1]);
ub=[zeros(120,1);ub];
ub(isnan(ub))=[];
[Y,I]=sort(ub,'descend');
event=50;

inde1=[u10(:,1),v10(:,1),v10(:,4)];
inde=zeros(size(inde1));
for i=1:3
    yy(1:14000)=inde1(1:14000,i);
    inde(1:14000,i)=smooth(yy,5);
end
clear yy inde1
data=zeros(11,3,event);
for i=1:event
    k=I(i);
    data(1:11,1:3,i)=inde(k-5:k+5,1:3);
end
datas(1:3,1:50)=nanmean(data);

[Y1,I1]=sort(datas(1,:),'descend');
[Y2,I2]=sort(datas(2,:),'descend');
% [Y3,I3]=sort(datas(3,:),'descend');
for i=1:25
    num1=(i-1)*11+1;
    num2=i*11;
    k=I(I1(i));
    outdata1(num1:num2)=k-5:k+5;
    k=I(I2(i));
    outdata3(num1:num2)=k-5:k+5;
end
for i=26:50
    num1=(i-26)*11+1;
    num2=(i-25)*11;
    k=I(I1(i));
    outdata2(num1:num2)=k-5:k+5;
    k=I(I2(i));
    outdata4(num1:num2)=k-5:k+5;
end

clear u10 v10
%%
filename1='C:\new\data\ERAdata\ERA_interim_1.0_1979_2017430_seaicecover_4times.nc';
lon=ncread(filename1,'longitude');
lat=ncread(filename1,'latitude');
lon=[lon;lon(1)];
lat=lat(1:61);
data=zeros(360,61,14000);
for i=1:14000
    num=(i-1)*4+1;
    data1=ncread(filename1,'ci',[1,1,num],[360,61,4]);
    data(:,:,i)=mean(data1,3);
    clear data1;
end
clear filename1
dataa=zeros(360,61,366,39);
num=0;
for i=1979:2016
    if mod(i,4)==0
        dataa(1:360,1:61,1:366,i-1978)=data(1:360,1:61,num+1:num+366);
        num=num+366;
    else
        dataa(1:360,1:61,1:59,i-1978)=data(1:360,1:61,num+1:num+59);
        dataa(1:360,1:61,60,i-1978)=nan;
        dataa(1:360,1:61,61:366,i-1978)=data(1:360,1:61,num+60:num+365);
        num=num+365;
    end
end
dataa(1:360,1:61,1:59,39)=data(1:360,1:61,num+1:num+59);
dataa(1:360,1:61,60,39)=nan;
dataa(1:360,1:61,61:121,39)=data(1:360,1:61,num+60:num+120);
dataa(1:360,1:61,122:366,39)=nan;
clear data
data1=dataa-repmat(nanmean(dataa,4),1,1,1,39);
data1=reshape(data1,[360,61,366*39]);
data1(:,:,isnan(data1(1,1,:))==1)=[];
datat=zeros(360,61,14000);
for i=1:360
    for j=1:61
        y(1:14000)=data1(i,j,:);
        datat(i,j,:)=detrend(y);
    end
end
clear i j y num
clear dataa data1

p=zeros(360,61,4);
f=zeros(360,61,4);
[~,f(:,:,1)]=composed_student( datat, datat(:,:,outdata1) );
[~,f(:,:,2)]=composed_student( datat, datat(:,:,outdata2) );
[~,f(:,:,3)]=composed_student( datat, datat(:,:,outdata3) );
[~,f(:,:,4)]=composed_student( datat, datat(:,:,outdata4) );
p(:,:,1)=nanmean(datat(:,:,outdata1),3);
p(:,:,2)=nanmean(datat(:,:,outdata2),3);
p(:,:,3)=nanmean(datat(:,:,outdata3),3);
p(:,:,4)=nanmean(datat(:,:,outdata4),3);

% for i=1:360
%     for j=1:61
%         for k=1:4
%             if (f(i,j,k)==0)
%                 p(i,j,k)=0;
%             end
%         end
%     end
% end
% clear i j k

[lon1,lat1]=meshgrid(lon,lat);
color=[0.142,0.0,0.85;0.097,0.112,0.97;0.16,0.342,1.0;0.24,0.531,1.0;0.34,0.692,1.0;...
0.46,0.829,1.0;0.6,0.92,1.0;0.74,0.978,1.0;1.0,1.0,1.0;1.0,1.0,1.0;1.0,0.948,0.74;1.0,0.84,0.6;...
1.0,0.676,0.46;1.0,0.472,0.34;1.0,0.24,0.24;0.97,0.155,0.21;0.85,0.085,0.187;0.65,0.0,0.13];
c = interp1(linspace(0,1,size(color,1)),color,linspace(0,1,40),'pchip');
v1=(-40:1:40);
a=20;
fig11 = figure(1);
set(fig11,'Unit','centimeters')  
set(fig11,'Position',[0,0,50,26])  
axes1 = axes('Parent',fig11,'Unit','centimeters','Position',[2 15 10 10]);
datati=[p(:,:,1);p(1,:,1)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati*100,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)
hold on
for i=1:20:360
    for j=1:1:4
        if(f(i,j,1)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:10:360
    for j=5:1:10
        if(f(i,j,1)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:4:360
    for j=11:1:20
        if(f(i,j,1)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:3:360
    for j=21:1:40
        if(f(i,j,1)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end

axes2 = axes('Parent',fig11,'Unit','centimeters','Position',[14 15 10 10]);
datati=[p(:,:,2);p(1,:,2)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati*100,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)
hold on
for i=1:20:360
    for j=1:1:4
        if(f(i,j,2)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:10:360
    for j=5:1:10
        if(f(i,j,2)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:4:360
    for j=11:1:20
        if(f(i,j,2)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:3:360
    for j=21:1:40
        if(f(i,j,2)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end

axes3 = axes('Parent',fig11,'Unit','centimeters','Position',[26 15 10 10]);
datati=[p(:,:,3);p(1,:,3)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati*100,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)
hold on
for i=1:20:360
    for j=1:1:4
        if(f(i,j,3)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:10:360
    for j=5:1:10
        if(f(i,j,3)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:4:360
    for j=11:1:20
        if(f(i,j,3)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:3:360
    for j=21:1:40
        if(f(i,j,3)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end

axes4 = axes('Parent',fig11,'Unit','centimeters','Position',[38 15 10 10]);
datati=[p(:,:,4);p(1,:,4)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati*100,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)
hold on
for i=1:20:360
    for j=1:1:4
        if(f(i,j,4)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:10:360
    for j=5:1:10
        if(f(i,j,4)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:4:360
    for j=11:1:20
        if(f(i,j,4)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end
for i=1:3:360
    for j=21:1:40
        if(f(i,j,4)==1)
            m_plot(lon(i),lat(j),'.','markersize',4,'color','k');
        end
    end
end

clear datat datai p f i j
%%
filename1='C:\new\data\ERAdata\ERA_interim_1.0_1979_2017530_temp_4times.nc';
lon=ncread(filename1,'longitude');
lat=ncread(filename1,'latitude');
lon=[lon;lon(1)];
lat=lat(1:61);
data=zeros(360,61,14000);
for i=1:14000
    num=(i-1)*4+1;
    data1=ncread(filename1,'t2m',[1,1,num],[360,61,4]);
    data(:,:,i)=mean(data1,3);
    clear data1;
end
clear filename1
dataa=zeros(360,61,366,39);
num=0;
for i=1979:2016
    if mod(i,4)==0
        dataa(1:360,1:61,1:366,i-1978)=data(1:360,1:61,num+1:num+366);
        num=num+366;
    else
        dataa(1:360,1:61,1:59,i-1978)=data(1:360,1:61,num+1:num+59);
        dataa(1:360,1:61,60,i-1978)=nan;
        dataa(1:360,1:61,61:366,i-1978)=data(1:360,1:61,num+60:num+365);
        num=num+365;
    end
end
dataa(1:360,1:61,1:59,39)=data(1:360,1:61,num+1:num+59);
dataa(1:360,1:61,60,39)=nan;
dataa(1:360,1:61,61:121,39)=data(1:360,1:61,num+60:num+120);
dataa(1:360,1:61,122:366,39)=nan;
clear data
data1=dataa-repmat(nanmean(dataa,4),1,1,1,39);
data1=reshape(data1,[360,61,366*39]);
data1(:,:,isnan(data1(1,1,:))==1)=[];
datat=zeros(360,61,14000);
for i=1:360
    for j=1:61
        y(1:14000)=data1(i,j,:);
        datat(i,j,:)=detrend(y);
    end
end
clear i j y num
clear dataa data1

filename1='C:\new\data\ERAdata\ERA_interim_1.0_1979_2017430_500hgt_4times.nc';
data=zeros(360,61,14000);
for i=1:14000
    num=(i-1)*4+1;
    data1=ncread(filename1,'z',[1,1,num],[360,61,4]);
    data(:,:,i)=mean(data1,3);
    clear data1;
end
data=data/9.8;
clear filename1
dataa=zeros(360,61,366,39);
num=0;
for i=1979:2016
    if mod(i,4)==0
        dataa(1:360,1:61,1:366,i-1978)=data(1:360,1:61,num+1:num+366);
        num=num+366;
    else
        dataa(1:360,1:61,1:59,i-1978)=data(1:360,1:61,num+1:num+59);
        dataa(1:360,1:61,60,i-1978)=nan;
        dataa(1:360,1:61,61:366,i-1978)=data(1:360,1:61,num+60:num+365);
        num=num+365;
    end
end
dataa(1:360,1:61,1:59,39)=data(1:360,1:61,num+1:num+59);
dataa(1:360,1:61,60,39)=nan;
dataa(1:360,1:61,61:121,39)=data(1:360,1:61,num+60:num+120);
dataa(1:360,1:61,122:366,39)=nan;
clear data
data1=dataa-repmat(nanmean(dataa,4),1,1,1,39);
data1=reshape(data1,[360,61,366*39]);
data1(:,:,isnan(data1(1,1,:))==1)=[];
dataz=zeros(360,61,14000);
for i=1:360
    for j=1:61
        y(1:14000)=data1(i,j,:);
        dataz(i,j,:)=detrend(y);
    end
end
clear i j y num
clear dataa data1

p=zeros(360,61,8);
f=zeros(360,61,8);
[~,f(:,:,1)]=composed_student( datat, datat(:,:,outdata1) );
[~,f(:,:,2)]=composed_student( datat, datat(:,:,outdata2) );
[~,f(:,:,3)]=composed_student( dataz, dataz(:,:,outdata1) );
[~,f(:,:,4)]=composed_student( dataz, dataz(:,:,outdata2) );
[~,f(:,:,5)]=composed_student( datat, datat(:,:,outdata3) );
[~,f(:,:,6)]=composed_student( datat, datat(:,:,outdata4) );
[~,f(:,:,7)]=composed_student( dataz, dataz(:,:,outdata3) );
[~,f(:,:,8)]=composed_student( dataz, dataz(:,:,outdata4) );
p(:,:,1)=nanmean(datat(:,:,outdata1),3);
p(:,:,2)=nanmean(datat(:,:,outdata2),3);
p(:,:,3)=nanmean(dataz(:,:,outdata1),3);
p(:,:,4)=nanmean(dataz(:,:,outdata2),3);
p(:,:,5)=nanmean(datat(:,:,outdata3),3);
p(:,:,6)=nanmean(datat(:,:,outdata4),3);
p(:,:,7)=nanmean(dataz(:,:,outdata3),3);
p(:,:,8)=nanmean(dataz(:,:,outdata4),3);

for i=1:360
    for j=1:61
        for k=1:2
            if (f(i,j,k)==0)
                p(i,j,k)=0;
            end
        end
        for k=5:6
            if (f(i,j,k)==0)
                p(i,j,k)=0;
            end
        end
    end
end
clear i j k

[lon1,lat1]=meshgrid(lon,lat);
color=[0.142,0.0,0.85;0.097,0.112,0.97;0.16,0.342,1.0;0.24,0.531,1.0;0.34,0.692,1.0;...
0.46,0.829,1.0;0.6,0.92,1.0;0.74,0.978,1.0;1.0,1.0,1.0;1.0,1.0,1.0;1.0,0.948,0.74;1.0,0.84,0.6;...
1.0,0.676,0.46;1.0,0.472,0.34;1.0,0.24,0.24;0.97,0.155,0.21;0.85,0.085,0.187;0.65,0.0,0.13];
c = interp1(linspace(0,1,size(color,1)),color,linspace(0,1,40),'pchip');
v1=(-20:0.5:20);
v2=(30:30:300);
v3=(-300:30:-30);
a=10;
axes5 = axes('Parent',fig11,'Unit','centimeters','Position',[2 2 10 10]);
datati=[p(:,:,1);p(1,:,1)]';
m_proj('stereographic','lon',60,'lat',90,'radius',60,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
hold on
datazi=[p(:,:,3);p(1,:,3)]';
dataz1=datazi;dataz2=datazi;
dataz1(dataz1<0)=0;dataz2(dataz2>0)=0;
m_contour(lon1,lat1,dataz1,'LevelList',v2,'LineStyle','-','color','k','LineWidth',1.2);
hold on
m_contour(lon1,lat1,dataz2,'LevelList',v3,'LineStyle','-.','color','k','LineWidth',1.2);
m_coast('line','color','k');
m_grid('xtick',12,'ytick',7,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)
% text(-1.05,0.95,'(a)','fontsize',23)
% research_rangee([65,85],[30,90],'-g',1.5)  

axes6 = axes('Parent',fig11,'Unit','centimeters','Position',[14 2 10 10]);
datati=[p(:,:,2);p(1,:,2)]';
m_proj('stereographic','lon',60,'lat',90,'radius',60,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
hold on
datazi=[p(:,:,4);p(1,:,4)]';
dataz1=datazi;dataz2=datazi;
dataz1(dataz1<0)=0;dataz2(dataz2>0)=0;
m_contour(lon1,lat1,dataz1,'LevelList',v2,'LineStyle','-','color','k','LineWidth',1.2);
hold on
m_contour(lon1,lat1,dataz2,'LevelList',v3,'LineStyle','-.','color','k','LineWidth',1.2);
m_coast('line','color','k');
m_grid('xtick',12,'ytick',7,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)
% text(-1.05,0.95,'(b)','fontsize',23)

axes7 = axes('Parent',fig11,'Unit','centimeters','Position',[26 2 10 10]);
datati=[p(:,:,5);p(1,:,5)]';
m_proj('stereographic','lon',60,'lat',90,'radius',60,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
hold on
datazi=[p(:,:,7);p(1,:,7)]';
dataz1=datazi;dataz2=datazi;
dataz1(dataz1<0)=0;dataz2(dataz2>0)=0;
m_contour(lon1,lat1,dataz1,'LevelList',v2,'LineStyle','-','color','k','LineWidth',1.2);
hold on
m_contour(lon1,lat1,dataz2,'LevelList',v3,'LineStyle','-.','color','k','LineWidth',1.2);
m_coast('line','color','k');
m_grid('xtick',12,'ytick',7,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)
% text(-1.05,0.95,'(b)','fontsize',23)
% research_rangee([65,85],[30,90],'-g',1.5)  

axes8 = axes('Parent',fig11,'Unit','centimeters','Position',[38 2 10 10]);
datati=[p(:,:,6);p(1,:,6)]';
m_proj('stereographic','lon',60,'lat',90,'radius',60,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
hold on
datazi=[p(:,:,8);p(1,:,8)]';
dataz1=datazi;dataz2=datazi;
dataz1(dataz1<0)=0;dataz2(dataz2>0)=0;
m_contour(lon1,lat1,dataz1,'LevelList',v2,'LineStyle','-','color','k','LineWidth',1.2);
hold on
m_contour(lon1,lat1,dataz2,'LevelList',v3,'LineStyle','-.','color','k','LineWidth',1.2);
m_coast('line','color','k');
m_grid('xtick',12,'ytick',7,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

set(fig11,'PaperUnits','centimeters','PaperPosition',[0 0 50 26])
print -f1 -r600 -dpng com_ubuv11_1;  
clear fig11 axes1 axes2 axes3 axes4

%%
clc;clear;clear all memory;
cro=5;
fid=fopen('C:\new\1wind_sst_sic\daily_line\U10_daily_bks_bk.txt','r');
u10=fscanf(fid,'%g',[14000,cro]);
fclose(fid);
clear fid
datay=zeros(366,39,cro);
num=0;
for i=1979:2016
    if mod(i,4)==0
        datay(1:366,i-1978,1:cro)=u10(num+1:num+366,1:cro);
        num=num+366;
    else
        datay(1:59,i-1978,1:cro)=u10(num+1:num+59,1:cro);
        datay(60,i-1978,1:cro)=nan;
        datay(61:366,i-1978,1:cro)=u10(num+60:num+365,1:cro);
        num=num+365;
    end
end
datay(1:59,39,1:cro)=u10(num+1:num+59,1:cro);
datay(60,39,1:cro)=nan;
datay(61:121,39,1:cro)=u10(num+60:num+120,1:cro);
datay(122:366,39,1:cro)=nan;
clear u10 num
datay=datay-repmat(nanmean(datay,2),1,39,1);
datay=reshape(datay,[366*39,cro]);
datay(isnan(datay(:,1))==1,:)=[];
u10=zeros(size(datay));
for i=1:cro
    y(1:14000)=datay(:,i);
    u10(:,i)=detrend(y);
end
clear datay i y

fid=fopen('C:\new\1wind_sst_sic\daily_line\V10_daily_bks_bk.txt','r');
v10=fscanf(fid,'%g',[14000,cro]);
fclose(fid);
clear fid
datay=zeros(366,39,cro);
num=0;
for i=1979:2016
    if mod(i,4)==0
        datay(1:366,i-1978,1:cro)=v10(num+1:num+366,1:cro);
        num=num+366;
    else
        datay(1:59,i-1978,1:cro)=v10(num+1:num+59,1:cro);
        datay(60,i-1978,1:cro)=nan;
        datay(61:366,i-1978,1:cro)=v10(num+60:num+365,1:cro);
        num=num+365;
    end
end
datay(1:59,39,1:cro)=v10(num+1:num+59,1:cro);
datay(60,39,1:cro)=nan;
datay(61:121,39,1:cro)=v10(num+60:num+120,1:cro);
datay(122:366,39,1:cro)=nan;
clear v10 num
datay=datay-repmat(nanmean(datay,2),1,39,1);
datay=reshape(datay,[366*39,cro]);
datay(isnan(datay(:,1))==1,:)=[];
v10=zeros(size(datay));
for i=1:cro
    y(1:14000)=datay(:,i);
    v10(:,i)=detrend(y);
end
clear datay i y

fid=fopen('C:\new\blocking\UB_1202_1004_2.txt','r');
ub=fscanf(fid,'%g',[38,213]);
fclose(fid);
clear fid
ub=ub';
ub=[zeros(153,38);ub];
ub=reshape(ub,[366*38,1]);
ub=[zeros(120,1);ub];
ub(isnan(ub))=[];
[Y,I]=sort(ub,'descend');
event=50;

inde1=[u10(:,1),v10(:,1),v10(:,4)];
inde=zeros(size(inde1));
for i=1:3
    yy(1:14000)=inde1(1:14000,i);
    inde(1:14000,i)=smooth(yy,5);
end
clear yy inde1
data=zeros(5,3,event);
for i=1:event
    k=I(i);
    data(1:5,1:3,i)=inde(k-2:k+2,1:3);
end
datas(1:3,1:50)=nanmean(data);

[Y1,I1]=sort(datas(1,:),'descend');
[Y2,I2]=sort(datas(2,:),'descend');
% [Y3,I3]=sort(datas(3,:),'descend');
for i=1:25
    num1=(i-1)*5+1;
    num2=i*5;
    k=I(I1(i));
    outdata1(num1:num2)=k-2:k+2;
    k=I(I2(i));
    outdata3(num1:num2)=k-2:k+2;
end
for i=26:50
    num1=(i-26)*5+1;
    num2=(i-25)*5;
    k=I(I1(i));
    outdata2(num1:num2)=k-2:k+2;
    k=I(I2(i));
    outdata4(num1:num2)=k-2:k+2;
end

clear u10 v10
%%

filename1='C:\new\data\ERAdata\ERA_interim_1.0_1979_2017430_sur_thermal_down_2times.nc';
% ncdisp(filename1,'/','full');
lon=ncread(filename1,'longitude');
lat=ncread(filename1,'latitude');
lon=[lon;lon(1)];
lat=lat(1:61);
data=zeros(360,61,14000);
for i=1:14000
    num=(i-1)*2+1;
    data1=ncread(filename1,'strd',[1,1,num],[360,61,2]);
    data(:,:,i)=mean(data1,3);
    clear data1;
end
clear filename1
data=data/(3*3600);
dataa=zeros(360,61,366,39);
num=0;
for i=1979:2016
    if mod(i,4)==0
        dataa(1:360,1:61,1:366,i-1978)=data(1:360,1:61,num+1:num+366);
        num=num+366;
    else
        dataa(1:360,1:61,1:59,i-1978)=data(1:360,1:61,num+1:num+59);
        dataa(1:360,1:61,60,i-1978)=nan;
        dataa(1:360,1:61,61:366,i-1978)=data(1:360,1:61,num+60:num+365);
        num=num+365;
    end
end
dataa(1:360,1:61,1:59,39)=data(1:360,1:61,num+1:num+59);
dataa(1:360,1:61,60,39)=nan;
dataa(1:360,1:61,61:121,39)=data(1:360,1:61,num+60:num+120);
dataa(1:360,1:61,122:366,39)=nan;
clear data
data1=dataa-repmat(nanmean(dataa,4),1,1,1,39);
data1=reshape(data1,[360,61,366*39]);
data1(:,:,isnan(data1(1,1,:))==1)=[];
datat=zeros(360,61,14000);
for i=1:360
    for j=1:61
        y(1:14000)=data1(i,j,:);
        datat(i,j,:)=detrend(y);
    end
end
clear i j y num
clear dataa data1

p=zeros(360,61,4);
f=zeros(360,61,4);
[~,f(:,:,1)]=composed_student( datat, datat(:,:,outdata1) );
[~,f(:,:,2)]=composed_student( datat, datat(:,:,outdata2) );
[~,f(:,:,3)]=composed_student( datat, datat(:,:,outdata3) );
[~,f(:,:,4)]=composed_student( datat, datat(:,:,outdata4) );
p(:,:,1)=nanmean(datat(:,:,outdata1),3);
p(:,:,2)=nanmean(datat(:,:,outdata2),3);
p(:,:,3)=nanmean(datat(:,:,outdata3),3);
p(:,:,4)=nanmean(datat(:,:,outdata4),3);

for i=1:360
    for j=1:61
        for k=1:4
            if (f(i,j,k)==0)
                p(i,j,k)=0;
            end
        end
    end
end
clear i j k

[lon1,lat1]=meshgrid(lon,lat);
color=[0.142,0.0,0.85;0.097,0.112,0.97;0.16,0.342,1.0;0.24,0.531,1.0;0.34,0.692,1.0;...
0.46,0.829,1.0;0.6,0.92,1.0;0.74,0.978,1.0;1.0,1.0,1.0;1.0,1.0,1.0;1.0,0.948,0.74;1.0,0.84,0.6;...
1.0,0.676,0.46;1.0,0.472,0.34;1.0,0.24,0.24;0.97,0.155,0.21;0.85,0.085,0.187;0.65,0.0,0.13];
c = interp1(linspace(0,1,size(color,1)),color,linspace(0,1,40),'pchip');
v1=(-60:2:60);
a=40;
fig11 = figure(2);
set(fig11,'Unit','centimeters')  
set(fig11,'Position',[0,0,50,26])  
axes1 = axes('Parent',fig11,'Unit','centimeters','Position',[2 15 10 10]);
datati=[p(:,:,1);p(1,:,1)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

axes2 = axes('Parent',fig11,'Unit','centimeters','Position',[14 15 10 10]);
datati=[p(:,:,2);p(1,:,2)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

axes3 = axes('Parent',fig11,'Unit','centimeters','Position',[26 15 10 10]);
datati=[p(:,:,3);p(1,:,3)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

axes4 = axes('Parent',fig11,'Unit','centimeters','Position',[38 15 10 10]);
datati=[p(:,:,4);p(1,:,4)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

clear datat datai p f i j

%%
filename='C:\new\data\ERAdata\vapor\197901.nc';
% ncdisp(filename,'/','full');
lon=ncread(filename,'longitude');
lat=ncread(filename,'latitude');
lon=[lon;lon(1)];
lat=lat(1:61);
data=ncread(filename,'tcwv');
data=data(1:360,1:61,:);
for j=2:9
    strm=num2str(j);
    filename=['C:\new\data\ERAdata\vapor\19790',strm,'.nc'];
    data1=ncread(filename,'tcwv');
    data1=data1(1:360,1:61,:);
    data=cat(3,data,data1);
end
for j=10:12
    strm=num2str(j);
    filename=['C:\new\data\ERAdata\vapor\1979',strm,'.nc'];
    data1=ncread(filename,'tcwv');
    data1=data1(1:360,1:61,:);
    data=cat(3,data,data1);
end
for i=1:37
    stry=num2str(i+1979);
    for j=1:9
        strm=num2str(j);
        filename=['C:\new\data\ERAdata\vapor\',stry,'0',strm,'.nc'];
        data1=ncread(filename,'tcwv');
        data1=data1(1:360,1:61,:);
        data=cat(3,data,data1);
    end
    for j=10:12
        strm=num2str(j);
        filename=['C:\new\data\ERAdata\vapor\',stry,strm,'.nc'];
        data1=ncread(filename,'tcwv');
        data1=data1(1:360,1:61,:);
        data=cat(3,data,data1);
    end
end
for j=1:4
    strm=num2str(j);
    filename=['C:\new\data\ERAdata\vapor\20170',strm,'.nc'];
    data1=ncread(filename,'tcwv');
    data1=data1(1:360,1:61,:);
    data=cat(3,data,data1);
end
clear filename data1 stry strm i j
dataa=zeros(360,61,366,39);
num=0;
for i=1979:2016
    if mod(i,4)==0
        dataa(1:360,1:61,1:366,i-1978)=data(1:360,1:61,num+1:num+366);
        num=num+366;
    else
        dataa(1:360,1:61,1:59,i-1978)=data(1:360,1:61,num+1:num+59);
        dataa(1:360,1:61,60,i-1978)=nan;
        dataa(1:360,1:61,61:366,i-1978)=data(1:360,1:61,num+60:num+365);
        num=num+365;
    end
end
dataa(1:360,1:61,1:59,39)=data(1:360,1:61,num+1:num+59);
dataa(1:360,1:61,60,39)=nan;
dataa(1:360,1:61,61:121,39)=data(1:360,1:61,num+60:num+120);
dataa(1:360,1:61,122:366,39)=nan;
clear data
data1=dataa-repmat(nanmean(dataa,4),1,1,1,39);
data1=reshape(data1,[360,61,366*39]);
data1(:,:,isnan(data1(1,1,:))==1)=[];
datat=zeros(360,61,14000);
for i=1:360
    for j=1:61
        y(1:14000)=data1(i,j,:);
        datat(i,j,:)=detrend(y);
    end
end
clear i j y num
clear dataa data1

p=zeros(360,61,4);
f=zeros(360,61,4);
[~,f(:,:,1)]=composed_student( datat, datat(:,:,outdata1) );
[~,f(:,:,2)]=composed_student( datat, datat(:,:,outdata2) );
[~,f(:,:,3)]=composed_student( datat, datat(:,:,outdata3) );
[~,f(:,:,4)]=composed_student( datat, datat(:,:,outdata4) );
p(:,:,1)=nanmean(datat(:,:,outdata1),3);
p(:,:,2)=nanmean(datat(:,:,outdata2),3);
p(:,:,3)=nanmean(datat(:,:,outdata3),3);
p(:,:,4)=nanmean(datat(:,:,outdata4),3);

for i=1:360
    for j=1:61
        for k=1:4
            if (f(i,j,k)==0)
                p(i,j,k)=0;
            end
        end
    end
end
clear i j k

[lon1,lat1]=meshgrid(lon,lat);
color=[0.142,0.0,0.85;0.097,0.112,0.97;0.16,0.342,1.0;0.24,0.531,1.0;0.34,0.692,1.0;...
0.46,0.829,1.0;0.6,0.92,1.0;0.74,0.978,1.0;1.0,1.0,1.0;1.0,1.0,1.0;1.0,0.948,0.74;1.0,0.84,0.6;...
1.0,0.676,0.46;1.0,0.472,0.34;1.0,0.24,0.24;0.97,0.155,0.21;0.85,0.085,0.187;0.65,0.0,0.13];
c = interp1(linspace(0,1,size(color,1)),color,linspace(0,1,40),'pchip');
v1=(-4:0.2:4);
a=4;
axes5 = axes('Parent',fig11,'Unit','centimeters','Position',[2 2 10 10]);
datati=[p(:,:,1);p(1,:,1)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

axes6 = axes('Parent',fig11,'Unit','centimeters','Position',[14 2 10 10]);
datati=[p(:,:,2);p(1,:,2)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

axes7 = axes('Parent',fig11,'Unit','centimeters','Position',[26 2 10 10]);
datati=[p(:,:,3);p(1,:,3)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

axes8 = axes('Parent',fig11,'Unit','centimeters','Position',[38 2 10 10]);
datati=[p(:,:,4);p(1,:,4)]';
m_proj('stereographic','lon',60,'lat',90,'radius',40,'rotangle',0);
m_contourf(lon1,lat1,datati,v1,'linestyle','none');
m_coast('line','color','k');
m_grid('xtick',12,'ytick',5,'yticklabel',[],'tickdir','out','linest','--');
caxis([-a a]);
colormap(c)

clear datat datai p f i j
set(fig11,'PaperUnits','centimeters','PaperPosition',[0 0 50 26])
print -f2 -r600 -dpng com_ubuv11_2;  
clear fig11 axes1 axes2 axes3 axes4
