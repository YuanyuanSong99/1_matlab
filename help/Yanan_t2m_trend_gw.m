clc;clear;clear all memory;
%%%%%%%处理气温数据
t2m_location='C:\new\data\ERAdata\Monthly\Surface2_2\';
filename=[t2m_location,'1979.nc'];
% ncdisp(filename,'/','full');
lon=ncread(filename,'longitude');
lat=ncread(filename,'latitude');
tim=ncread(filename,'time');
a=length(lon);b=(length(lat)+1)/2;c=length(tim);
clear lat lon tim
data=ncread(filename,'t2m',[1,1,1],[a,b,c]);
clear filename
for i=1980:2017
    filename=[t2m_location,num2str(i),'.nc'];
    data1=ncread(filename,'t2m',[1,1,1],[a,b,c]);
    data=cat(3,data,data1);
    clear filename data1
end
clear i
filename=[t2m_location,'2018.nc'];
data1=ncread(filename,'t2m',[1,1,1],[a,b,4]);
data=cat(3,data,data1);
clear filename data1
clear lon lat tim a b c
%%%%%%%%%%
t2m_location='C:\new\data\ERA20C\monthly\surface2_2\';
filename=[t2m_location,'1900.nc'];
% ncdisp(filename,'/','full');
lon=ncread(filename,'longitude');
lat=ncread(filename,'latitude');
tim=ncread(filename,'time');
a=length(lon);b=(length(lat)+1)/2;c=length(tim);
clear lat lon tim
datac=ncread(filename,'t2m',[1,1,1],[a,b,c]);
clear filename
for i=1901:2010
    filename=[t2m_location,num2str(i),'.nc'];
    data1=ncread(filename,'t2m',[1,1,1],[a,b,c]);
    datac=cat(3,datac,data1);
    clear filename data1
end
clear i
clear filename data1
clear lon lat tim a b c
clear t2m_location
%%%%%%%%%%%%%%%%%%%%%%%%
a=180;b=46;c=118;
datac_data=zeros(size(datac));
for i=1:a
    for j=1:b
        x(1:384)=datac(i,j,79*12+1:111*12);
        y(1:384)=data(i,j,1:32*12);
        b1=polyfit(x,y,1);
        datac_data(i,j,1:111*12)=datac(i,j,1:111*12)*b1(1)+b1(2);
        clear b1 x y
    end
end  
data_new=cat(3,datac_data(:,:,3:79*12),data(:,:,1:39*12+2));
clear i j
clear datac data datac_data
%%%%%%%%%%%%%%%%%%%%
datat=reshape(data_new,[a,b,12,c]);
clear data_new
sat=zeros(a,b,4,c);
mam=[1,2,3];jja=[4,5,6];son=[7,8,9];djf=[10,11,12];
sat(1:a,1:b,1,1:c)=nanmean(datat(:,:,mam,:),3);
sat(1:a,1:b,2,1:c)=nanmean(datat(:,:,jja,:),3);
sat(1:a,1:b,3,1:c)=nanmean(datat(:,:,son,:),3);
sat(1:a,1:b,4,1:c)=nanmean(datat(:,:,djf,:),3);
clear mam jja son djf
clear datat
%%%%%%%%%%%%%%%%%%%%%
fid=fopen('C:\new\4EAcooling\0809_cmip\globalmean_CMIP5_T18502100.txt','r');
gw2=fscanf(fid,'%g');
fclose(fid);
clear fid ans
gw=gw2(51:168);%1900-2017
sat_gw=zeros(size(sat));
for i=1:a
    for j=1:b
        for k=1:4
            x(1:c)=gw;
            y(1:c)=sat(i,j,k,1:c);
            b1=polyfit(x,y,1);
            sat_gw(i,j,k,1:c)=y-x*b1(1)-b1(2);
            clear x y b1  
        end
    end
end
clear i j k
clear gw gw2 sat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b,c,d]=size(sat_gw);

t1=1979:2017;
p1=zeros(a,b,c,2);
y=zeros(39,1);
for i=1:a
    for j=1:b
        for k=1:c
            y(:)=sat_gw(i,j,k,(1980-1900):(2018-1900));
            p1(i,j,k,:)=polyfit(t1,y',1);
        end
    end
end
clear i j k y

t2=2000:2013;
p2=zeros(a,b,c,2);
y=zeros(14,1);
for i=1:a
    for j=1:b
        for k=1:c
            y(:)=sat_gw(i,j,k,(2001-1900):(2014-1900));
            p2(i,j,k,:)=polyfit(t2,y',1);
        end
    end
end
clear i j k y

lon=(-90:2:270);
lat=(90:-2:0);
[lon1,lat1]=meshgrid(lon,lat);
rangelon=[-90,270];
rangelat=[0,90];
v1=(-10:0.2:10);
load C:\new\4EAcooling\0805_ttrend\t2m_trend.mat
cc = interp1(linspace(0,1,size(mycmap,1)),mycmap,linspace(0,1,80),'pchip');
clear mycmap

strname={'MAM','JJA','SON','DJF'};
for i=1:4
    fig11 = figure(1);
    set(fig11,'Unit','centimeters')  
    set(fig11,'Position',[0,0,33,18])  
    axes1 = axes('Parent',fig11,'Unit','centimeters','Position',[2 1 30 16]);
    datat=[p1(136:180,:,i,1);p1(1:136,:,i,1)]';
    m_proj('Equidistant','lon',rangelon,'lat',rangelat);
    m_contourf(lon1,lat1,datat*10,v1,'linestyle','none');
    m_coast('line','color',[0.5 0.5 0.5],'LineWidth',1.5);
    m_grid('xtick',13,'ytick',10,'tickdir','out','fontsize',16);
    caxis([-4 4])
    colormap(cc)
    colorbar('SouthOutside','position',[2/33 1.5/18 30/33 0.5/18],'fontsize',16)
    text(-1.8,2.4,strcat(strname(i),' surface temperature trends (1979–2017)'),'fontsize',24)
    research_rangee([40,60],[60,120],'-k',1.5)
    clear axes1
    set(fig11,'PaperUnits','centimeters','PaperPosition',[0 0 33 18])
    print(fig11,'-dpng','-r300',['C:\new\4EAcooling\0805_ttrend\trendtu\tad_trend19792017_',num2str(i),'.png']);
    clear fig11
    close(1)
    
    fig11 = figure(1);
    set(fig11,'Unit','centimeters')  
    set(fig11,'Position',[0,0,33,18])  
    axes1 = axes('Parent',fig11,'Unit','centimeters','Position',[2 1 30 16]);
    datat=[p2(136:180,:,i,1);p2(1:136,:,i,1)]';
    m_proj('Equidistant','lon',rangelon,'lat',rangelat);
    m_contourf(lon1,lat1,datat*10,v1,'linestyle','none');
    m_coast('line','color',[0.5 0.5 0.5],'LineWidth',1.5);
    m_grid('xtick',13,'ytick',10,'tickdir','out','fontsize',16);
    caxis([-4 4])
    colormap(cc)
    colorbar('SouthOutside','position',[2/33 1.5/18 30/33 0.5/18],'fontsize',16)
    text(-1.8,2.4,strcat(strname(i),' surface temperature trends (2000–2013)'),'fontsize',24)
    research_rangee([40,60],[60,120],'-k',1.5)
    clear axes1
    set(fig11,'PaperUnits','centimeters','PaperPosition',[0 0 33 18])
    print(fig11,'-dpng','-r300',['C:\new\4EAcooling\0805_ttrend\trendtu\tad_trend20002013_',num2str(i),'.png']);
    clear fig11
    close(1)
end


%%
%%%%%%%%%处理温度场
clc;clear;clear all memory;
filename='C:\new\data\Ncep\monthly\surface\air.mon.mean.nc';
% ncdisp(filename,'/','full');
lon1=ncread(filename,'lon');
lat1=ncread(filename,'lat');
data=ncread(filename,'air');
data=data(:,1:37,:);
clear filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename='C:\new\data\Ncep\monthly\20C\air.sig995.mon.mean.nc';
% ncdisp(filename,'/','full');
lon2=ncread(filename,'lon');
lat2=ncread(filename,'lat');
datac=ncread(filename,'air');
datac=datac-273.15;
clear filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%
lon3=double(lon2);
lat3=double(lat2);
[lon3,lat3]=meshgrid(lon3,lat3);
lon4=double(lon1);
lat4=double(lat1);
[lon4,lat4]=meshgrid(lon4,lat4);

datac2=double(datac);
for i=1:length(datac)
    datac3(:,:,i)=griddata(lon3',lat3',datac2(:,:,i),lon4',lat4');
end
clear i lat1 lat2 lat3 lat4 lon1 lon2 lon3 lon4
clear datac2 datac
datac=datac3(:,1:37,:);
clear datac3
%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=144;b=37;
datac_data=zeros(size(datac));
for i=1:a
    for j=1:b
        x(1:(2014-1947)*12)=datac(i,j,(1948-1851)*12+1:(2014-1850)*12);
        y(1:(2014-1947)*12)=data(i,j,1:(2014-1947)*12);
        b1=polyfit(x,y,1);
        datac_data(i,j,:)=datac(i,j,:)*b1(1)+b1(2);
        clear b1 x y
    end
end  
data_new=cat(3,datac_data(:,:,(1900-1851)*12+3:(1948-1851)*12),data(:,:,1:(2018-1948)*12+2));
clear i j
clear datac data datac_data
%%%%%%%%%%%%%%%%%%%%
c=118;
datat=reshape(data_new,[a,b,12,c]);
clear data_new
%%%%%%%%%%%%%%%%%%%%%
sat=zeros(a,b,4,c);
mam=[1,2,3];jja=[4,5,6];son=[7,8,9];djf=[10,11,12];
sat(1:a,1:b,1,1:c)=nanmean(datat(:,:,mam,:),3);
sat(1:a,1:b,2,1:c)=nanmean(datat(:,:,jja,:),3);
sat(1:a,1:b,3,1:c)=nanmean(datat(:,:,son,:),3);
sat(1:a,1:b,4,1:c)=nanmean(datat(:,:,djf,:),3);
clear mam jja son djf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('C:\new\4EAcooling\0809_cmip\globalmean_CMIP5_T18502100.txt','r');
gw2=fscanf(fid,'%g');
fclose(fid);
clear fid
x=gw2(51:168)';%1900-2017
sat_gw=zeros(a,b,4,c);
for k=1:4
    for i=1:a
        for j=1:b     
            y(1:c)=sat(i,j,k,1:c);
            b1=polyfit(x,y,1);
            sat_gw(i,j,k,1:c)=y-x*b1(1)-b1(2);
            clear y b1  
        end
    end
end
clear i j k
clear x gw2 sat datat
clear a b c
[a,b,c,d]=size(sat_gw);

t1=1993:2013;
p1=zeros(a,b,c,2);
y=zeros(21,1);
for i=1:a
    for j=1:b
        for k=1:c
            y(:)=sat_gw(i,j,k,(1994-1900):(2014-1900));
            p1(i,j,k,:)=polyfit(t1,y',1);
        end
    end
end
clear i j k y

t2=2000:2013;
p2=zeros(a,b,c,2);
y=zeros(14,1);
for i=1:a
    for j=1:b
        for k=1:c
            y(:)=sat_gw(i,j,k,(2001-1900):(2014-1900));
            p2(i,j,k,:)=polyfit(t2,y',1);
        end
    end
end
clear i j k y

lon=(-90:2.5:270);
lat=(90:-2.5:0);
[lon1,lat1]=meshgrid(lon,lat);
rangelon=[-90,270];
rangelat=[0,90];
v1=(-10:0.2:10);
load C:\new\4EAcooling\0805_ttrend\t2m_trend.mat
cc = interp1(linspace(0,1,size(mycmap,1)),mycmap,linspace(0,1,80),'pchip');
clear mycmap

strname={'MAM','JJA','SON','DJF'};
for i=1:4
    fig11 = figure(1);
    set(fig11,'Unit','centimeters')  
    set(fig11,'Position',[0,0,33,18])  
    axes1 = axes('Parent',fig11,'Unit','centimeters','Position',[2 1 30 16]);
    datat=[p1(109:144,:,i,1);p1(1:109,:,i,1)]';
    m_proj('Equidistant','lon',rangelon,'lat',rangelat);
    m_contourf(lon1,lat1,datat*10,v1,'linestyle','none');
    m_coast('line','color',[0.5 0.5 0.5],'LineWidth',1.5);
    m_grid('xtick',13,'ytick',10,'tickdir','out','fontsize',16);
    caxis([-4 4])
    colormap(cc)
    colorbar('SouthOutside','position',[2/33 1.5/18 30/33 0.5/18],'fontsize',16)
    text(-1.8,2.4,strcat(strname(i),' surface temperature trends (1993–2013)'),'fontsize',24)
%     research_rangee([40,60],[60,120],'-k',1.5)
    clear axes1
    set(fig11,'PaperUnits','centimeters','PaperPosition',[0 0 33 18])
    print(fig11,'-dpng','-r300',['C:\new\4EAcooling\0805_ttrend\trendtu\tad_ncep_trend19932013_',num2str(i),'.png']);
    clear fig11
    close(1)
    
    fig11 = figure(1);
    set(fig11,'Unit','centimeters')  
    set(fig11,'Position',[0,0,33,18])  
    axes1 = axes('Parent',fig11,'Unit','centimeters','Position',[2 1 30 16]);
    datat=[p2(109:144,:,i,1);p2(1:109,:,i,1)]';
    m_proj('Equidistant','lon',rangelon,'lat',rangelat);
    m_contourf(lon1,lat1,datat*10,v1,'linestyle','none');
    m_coast('line','color',[0.5 0.5 0.5],'LineWidth',1.5);
    m_grid('xtick',13,'ytick',10,'tickdir','out','fontsize',16);
    caxis([-4 4])
    colormap(cc)
    colorbar('SouthOutside','position',[2/33 1.5/18 30/33 0.5/18],'fontsize',16)
    text(-1.8,2.4,strcat(strname(i),' surface temperature trends (2000–2013)'),'fontsize',24)
%     research_rangee([40,60],[60,120],'-k',1.5)
    clear axes1
    set(fig11,'PaperUnits','centimeters','PaperPosition',[0 0 33 18])
    print(fig11,'-dpng','-r300',['C:\new\4EAcooling\0805_ttrend\trendtu\tad_ncep_trend20002013_',num2str(i),'.png']);
    clear fig11
    close(1)
end
