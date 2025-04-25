clc,clear,close all;
load('E:\1_matlab\ERA_TM_UB\MatFile\lat.mat');
load('E:\1_matlab\ERA_TM_UB\MatFile\lon.mat');
load('E:\1_matlab\ERA_TM_UB\MatFile\z_ndjfm_daily.mat');
z500nm = calNDJFM(z_ndjfm_daily);
%%
for s = 9:39
    s
    [TM1{1,s},eventDD11{1,s}]=TM90(z_ndjfm_daily{1,s},z500nm{1,s},lonData,latData,30,90,3);
end
num1=0;
for s = 1:39
    num1 = num1+size(eventDD11{1,s},1);
end
%% 
for s = 9:39
    eventDD11{1,s}(find(eventDD11{1,s}(:,1)<31),:)=[];
    eventDD11{1,s}(find(eventDD11{1,s}(:,1)>(size(TM1{1,s},2)-31)),:)=[];
end
num2=0;
for s = 1:39
    num2 = num2+size(eventDD11{1,s},1);
end
%%
[TMt,eventDDt,c,ub,ubb,ubbb,ub1,ube,event1]=TM90(z_djf_daily{1,2},dz{1,2},lonData,latData,30,90,3);
%% µ÷ÊÔ
close all;
lon=lonData(73:125);
for s=1:39  
    figure('position',[30 30 1800 900])
    if size(eventDD11{1,s},1) > 0;
        figure(s);
        
      for t1=1:size(eventDD11{1,s}(:,1),1) %UBÊÂ¼þ
          subplot(1,size(eventDD11{1,s}(:,1),1),t1);
m_proj('equidistant Cylindrical','lat',[0 90],'long',[0 130]);
m_coast('patch',[.8 .8 .8]);
m_grid('box','none','linestyle','none','backcolor',[.9 .99 1]);
hold on

%cc=1:0.5:25;
%[C,h]=m_contourf(lon,latData,bl_events_freq'*100,cc,'linestyle','none');
m_contour(lonData,latData,z500nm{1,s}(:,:,eventDD11{1,s}(t1,1))',[-600:40:600]);

colorbar;
colormap(jet);
caxis([-400,400]);
title([num2str(s+1978),'     500hpa z-anomaly      ']);
      end
      saveas(gcf,['e:/sea ice/figures/testChxd/',num2str(s+1978),'z_500.png']); 
    end
    

end
close all;
%%
function [TM,eventDD,c,ub,ubb,ubbb,ub1,ube,event1]=TM90(hgt,hgta,lon,lat,lonW,lonE,DD)

%__________________________________________________________________________

% input
% hgt  is total z  [360,46,timesteps] EQto90N
% hgta is total za [360,46,timesteps] EQto90N
% lon  is [1:1:360]
% lat  is [0:2:90]
% For Ural Blocking, lonW is 30, lonE is 90
% DD   is the duration threshold, commonly 3.

% output
% TM is the blocking frequency for each longitude each day. 1 = occurance.
% eventDD is the features of blocking event: 
%  eventDD.1: peakday (with the largest anomaly). The lag0 to composite.
%  eventDD.2: height anomaly (intensity, approximately)
%  eventDD.3: duration (given by this algorithms), but sometimes we better have a look at its height anomalies.

%__________________________________________________________________________


% 5-day smooth to fliter the small and synoptic eddies
for i=3:size(hgt,3)-2;
    hgt(:,:,i)= (hgt(:,:,i-2)+ hgt(:,:,i-1)+ hgt(:,:,i)+ hgt(:,:,i+1)+ hgt(:,:,i+2))./5;
end

%main body
y1=find(lat==80);y1n=find(lat==85);y1s=find(lat==75);
y2=find(lat==60);y2n=find(lat==65);y2s=find(lat==55);
y3=find(lat==40);y3n=find(lat==45);y3s=find(lat==35);
     for day=1:size(hgt,3);
         day;
         for i=1:size(hgt,1);
         GHGS(i,day,1)=(hgt(i,y2n,day)-hgt(i,y3n,day))/20;
         GHGN(i,day,1)=(hgt(i,y1n,day)-hgt(i,y2n,day))/20;
         GHGS(i,day,2)=(hgt(i,y2,day)-hgt(i,y3,day))/20;
         GHGN(i,day,2)=(hgt(i,y1,day)-hgt(i,y2,day))/20;
         GHGS(i,day,3)=(hgt(i,y2s,day)-hgt(i,y3s,day))/20;
         GHGN(i,day,3)=(hgt(i,y1s,day)-hgt(i,y2s,day))/20;      
        if ((GHGS(i,day,1)>0) && (GHGN(i,day,1)<-10))||((GHGS(i,day,2)>0) && (GHGN(i,day,2)<-10))||((GHGS(i,day,3)>0) && (GHGN(i,day,3)<-10)) ;    
               TM(i,day)=1;
               amp(i,day)=max(hgta(i,:,day));
               a=find(hgta(i,:,day)==amp(i,day));
               if length(a)==1;
                   lat(i,day)=a;
               else
                   lat(i,day,year)=sum(a)/length(a);
                    lat(i,day)=max(a);
               end
               lat(i,day)=(lat(i,day)-1)*1;
            end
        end    
         
     end
     
     for t=1:size(TM,2); 
     a=max(amp(:,t));
        b=find(amp(:,t)==a);
        if b==1;
            c(t,1)=t;c(t,2)=amp(b,t);c(t,3)=-182.5+b*2.5;
            c(t,4)=lat(b,t)*2.5;c(t,5)=1;
        else
            b=max(b);
            c(t,1)=t;c(t,2)=amp(b,t);c(t,3)=-182.5+b*2.5;
            c(t,4)=lat(b,t)*2.5;c(t,5)=1;
        end
     if a==0;
        c(t,1)=0;c(t,2)=0;c(t,5)=0; 
     end
     end;
    

% event 
ub=zeros(size(c,1),5); 
for i=1:size(c,1);
    if (c(i,3)>=lonW) && (c(i,3)<=lonE)  
        ub(i,:)=c(i,:);
    end
end
count=1;
for i=1+1:size(ub,1)-1;
    if ub(i,5)==0 && ub(i-1,5)==0 && ub(i+1,5)==0;
    else
        ubb(count,:)=ub(i,:);        count=count+1;
    end
end
%
count=1;
for i=1:size(ubb,1)-1;
    if ubb(i,1)==ubb(i+1,1)-1;
        ubbb(count,:)=ubb(i,:);        count=count+1;
    else
        ubbb(count,:)=ubb(i,:);        count=count+1;
        ubbb(count,:)=[0,0,0,0,0];     count=count+1;
    end
end
ubbb(count,:)=ubb(end,:);   
%
count=1;
for i=1+1:size(ubbb,1)-1;
    if ubbb(i,5)==0 && ubbb(i-1,5)==0 && ubbb(i+1,5)==0;
    else
        ub1(count,:)=ubbb(i,:);         count=count+1;
    end
end
count=1;
for i=1+1:size(ub1,1)-1;
    if ub1(i,5)~=0 || ub1(i+1,5)~=0 ;
                ube(count,:)=ub1(i,:);        count=count+1;
    end
end

ube(end+1,:)=0;
count=1;count1=1;
aa = find(ube(:,1));
for i=aa(1):size(ube,1);
    if ube(i,5)~=0;
        event1(count,:)=ube(i,:);
        count=count+1;
    else ube(i,5)==0;
        a=max(event1(:,2));
        b=find (event1(:,2)==a);
        peakday(count1,:)=event1(b,:); 
        peakday(count1,5)=count-1; 
        count1=count1+1;clear event1;count=1;
    end
end   
peakday
eventDD=peakday(find(peakday(:,5)>=DD),:); 
eventDD=eventDD(:,[1,2,5]);   

      
end    % function end

function [res] = calNDJFM(var)
l=0;
k=0;
yrs = size(var,2);
for s=1:yrs
    fe=size(var{1,s},3);
    if fe == 152 % notice feb 28or29
        l=l+1;
        zf(:,:,l)=var{1,s}(:,:,121); % 29th
        zq(:,:,1:120,l)=var{1,s}(:,:,1:120); % before 29th
        zh(:,:,1:31,l)=var{1,s}(:,:,122:152); % after 29th
    else
        k=k+1;
        zz(:,:,:,k)=var{1,s};
    end
end
zfm=nanmean(zf,3);
zq1=cat(4,zq,zz(:,:,1:120,:));
zqm=nanmean(zq1,4);
zh1=cat(4,zh,zz(:,:,121:151,:));
zhm=nanmean(zh1,4);
for s=1:yrs
    fe=size(var{1,s},3);
    if fe == 152
        res{1,s}(:,:,1:120)=var{1,s}(:,:,1:120)-zqm; % anomaly
        res{1,s}(:,:,121)=var{1,s}(:,:,121)-zfm;
        res{1,s}(:,:,122:152)=var{1,s}(:,:,122:152)-zhm;
    else 
        res{1,s}(:,:,1:120)=var{1,s}(:,:,1:120)-zqm;
        res{1,s}(:,:,121:151)=var{1,s}(:,:,121:151)-zhm;
    end
end
end
   
