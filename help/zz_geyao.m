%阻塞挑选 渣渣熊
clear all
%load f:\pna\hgt\data\hdata cc hh ddd
load  f:\pvy1\era\data\hdata 
hh=hh/9.8;
hhh(1:72,:,:)=hh(73:144,:,:);
hhh(73:144,:,:)=hh(1:72,:,:);
[nx ny nt]=size(hh);
for i=1:nx
    for j=1:ny
            hgtave(i,j,:)=smooth(hhh(i,j,:));  
        end
end
lat=90:-2.5:0;
lon=0:2.5:357.5;
zuselat=[3 5 7 11 13 15 19 21 23];
%zuselat=[9 11 13 15 17 19 21 23 25];
%zuselat=[3 5 7 9 11 13 15 17 19];%要用到的纬度
%经度选择150-240 180-240 150-260
for a=1:9
    h(:,a,:)=hgtave(17:33,zuselat(a),:);
end    
[nx ny nt]=size(h);
%TM指数判断
p=zeros(nx,nt);
for t=1:nt
    for l=1:nx
g1=h(l,1,t)-h(l,4,t);
g2=h(l,2,t)-h(l,5,t);
g3=h(l,3,t)-h(l,6,t);
if g1<-200&h(l,4,t)>h(l,7,t)%△为5
    p(l,t)=1;
elseif g2<-200&h(l,5,t)>h(l,8,t)%△为0
    p(l,t)=1;
elseif g3<-200&h(l,6,t)>h(l,9,t)%△为-5
    p(l,t)=1;
end
    end
end
%%%%%
p_2=sum(p,1);
p_3=zeros(1,nt);
for i=1:nt
if p_2(i)>=5
    p_3(i)=1;
end
end
pp=zeros(1,nt);
pp(1)=0;pp(2)=0;pp(nt)=0;pp(nt-1)=0;
for k=3:nt-2
    if p_3(k)==1
        if p_3(k-1)==0&p_3(k+1)==1&p_3(k+2)==1
            pp(k)=1;
        elseif p_3(k-1)==1&p_3(k+1)==1
            pp(k)=1;
        elseif p_3(k+1)==0&p_3(k-1)==1&p_3(k-2)==1
            pp(k)=1;
        end
    end
end    
%****************************
k1=0;k2=0;
for i=2:nt
    if pp(i-1)==0&pp(i)==1&pp(i+1)==1
       % pp(i)=2;
        k1=k1+1;
        kaishi(k1)=i;
    elseif pp(i-1)==1&pp(i)==1&pp(i+1)==0
        %pp(i)=3;
        k2=k2+1;
        jieshu(k2)=i;
    end
end
zhouqi=jieshu-kaishi+1;
%判断一下lag0，ghgs最大的那一天
ghgs=zeros(3,nx,nt);
for t=1:nt
    if pp(t)==1;
      for r=1:nx
        ghgs(1,r,t)=h(r,4,t)-h(r,7,t);
        ghgs(2,r,t)=h(r,5,t)-h(r,8,t);
        ghgs(3,r,t)=h(r,6,t)-h(r,9,t);
      end
    end
end
Z=squeeze(max(max(ghgs)));
a=size(kaishi);
for i=1:a(2)
    [m1,m2]=max(Z(kaishi(i):jieshu(i)));
    lag0(i)=kaishi(i)+m2-1;
end
zuse=ddd(lag0,:);
n=1;
a=size(lag0)
for i=1:(a(2)-1)
    c=lag0(i+1)-lag0(i);
    if c<=10
        dz(n)=i;
        dz(n+1)=i+1;
        zz(n,:)=zuse(i,:);
        zz(n+1,:)=zuse(i+1,:);
        n=n+2;
    end
end
dz=unique(dz);
a=size(dz);
for i=-a(2):-1
    zuse(dz(-i),:)=[];
end
a=size(zuse);
for i=-a(1):-1
    if zuse(-i,2)==3|zuse(-i,2)==11
        zuse(-i,:)=[];
    end
end
save f:\pvy1\mat\zuse zuse


        










    
