sst=importdata('e:/1硕士一年级/气候统计/homework/EEMD/shuju.txt');
%for k=1:100
    
    rslt(:,:)=eemd(sst(:,2),0.2,10000);
%end

plot(rslt(:,2));
hold on
plot(rslt(:,3)-1);
plot(rslt(:,4)-2);
plot(rslt(:,5)-3);
plot(rslt(:,6)-4);
plot(rslt(:,7)-5);
title('EEMD');
saveas(gcf,'e:/1硕士一年级/气候统计/homework/EEMD/EEMD.png')
close all;

c1=rslt(:,4);
kk=0;
for i=1:(length(c1)-1)
    if c1(i)*c1(i+1)<0
        kk=kk+1;
    end
end
T1=length(c1)*2/kk;

component=rslt(:,2:7);
contribution=std(component).^2/sum(std(component).^2)*100;

sigline99=significanceIMF(rslt);
lnt=sort(sigline99(:,1));
lne=sort(sigline99(:,2));
l=fliplr(lne');
plot(lnt,l);


