%function [ssp,ss_c,t,ttt,kkk,r,freq]=power_spec(yi,m,a,tag,pro)    
% copyright by Jasion Hu,Qingdao,China. 2004.5.1
% ----input:
% ---yi:  is the time serial, a vector.
% ---m :  is the length of correlation. m 相关的长度
% name is the name of the time series, used for plotting the figure; 
% ---if tag=0, it will  use the function semilogy, else use plot, %otherwise, no figures are made.
%     and the the lowest frequency(the first one)has been ignored; 
% ---pro is the degree of confidence, the default value is 0.95;pro是置信度，默认值是0.95

% ----output:
% ssp is the spectrum power of the input data;   ssp输入数据的频谱功率
% ss_c is the spectrum power of standard red noise; ss_c频谱功率的红噪声标准
%--- ttt:  contains the remarkable periods with red noise test;ttt包含了红噪声检验的显着周期
%--- kkk:  contains the index of the  remarkable periods;  kkk包含显著周期的指数
%--- r is the lag correlation,if the r(2,1) is lower than zero or is very close to zero, 
%     you maybe use white noise test,otherwise you may use the red noise;r是滞后相关，如果r(2,1)小于0，或是接近于0，这时要用白噪声检验
%--- freq:  frquency or the period
% from the figure, where the value of real power is larger than the standard one,
% the frquency or the period is remarkable.
%--- ss_c : 谱检验的标准置信区间
%--- ssg_r: 红谱检验
%--- path:  the path of saving figure

function [ssp,t,ttt,kkk,r,freq,ss_c]=power_spectrum(yi,m,tag,name)    
%function [ssp,t,ttt,kkk,r,freq,ssg_r]=power_spectrum(yi,m,tag,name)    
[dm dn]=size(yi);
if dm>dn
   n=size(yi,1);
elseif dn>dm
   n=dn;yi=yi';
end
s=std(yi);
yy=mean(yi);
r=zeros(m+1,1);
for j=0:m
    for i=1:n-j
        r(j+1,1)=r(j+1,1)+((yi(i,1)-yy)*(yi(i+j,1)-yy)/(s*s))/(n-j);
    end
end

% --- chu pu -----
for k=0:m
    t(k+1,1)=2*m/(k+eps);
    freq(k+1,1)=1/t(k+1,1);
end

% --- smooth ----
ss(1,1)=sum(r(2:m,1))/m+(r(1,1)+r(m+1,1))/(2*m);
for i=2:m
    ssm=0;
for j=1:m-1
    ssm=ssm+r(j+1,1)*cos((i-1)*pi*j/m);
end
     ss(i,1)=(r(1,1)+2*ssm+r(m+1,1)*cos((i-1)*pi))/m;
end

ssm=0;
for j=1:m-1
    ssm=ssm+r(j+1,1)*(-1)^j/m;
end
ss(m+1,1)=ssm+(r(1,1)+(-1)^m*r(m+1,1))/(2*m);

% --use the window smooth
ssp(1,1)=(ss(1,1)+ss(2,1))/2;
ssp(m+1,1)=(ss(m,1)+ss(m+1,1))/2;
for k=1:m-1
    ssp(k+1,1)=sum(ss(k:k+2,1))/4+ss(k+1,1)/4;
end

% ---- verify 
% white noise test  
% red noise test
    sss=(ssp(1,1)+ssp(m+1,1))/(2*m)+sum(ssp(2:k),1)/m;
    for k=0:m
        ssg_r(k+1,1)=sss*(1-r(2,1)^2)/(1+r(2,1)^2-2*r(2,1)*cos(pi*k/m));
    end;
    
% ---- the free degree（自由度）
   v=(2*n-m/2)/m;
   p=chi2inv(0.95,v);
   
   ss_c=ssg_r*p/v;   %-----谱检验的标准置信区间
   num_t=0;
   for hh=2:size(ss_c,1)
       if abs(ssp(hh,1))>=ss_c(hh,1)
           num_t=num_t+1;
           ttt(num_t,1)=t(hh,1);
           kkk(num_t,1)=hh-1;
       else
           ttt=0;
           kkk=0;
       end;
   end
 k=0:m;
 k=k';
%{
%fig1=strcat(path,'powerspeT.emf');%fig2=strcat(path,'powerspeK.emf');
fig1=strcat('powerspeTK_',name,'.emf');
ax=2;
figure(11);
if tag==1
   % h=plot(freq(2:end,1),abs(ssp(2:end,1)),'m',freq(2:end,1),ss_c(2:m+1,1),'r:');
   % xlabel('Frequency');ylabel('Power');
   %subplot(211)
   h=plot(1./freq(2:end,1),abs(ssp(2:end,1)),'b','Linewidth',1.5);hold on;
   plot(1./freq(2:end,1),ss_c(2:m+1,1),'r:','Linewidth',1.5);
   plot(1./freq(2:end,1),ssg_r(2:m+1,1),'m:','Linewidth',1.5);hold off;
   xlim([2 fix(length(ssp(:,1))/ax)]);
   xlabel('Periods/(month)');ylabel('Power');
   title(['The estimate of power spectra']);
   legend('Real','Standard','red noise');
   set(gcf,'color','w');
   %{
   subplot(212); 
   plot(k(2:end,1),abs(ssp(2:end,1)),'r','Linewidth',1.5); hold on;
   plot(k(2:end,1),ss_c(2:m+1,1),'b:','Linewidth',1.5); hold off;
   xlim([2 fix(length(ssp(:,1))/ax)]);
   xlabel('K');ylabel('Power');
   title(['The estimate of power spectra']);
   legend('Real','Standard');
   set(gcf,'color','w');
  % saveas(gcf,fig1);
   %}
else
%    h=semilogy(freq(2:end,1),abs(ssp(2:end,1)),'m',freq(2:end,1),ss_c(2:m+1,1),'r:');
%    xlabel('Frequency');ylabel('Power');
   h=semilogy(1./freq(2:end,1),abs(ssp(2:end,1)),'m',1./freq(2:end,1),ss_c(2:m+1,1),'r:');
   xlim([2 fix(length(ssp(:,1))/2)]);
   xlabel('Periods');ylabel('Power');
   title(['The estimate of power spectra']);
   legend('Real','Standard');
   set(gcf,'color','w');
   saveas(gcf,fig1);
   subplot(212) %figure(2);
   semilogy(k(2:end,1),abs(ssp(2:end,1)),'m',k(2:end,1),ss_c(2:m+1,1),'r:');
   xlabel('K');ylabel('Power');
   title(['The estimate of power spectra']);
   legend('Real','Standard'); 
   set(gcf,'color','w');
 %  saveas(gcf,fig2);
end
%}