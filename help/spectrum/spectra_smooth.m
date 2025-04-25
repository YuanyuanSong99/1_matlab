%Spactre function;
%对随即信号进行功率谱分析，用MATLAB的标准FFT格式
%[output,F,T,SS,maxl,minl]=spactre(data,dt,num,level)
% input :
% data  ------input data vector;
% dt   ------minmum data frequency;
% num   ----  smooth data number;
% level   ------ significant level;
% output :
% S------data spactrum;
% F,T   -------frequency and period serial;
% maxl,minl ------up and down significant level;


function [S,F,T,SS,maxl,minl]=spactre(data,dt,num,level)

N=max(size(data));          
dw=2*pi/(N*dt);
df=dw/(2*pi);

X=data-mean(data);   %standard
A=fft(X);              %F氏变换
for m=1:N
    F(m)=df*m;
    T(m)=1./F(m);
    S(m)=dt*(abs(A(m)))^2/(2*pi*N);
end

%--------对粗谱进行平滑 做插值
 %numn=2;   %对多少个数进行滑动平均
 if ( num==0 || isempty(num) )
     SS=S;
 else 
     %----method 1--
     SS=smooth(S,num);
   %{
  %---method  2--  
   mm=N/(num);    %得到平均值的个数 
   kk=0; 
   sum=0;
  for i=1:mm
      kk=kk+1;
      for j=1:num
          sum=sum+S((i-1)*num+j);
      end
      y(kk)=sum/num;
      atic(kk)=num*kk-num/2;
      sum=0;
  end  
 %下面做插值拟合
 q1=2*atic(length(atic))-atic(length(atic)-1);
 qq1=2*atic(1)-atic(2);
 data1=[y(2) y y(max(size(y))-1)];da1x=[qq1 atic q1];
 %周期条件
 pp1=csape(da1x,data1,[0 0]); 
 SS=fnval(pp1,1:N);  
 %}
 end
 
%-------significant level and freedom------
if isempty(level)
  maxl(size(SS))=0;
  minl(size(SS))=0;
else
  u=2*num;          %自由度80；
  B=level;          %置信区间；
  a=chi2inv(1-B,u);
  b=chi2inv(B,u);
  maxl=u*SS/a;
  minl=u*SS/b;
end

