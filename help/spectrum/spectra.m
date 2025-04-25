%Spactre function;
%对随即信号进行功率谱分析，用MATLAB的标准FFT格式
%[output,F,T]=spactre(data,dt)
% input :
% data  ------input data vector;
% dt   ------minmum data frequency;
% output :
% S------data spactrum;功率谱，功率谱密度，单位：mm^2*month (h^2*t) 
% F,T   -------frequency and period serial;
% 没有作置信检验

function [S,F,T]=spectra(data,dt)

N=max(size(data));  
dt=1/12;       %将月均的周期单位换为年周期单位
dw=2*pi/(N*dt);
df=dw/(2*pi);

X=data-mean(data);   %standard
A=fft(X);              %F氏变换
for m=1:N
    F(m)=df*m;
    T(m)=1./F(m);
    S(m)=dt*(abs(A(m)))^2/(2*pi*N);
end

