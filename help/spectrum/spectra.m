%Spactre function;
%���漴�źŽ��й����׷�������MATLAB�ı�׼FFT��ʽ
%[output,F,T]=spactre(data,dt)
% input :
% data  ------input data vector;
% dt   ------minmum data frequency;
% output :
% S------data spactrum;�����ף��������ܶȣ���λ��mm^2*month (h^2*t) 
% F,T   -------frequency and period serial;
% û�������ż���

function [S,F,T]=spectra(data,dt)

N=max(size(data));  
dt=1/12;       %���¾������ڵ�λ��Ϊ�����ڵ�λ
dw=2*pi/(N*dt);
df=dw/(2*pi);

X=data-mean(data);   %standard
A=fft(X);              %F�ϱ任
for m=1:N
    F(m)=df*m;
    T(m)=1./F(m);
    S(m)=dt*(abs(A(m)))^2/(2*pi*N);
end

