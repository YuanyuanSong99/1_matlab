clc,clear,close all;
load('d:/R2014b/bin/2sss.mat');
wp=sss(140:180,74:89,:);  %10S-5N, 140E-180
spcz=sss(160:220,60:74,:);  %24S-10S, 160E-140W
ect=sss(220:230,82:86,:);   %2S-2N,120W-110W
m_wp=nanmean(nanmean(wp));
m_spcz=nanmean(nanmean(spcz));
m_ect=nanmean(nanmean(ect));
w0=permute(m_wp,[3 2 1]);
s0=permute(m_spcz,[3 2 1]);
e0=permute(m_ect,[3 2 1]);
%去趋势（最小二乘法）-------------------------------------
A=[1:1:804];
X=[ones(804,1) A'];
b=regress(w0,X);
c=regress(s0,X);
d=regress(e0,X);
w1=b(1)+b(2)*A;
s1=c(1)+c(2)*A;
e1=d(1)+d(2)*A;
figure(2);
scatter(A,e0);
hold on
plot(A,e1);
w2=w0-w1';
s2=s0-s1';
e2=e0-e1';
%8年低通滤波---------------------------------------------
[w3,mw,hw,Cxw,fw] = lanczosfilter(w2,1,1/(8*12),[],'low');
[s3,ms,hs,Cxs,fs] = lanczosfilter(s2,1,1/(8*12),[],'low');
%8年滑动平均
n = 96;
sm = filter(ones(1,n)/n,1,e2);
sm = sm(n:end);
%figure(2);
%去趋势后直接做小波分析就很好
%小波分析-----------------------------------------------------wavelet
sst = e2;
%------------------------------------------------------ Computation

% normalize by standard deviation (not necessary, but makes it easier
% to compare with plot on Interactive Wavelet page, at
% "http://paos.colorado.edu/research/wavelets/plot/"
variance = std(sst)^2;
sst = (sst - mean(sst))/sqrt(variance) ;

n = length(sst);
dt = 1/12 ;
time = [0:length(sst)-1]*dt + 1951.0 ;  % construct time array
xlim = [1950,2016];  % plotting range
pad = 1;      % pad the time series with zeroes (recommended)
dj = 1/12;    % this will do 4 sub-octaves per octave
s0 = 4*dt;    % this says start at a scale of 6 months
j1 = 7/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.72;  % lag-1 autocorrelation for red noise background
mother = 'Morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum

% Significance levels: (variance=1 for the normalized SST)
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

% Global wavelet spectrum & significance levels:
global_ws = variance*(sum(power')/n);   % time-average over all times
dof = n - scale;  % the -scale corrects for padding at edges
global_signif = wave_signif(variance,dt,scale,1,lag1,-1,dof,mother);

% Scale-average between El Nino periods of 2--8 years
avg = find((scale >= 2) & (scale < 40));
Cdelta = 0.776;   % this is for the MORLET wavelet
scale_avg = (scale')*(ones(1,n));  % expand scale --> (J+1)x(N) array
scale_avg = power ./ scale_avg;   % [Eqn(24)]
scale_avg = variance*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
scaleavg_signif = wave_signif(variance,dt,scale,2,lag1,-1,[2,7.9],mother);

whos

%------------------------------------------------------ Plotting

%--- Plot time series
subplot('position',[0.1 0.65 0.65 0.3])
plot(time,sst)
hold on;

%plot(time,s1,'k');
%legend('原数据','10年滑动平均',1);
set(gca,'XLim',xlim(:))
xlabel('Time (year)')
ylabel('sss')
title('a) 赤道冷舌区SSS去趋势时间序列')

hold off

%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.1 0.65 0.4])
levels = [0.025,0.05,0.1,0.2,0.4,0.8,1.6,3.2,6.4,12.8,25.6] ;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
contourf(time,log2(period),log2(power),log2(levels));  %*** or use 'contourfill'
%imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
colormap('jet')
xlabel('Time (year)')
ylabel('Period (years)')
title('b) SSS Wavelet Power Spectrum')
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks,'FontSize',11)
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(time,log2(period),sig95,[-99,1],'k');
hold on
% cone-of-influence, anything "below" is dubious
plot(time,log2(coi),'k')
%localMaximum
%maxpower=localMaximum(log2(power),[5 5]);

hold off

%--- Plot global wavelet spectrum
subplot('position',[0.77 0.1 0.2 0.4])
plot(global_ws,log2(period))
hold on
plot(global_signif,log2(period),'--')
%localMaximum
x=findpeaks(global_ws);

y(1)=global_ws(46);
y(2)=global_ws(68);
y(3)=global_ws(78);
per(1)=period(46);
per(2)=period(68);
per(3)=period(78);
%text(x(1),log2(y(1)),'0.5','FontSize',15);
%text(y(2),log2(per(2)),'16','FontSize',12);
%text(y(3),log2(per(3)),'29','FontSize',12);
%text(x(4),log2(y(4)),'3','FontSize',15);
%text(x(5),log2(y(5)),'5','FontSize',15);
%text(x(6),log2(y(6)),'9','FontSize',15);



hold off
xlabel('Power ')
title('c) Global Wavelet Spectrum')
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel','')
set(gca,'XLim',[0,1.25*max(global_ws)])



saveas(gcf,'e:/1硕士一年级/气候统计/homework/PPTwork/ECT_SSS去趋势小波分析.png')
% end of code



