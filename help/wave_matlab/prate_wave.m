%WAVETEST Example Matlab script for WAVELET, using NINO3 SST dataset
%
% See "http://paos.colorado.edu/research/wavelets/"
% Written January 1998 by C. Torrence
%
% Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
%   changed all "log" to "log2", changed logarithmic axis on GWS to
%   a normal axis.
clc,clear,close all;
ncid=netcdf.open('F:/graduation/data/monthlyprecipitation1950-2016.nc','NOWRITE');
ncdisp('F:/graduation/data/monthlyprecipitation1950-2016.nc');
prate=ncread('F:/graduation/data/monthlyprecipitation1950-2016.nc','PRATE');
lon=ncread('F:/graduation/data/monthlyprecipitation1950-2016.nc','LONN95_97');
lat=ncread('F:/graduation/data/monthlyprecipitation1950-2016.nc','LAT');
p1=prate(161:193,42:53,:);%暖池区
p11=nanmean(nanmean(p1,1),2);
p12=[1:1:804];
p12(:)=p11(1,1,:);
%8年低通滤波
%[p1,c,h,Cx,f] = lanczosfilter(p,1,1/(8*12),[],'low');


%n = 96;
%sm = filter(ones(1,n)/n,1,s);
%sm = sm(n:end);
%plot(sm);


%qda2=detrend(sm);
 


sst = p12;

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
subplot('position',[0.1 0.75 0.65 0.2])
plot(time,sst)
set(gca,'XLim',xlim(:))
xlabel('Time (year)')
ylabel('降水率')
title('a) 西赤道太平洋暖池降水率时间序列')
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
title('b) 降水 Wavelet Power Spectrum')
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
maxpower=localMaximum(log2(power),[5 5]);
%x=[1.9643e+03 1.9973e+03 2.0075e+03];
%y=[9.2659 24.7371 11.0191];
%hold on
%scatter(x,log2(y),'+');
hold off

%--- Plot global wavelet spectrum
subplot('position',[0.77 0.1 0.2 0.4])
plot(global_ws,log2(period))
hold on
plot(global_signif,log2(period),'--')
%localMaximum
x=[2.3515e-10 2.5603e-10];
y=[0.4870 0.9740];

text(x(1),log2(y(1)),'0.5','FontSize',15);
text(x(2),log2(y(2)),'1','FontSize',15);


hold off
xlabel('Power ')
title('c) Global Wavelet Spectrum')
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel','')
set(gca,'XLim',[0,1.25*max(global_ws)])



saveas(gcf,'f:/graduation/data/EN4/西赤道太平洋暖池降水率直接小波分析.png')
%saveas(gcf,'f:/graduation/data/EN4/热带太平洋sss直接小波分析.png')
% end of code

