%WAVETEST Example Matlab script for WAVELET, using NINO3 SST dataset
addpath G:\1_matlab\help\wave_matlab\;

sst = DIzd_0; str = 'PmAI';

%------------------------------------------------------ Computation
variance = std(sst)^2;
sst = (sst - mean(sst))/sqrt(variance) ;

n = length(sst);
dt = 1 ;
time = [0:length(sst)-1]*dt + 1940.0 ;  % construct time array
xlim = [1940,2020];  % plotting range
pad = 1;      % pad the time series with zeroes (recommended)
dj = 1/4;    % this will do 4 sub-octaves per octave
s0 = 2*dt;    % 
j1 = -1;    % this says do 7 powers-of-two with dj sub-octaves each
cor1 = corrcoef(sst(1:end-1),sst(2:end));
lag1 = cor1(1,2);  % lag-1 autocorrelation for red noise background
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
close all;
Fig = figure('position',[10 50 800 600]);
%--- Plot time series
subplot('position',[0.1 0.65 0.65 0.3])
plot(time,sst)
hold on;

set(gca,'XLim',xlim(:),'FontSize',11)
xlabel('Time (year)')
ylabel('Normalized T_7_0_0')
title(['a) ',str])

%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.1 0.65 0.4])
levels = [0.025,0.05,0.1,0.2,0.4,0.8,1.6,3.2,6.4,12.8,25.6] ;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
contourf(time,log2(period),log2(power),log2(levels));  %*** or use 'contourfill'
%imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
colormap('jet')
xlabel('Time (year)')
ylabel('Period (years)')
title('b) Wavelet Power Spectrum')
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

%--- Plot global wavelet spectrum
subplot('position',[0.77 0.1 0.2 0.4])
plot(global_ws,log2(period))
hold on
plot(global_signif,log2(period),'--')
%localMaximum
%ans=findpeaks(global_ws);
x=[0.3346    0.1484    0.2113    0.2114    0.5477    0.4187];
y=[0.5159    1.0319    2.1865    3.0921    5.2003    9.2659];
%scatter(x,log2(y),'k','o');
% text(x(1),log2(y(1)),'0.5','FontSize',15);
% text(x(2),log2(y(2)),'1','FontSize',15);
% text(x(3),log2(y(3)),'2','FontSize',15);
% text(x(4),log2(y(4)),'3','FontSize',15);
% text(x(5),log2(y(5)),'5','FontSize',15);
% text(x(6),log2(y(6)),'9','FontSize',15);

xlabel('Power ')
title('c) Global Wavelet Spectrum')
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel','','FontSize',11)
set(gca,'XLim',[0,1.25*max(global_ws)])

print(Fig,['G:\figures\IAP\Yearly\20230911_IPO_SouthernOcean_46_61S\wavelet_',str,'.png'],'-dpng','-r300')


