function [obs,omega,rtau,ps_obs,ps_noise,corr,corrp]=con_spectrum(obs,fs,rtau1,preproc,m,noiseind)
%% usage: [omega,ps_obs,ps_noise]=con_spetrum(obs,fs,ratu1,preproc,m,noiseind)
%%input::    obs: observation time series
%%           m  : maximum lag time scale, default value is N/5
%%           noiseind: chosen noise spectrum type; 1, red noise spectrum,
%%              2, white nosie spectrum, 0, judged according to
%%              autocorrelation
%%           fs: sampling time interval
%%           preproc: preprocess indicator 1, anomalized; 2, normalized 
%%           rtau1: critical value for rtau(1)
%% output:: omega: circular frequency
%%          rtau: autocorrelation coefficients
%%          ps_obs: power spetrum of observations, it has two columns, the
%%          first column is gross spectrum, the second column is smoothed
%%          spectrum
%%          ps_noise: power spectrum of noise signals
%%          corr: lag 1 autocorrelation coefficients
%%          corrp: significance level of corr
%% created by Zhai Fangguo on July 7, 2010
error(nargchk(2,6,nargin));
if nargin<3   %% no preproc,m and noiseind
    rtau1=0.001;
    preproc=1;  %% anomalized
    m=5;   %% default value
    noiseind=1;   %%default red spectrum 
elseif nargin<4  %% no m and noiseind
    preproc=1;
    m=5;
    noiseind=1;
elseif nargin<5  %% no noiseind
    m=5;
    noiseind=1;
elseif nargin<6
    noiseind=1;
end

obsnum=length(obs);
m=floor(obsnum/m);  %% maximum lag time scale

%% preprocess the observation data
if preproc==1  %% anomalized
    obs=obs-mean(obs);
    sf=std(obs);
elseif preproc==2   %% normalized
    obs=(obs-mean(obs))/std(obs);
    sf=1;
end

%% determine the circular frequency
omega=pi*[0:m]/(m*fs);

%% calculate the autocorrelation
rtau=nan(m+1,1);
for i=0:m
    rtau(i+1)=sum(obs(1:obsnum-i).*obs(1+i:obsnum))/(obsnum);
end

bl=ones(m+1,1);
bl(1)=1/2;
bl(m+1)=1/2;

%% calculate the observation spectrum
ps_obs=nan(m+1,2);
% first the gross spectrum
for i=0:m
    ps_obs(i+1,1)=rtau(1)+2*sum(rtau(2:m).*cos(pi*i*[1:m-1]'/m))+rtau(m+1)*cos(i*pi);
end
% second the smoothed spectrum
ps_obs(1,2)=(ps_obs(1,1)+ps_obs(2,1))/2;
ps_obs(2:m,2)=ps_obs(1:m-1,1)/4+ps_obs(2:m,1)/2+ps_obs(3:m+1,1)/4;
ps_obs(m+1,2)=(ps_obs(m,1)+ps_obs(m+1,1))/2;

%% the following part is to calculate the noise spectrum for the purpose of
%% testing 
freed=(2*obsnum-3/2*m)/m;   %% freedom degree 
chi005=freed+0.85+1.645*sqrt(2*freed-1);

[r,p]=corrcoef(obs(1:obsnum-1),obs(2:obsnum));
corr=r(1,2);
corrp=p(1,2);

% determine the lag 1 autocorrelation coefficient
if (corr<rtau1&&noiseind==0)||noiseind==2  %% white spectrum
    ps_noise=sf^2*chi005/freed*ones(m+1,1);
elseif (corr>=rtau1&&noiseind==0)||noiseind==1  %% red spectrum
    meane=(ps_obs(1,2)+ps_obs(m+1,2))/2/m+sum(ps_obs(2:m,2))/m;
    ps_noise=meane*chi005/freed*(1-corr^2)./(1+corr^2-2*corr*cos(pi*[0:m]'/m));
end


