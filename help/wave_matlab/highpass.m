function y = highpass(x,Fc,Fs,n);
% HIGHPASS Zero-phase shift, high-pass digital filter.
%
% y = highpass(x,Fc,Fs,n) applies a single pole, zero phase shift, high-pass 
% filter with a 3-db point at Fc to the input signal x. The sampling frequency
% Fs and the filter order n are optional input parameters. If x is a matrix, 
% the function works on the columns of x.
%
% Input parameters:
%    x  input vector or matrix
%    Fc corner frequency 截止频率w=pi/T
%    Fs sampling frequency [Hz, default 256]样本频率f=1/样本数
%    n  filter order [default 1]
%
% Output parameters:
%    y  filtered vector or matrix


% Fabian Wolk 2001/02/21
% 2001 Alec Electronics Co., Ltd.
% #Revision 0.0: 2001/08/09#

error(nargchk(2,4,nargin));
if nargin<4, n = 1; end
if nargin<3, Fs = TMFS;end
if isempty(Fs), Fs = TMFS; end

Fny = Fs/2; % Nyquist frequency
[b,a] = butter(n,Fc/Fny,'high');
[r,c] = size(x);
y = zeros(r,c);
for i = 1:c
   y(:,i) = filtfilt(b,a,x(:,i));
end
