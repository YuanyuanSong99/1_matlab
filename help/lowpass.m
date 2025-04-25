function y=lowpass(x,lpf)
%function y=lowpass(x,lpf) lowpasses time series x.
%all frequencies higher than lpf in cycles per time step
%are removed.

%D. Rudnick May 20, 1994.
if size(x,1)==1
    x=x';
end
[len,num]=size(x);
% span=len-1;
% slope=(x(len,:)-x(1,:))/span;
% off=(x(1,:)*len-x(len,:))/span;
% match=(1:len)'*slope+ones(len,1)*off;
match=x-detrend(x);
xi=fft(x-match);
mlpf=ceil(lpf*len);
xi(mlpf+1:len-mlpf+1,:)=zeros(size(xi(mlpf+1:len-mlpf+1,:)));
y=ifft(xi)+match;
