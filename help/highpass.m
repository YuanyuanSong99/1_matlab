function y=highpass(x,hpf)
%function y=highpass(x,hpf) highpasses time series x.
%all frequencies lower than hpf in cycles per time step
%are removed.

%D. Rudnick May 23, 1994.
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
mhpf=floor(hpf*len);
xi(1:mhpf+1,:)=zeros(mhpf+1,num);
xi(len-mhpf+1:len,:)=zeros(mhpf,num);
y=real(ifft(xi));
