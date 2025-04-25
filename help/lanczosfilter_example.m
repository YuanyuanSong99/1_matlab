%---------------------------------------------20-90day bandpass filter,
% 2D
dT=1;
ufilter(1:length(lat),1:length(lon),1:length(L)) =nan;
for ii=1:length(lat);
    for jj=1:length(lon);
    L =find(~isnan(uvel(jj,ii,:)));
    if length(L)
    [a,coef,window,Cx,Ff] = lanczosfilter(squeeze(uvel(jj,ii,L)),dT,1/20,[],'low');
    [ufilter(jj,ii,L),coef,window,Cx,Ff] = lanczosfilter(a,dT,1/90,[],'high');
    clear a coef window Cx Ff;
    end
    clear L
    end
end 

