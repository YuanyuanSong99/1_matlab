%--------------------------
% Calculate vorticity
% By Linlin Zhang, IOCAS, Ocotber, 11, 2012.
%--------------------------
function vor=zh_curl(u1,v1,lon,lat)
[u1x,u1y]=sun_grad2(u1,lon,lat);
[v1x,v1y]=sun_grad2(v1,lon,lat);
vor=v1x-u1y;
