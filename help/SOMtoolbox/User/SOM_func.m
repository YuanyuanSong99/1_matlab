function [S,hits,sM,sD,bmus,Qe,Te]=SOM_func(data,a,b,ab,num)
sD = data';
A=size(sD);
nt = A(1);
ns = A(2);
% sD = som_data_struct(D);
clear data

% ----- normalize data
sD = som_normalize(sD,'var');
% clear D

% ----- initialize the prototype vectors
msize = [a b]; units = ab;
a=1;
sM = som_lininit(sD,'msize',msize,'munits',units,'rect','sheet');

% ----- map training for 5*60 steps
 bmus = zeros(nt,num); qe = zeros(num,1); te = zeros(num,1);
for m = 1:num 
    m;
    
    % --- batch training 
    sM = som_batchtrain(sM, sD,'radius',[1 1],'ep','trainlen',1,'tracking',1);  
          
    % --- the Best-Matching Unit  
    bmu = som_bmus(sM, sD);
    bmus(:,m) = bmu(:,1);
    
    % --- self-organizing map quality
    % Qe: average quantization error between data vectors and their BMUs on the map. 
    % Te: percentage of data vectors for which the first- and second-BMUs are not adjacent units.    
    [Qe,Te] = som_quality(sM, sD);
    qe(m,1) = Qe
    te(m,1) = Te

end

save('H:\station_analyse\SOM\analysis\bmus.mat','bmus'); 

%% Storage the map data 
% Inserte the default values into the original position
% S(nxy,mode) is the final spatial mode.
S = zeros(ns,units); S0 = sM.codebook.';

 m = 1;
 for i = 1:ns
     if (isnan(S(i,1)));
     else
         S(i,:) = S0(m,:);
         m = m + 1;
    end
 end
 clear m i i_nan S0

%fid = fopen('H:\SOM\analysis\som_units.dat','wb');
% fwrite(fid,S,'float');
% fclose(fid);
hits=som_hits(sM,sD)
% U=som_umat(sM)
% Um=U(1:2:end,1:2:end)









