clear all;
close all;
clc;

% ----- Load the original data, slaa 
nxy = 122*50; nt = 264;
fid = fopen('E:\data5\adt_sla\sla_monthly_low.dat','rb');
sla0 = fread(fid,[nxy,nt],'float');
fclose(fid);
sla0 = reshape(sla0,[122,50,264]);

% %  Choose the study area --- lon = [140.875,160.125]; lat=[31.875,38.125];
% ns = 78*26;
% sla = sla0(45:122,17:42,:);
% data = reshape(sla,[ns,nt]);
% clear sla0 sla

%  Choose the study area ---£¨141.125-153.875E£¬32.125-37.875N£©
ns = 48*24;
sla = sla0(46:93,18:41,:);
data = reshape(sla,[ns,nt]);
clear sla0 sla

% ----- Remove the default value before the SOM analysi
[i_nan,it] = find(isnan(data)); 
for i = 1:length(i_nan)  
    data(i_nan(i),:) = NaN ;
end  
for i = ns:-1:1
    if ( isnan(data(i,1)) ); 
    data(i,:)=[];
    end
end
clear it i

%% SOM analysis
% ----- parameters
% The values which are unambiguous (marked with '*') can be given without the
% preceeding argID.
%   'init'       *(string) initialization: 'randinit' or 'lininit' (default)
%   'algorithm'  *(string) training: 'seq' or 'batch' (default) or 'sompak'
%   'munits'      (scalar) the preferred number of map units
%   'msize'       (vector) map grid size
%   'mapsize'    *(string) do you want a 'small', 'normal' or 'big' map
%                          Any explicit settings of munits or msize override this.
%   'lattice'    *(string) map lattice, 'hexa' or 'rect'
%   'shape'      *(string) map shape, 'sheet', 'cyl' or 'toroid'
%   'neigh'      *(string) neighborhood function, 'gaussian', 'cutgauss',
%                          'ep' or 'bubble'
%   'topol'      *(struct) topology struct
%   'som_topol','sTopol' = 'topol'
%   'mask'        (vector) BMU search mask, size dim x 1
%   'name'        (string) map name
%   'comp_names'  (string array / cellstr) component names, size dim x 1
%   'tracking'    (scalar) how much to report, default = 1
%   'training'    (string) 'short', 'default', 'long'
%                 (vector) size 1 x 2, first length of rough training in epochs, and then length of finetuning in epochs
%   'radius'       neighborhood radiuses: [3 1],[1 1],[0.1 0.01] or others; ... 

% ----- Construct data
sD = data.';
% sD = som_data_struct(D);
clear data

% ----- normalize data
% sD = som_normalize(sD,'var');
% clear D

% ----- initialize the prototype vectors
msize = [2 2]; units = 4;
sM = som_lininit(sD,'msize',msize,'munits',units,'rect','sheet');  

% ----- map training for 5*60 steps
num = 20; bmus = zeros(nt,num); qe = zeros(num,1); te = zeros(num,1);
for m = 1:num 
    
    % --- batch training 
    sM = som_batchtrain(sM, sD,'radius',[1 1],'ep','trainlen',5,'tracking',1);  
          
    % --- the Best-Matching Unit  
    bmu = som_bmus(sM, sD);
    bmus(:,m) = bmu(:,1);
    
    % --- self-organizing map quality
    % Qe: average quantization error between data vectors and their BMUs on the map. 
    % Te: percentage of data vectors for which the first- and second-BMUs are not adjacent units.    
    [Qe,Te] = som_quality(sM, sD);
    qe(m,1) = Qe;
    te(m,1) = Te;

    clear bmu Te Qe

end

save('bmus.mat','bmus'); 

%% Storage the map data 
% Inserte the default values into the original position
% S(nxy,mode) is the final spatial mode.
S = zeros(ns,units); S0 = sM.codebook.';

for i = 1:length(i_nan)
    S(i_nan(i),:) = NaN ;
end
clear i

m = 1;
for i = 1:ns
    if (isnan(S(i,1)));
    else
        S(i,:) = S0(m,:);
        m = m + 1;
    end
end
clear m i i_nan S0

% save('som_units.mat','S');

fid = fopen('som_units.dat','wb');
fwrite(fid,S,'float');
fclose(fid);








