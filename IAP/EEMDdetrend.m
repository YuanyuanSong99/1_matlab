clc,clear,close all;
addpath E:\1_matlab\help;
datadir='E:\data\IAP\temperature\Yearly\'; %指定批量数据所在的文件夹
filelist=dir([datadir,'IAP*.h5']); %指定批量数据的类型
h5disp([datadir,filelist(1).name]);
filetemp=[datadir,filelist(1).name]; %查看你要读取的文件的编号。filelist（1）.name在window下为第一个标号的数据
lonData = h5read(filetemp,'/lon');
latData = h5read(filetemp,'/lat');
depthData = h5read(filetemp,'/level');
% save('MatFile/lonData.mat','lonData');
% save('MatFile/latData.mat','latData');
% save('MatFile/depthData.mat','depthData');
% anomaly 
k=length(filelist);

for s=1:k;   
    s
    filename1=[datadir,filelist(s).name];
    Temp(:,:,:,s) = ncread(filename1,'/temp in degC'); %读入变量
end
%
Tempa = double(Temp - nanmean(Temp(:,:,:,42:71),4)); % remove climatology from 1981-2010 
%% nonlinear detrend and filter 8 years each grid
addpath E:\1_matlab\help\EEMD_code\;
Tempa1 = double(permute(Tempa,[4 1 2 3]));
Tempa1 = reshape(Tempa1,[81,360*180*41]);
dT = 1; cf = 1/8; % 8-yr filter
Tempadnl = nan(81,2656800);
Tempadnlf = nan(81,2656800);
parfor i = 1:2656800
    i
    if sum(ismissing(Tempa1(:,i))) == 0 
    Tempeemd = eemd(Tempa1(:,i),0.2,50); % EEMD
    Tempadnl(:,i) = Tempa1(:,i)-Tempeemd(:,7); % nonlinear detrend
    Tempadnlf(:,i) = lanczosfilter(Tempadnl(:,i),dT,cf,[],'low'); % 8 year filtered
    end
end
Tempadnl = permute(reshape(Tempadnl,81,360,180,41),[2 3 4 1]);  % nonlinear detrended      
save('MatFile/Tempadnl.mat','Tempadnl');
Tempadnlf = permute(reshape(Tempadnlf,[81,360,180,41]),[2 3 4 1]);   % nonlinear detrended and 8-year filtered  
save('MatFile/Tempadnlf.mat','Tempadnlf');
