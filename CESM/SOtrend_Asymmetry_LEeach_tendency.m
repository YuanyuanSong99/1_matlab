% Data 
clc,clear,close all;
addpath(genpath('/Volumes/Togo4T/1_matlab/help'));
load("MatFile/lonData.mat");
load("MatFile/latData.mat");
load("MatFile/depthData.mat");
nlon = length(lonData); nlat = length(latData);
nlev = length(depthData);
cp = 4096 % J (kg oC)
ro = 1025 % kg m-3
depthstr = '0-700m';
lats_int = [35:55]; % 55S-35S
lats_so = [1:60]; % 90S-30S
intlev = 37; intdep = 700; % 700m
dy = 111*1000; % each latitude
clear dx
for j = 1:length(latData(lats_so));
    dx(1,j) = sw_dist([latData(j),latData(j)],[lonData(1),lonData(2)],'km')*1000; % 每个纬度上，经度之间距离不一样
end
% dweit = depthData(2:intlev-1)-depthData(1:intlev-2);
levs = [1:intlev];
deps = [depthData(1:intlev-1);intdep];
dz = gradient(deps);
%% tendency
datadir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*192001-200512.nc']); %指定批量数据的类型
% rcp85
datadir5='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
filelist5=dir([datadir5,'*200601*.nc']); %指定批量数据的类型
clear s*_0
for s=1:40
    s 
    filename1=[datadir1,filelist1(s).name];
    ncid=netcdf.open(filename1,'NC_NOWRITE');
    Temp1 = ncread(filename1,'TEMP');
    filename5=[datadir5,filelist5(s).name];
    ncid=netcdf.open(filename5,'NC_NOWRITE');
    Temp2 = ncread(filename5,'TEMP');  
    Temp = cat(4,Temp1(:,lats_so,1:intlev,41:86),Temp2(:,lats_so,1:intlev,1:16));
    clear Temp1 Temp2
    % 0-700m OHC
    % VV = ones(length(lonData),length(latData),length(depthData)); % calculate volume
    Tempint = Temp(:,:,intlev-1,:)+(Temp(:,:,intlev,:)-Temp(:,:,intlev-1,:))*(intdep-depthData(intlev-1))/(depthData(intlev)-depthData(intlev-1));
    Tempz = cat(3,Temp(:,:,1:intlev-1,:),Tempint).*permute(dz,[3 2 1]);
    % VVint = VV(:,:,intlev-1,:)+(VV(:,:,intlev,:)-VV(:,:,intlev-1,:))*(intdep-depthData(intlev-1))/(depthData(intlev)-depthData(intlev-1));
    % VVz = cat(3,VV(:,:,1)*5,VV(:,:,2:intlev-1).*(permute(dweit,[3 2 1])),VVint*(intdep-depthData(intlev-1)));
    Tempzy = Tempz*dy; % 111 km / latitude
    % VVzy = VVz*111*1000;
    clear Tempzyx VVzyx
    Tempzyx = Tempzy.*dx; 
    save(['/Volumes/Togo4T/1_matlab/CESM/MatFile/OHC_0_700m_30S/Tempzyx_',num2str(s)],'Tempzyx');
    % VVzyx = VVzy.*dx*1000; % m
    Tzyxsub_r_0 = cat(1,Tempzyx(151:360,:,:,:),Tempzyx(1:150,:,:,:));
    % VVzyx(find(isnan(Tempzyx(:,:,:,1)) == 1)) = nan;
    % Vzyxsub_r_0 = cat(1,VVzyx(151:360,:,:),VVzyx(1:150,:,:));
    spac_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:150,lats_int,levs,:),1),2),3));
    % spacV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(1:150,lats,levs),1),2),3));
    % spac0V = spac_0./spacV;
    sia_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(151:360,lats_int,levs,:),1),2),3));
    % siaV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(151:360,lats,levs),1),2),3));
    % sia0V = sia_0./siaV;
    % sdiff = spac0V-sia0V;
    sat_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(151:230,lats_int,levs,:),1),2),3));
    % satV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(151:230,lats,levs),1),2),3));
    % sat0V = sat_0./satV;
    sio_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(231:360,lats_int,levs,:),1),2),3));
    % sioV = cp*ro*squeeze(nansum(nansum(nansum(Vzyxsub_r_0(231:360,lats,levs),1),2),3));
    % sio0V = sio_0./sioV;
    sall_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(:,lats_int,levs,:),1),2),3));
end
clear Temp* Tzyxsub_r_0
%% Tendency of OHC (minus 1960) ZJ
for s = 1:40
    load(['MatFile/OHC_0_700m_30S/Tempzyx_',num2str(s),'.mat'])
    Tzyxsub_r_0 = cat(1,Tempzyx(151:360,:,:,:),Tempzyx(1:150,:,:,:));
    spac_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(1:150,lats_int,levs,:),1),2),3));
    sia_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(151:360,lats_int,levs,:),1),2),3));
    sat_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(151:230,lats_int,levs,:),1),2),3));
    sio_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(231:360,lats_int,levs,:),1),2),3));
    sall_0(:,s) = cp*ro*squeeze(nansum(nansum(nansum(Tzyxsub_r_0(:,lats_int,levs,:),1),2),3));
end
figure('Position',[100 100 800 500]);
hold on
for s = 1:40    
    tdpac(:,s) = (spac_0(:,s)-spac_0(1,s))/10^21;
    tdia(:,s) = (sia_0(:,s)-sia_0(1,s))/10^21;
    tdat(:,s) = (sat_0(:,s)-sat_0(1,s))/10^21;
    tdio(:,s) = (sio_0(:,s)-sio_0(1,s))/10^21;
    tdall(:,s) = (sall_0(:,s)-sall_0(1,s))/10^21;
    plot(tdia(:,s));
end
%% heat forcing (minus 1960) ZJ
% NHFLX
datadir1='/Volumes/CESM-post2/CESM-post/LE/FLNS/yearly/'; %指定批量数据所在的文件夹
filelist1=dir([datadir1,'*192001-200512.nc']); %指定批量数据的类型
datadir2='/Volumes/CESM-post2/CESM-post/LE/FSNS/yearly/'; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*192001-200512.nc']); %指定批量数据的类型
datadir3='/Volumes/CESM-post2/CESM-post/LE/LHFLX/yearly/'; %指定批量数据所在的文件夹
filelist3=dir([datadir3,'*192001-200512.nc']); %指定批量数据的类型
datadir4='/Volumes/CESM-post2/CESM-post/LE/SHFLX/yearly/'; %指定批量数据所在的文件夹
filelist4=dir([datadir4,'*192001-200512.nc']); %指定批量数据的类型
% rcp85
datadir5='/Volumes/CESM-post2/CESM-post/LE/FLNS/yearly/'; %指定批量数据所在的文件夹
filelist5=dir([datadir5,'*200601*.nc']); %指定批量数据的类型
datadir6='/Volumes/CESM-post2/CESM-post/LE/FSNS/yearly/'; %指定批量数据所在的文件夹
filelist6=dir([datadir6,'*200601*.nc']); %指定批量数据的类型
datadir7='/Volumes/CESM-post2/CESM-post/LE/LHFLX/yearly/'; %指定批量数据所在的文件夹
filelist7=dir([datadir7,'*200601*.nc']); %指定批量数据的类型
datadir8='/Volumes/CESM-post2/CESM-post/LE/SHFLX/yearly/'; %指定批量数据所在的文件夹
filelist8=dir([datadir8,'*200601*.nc']); %指定批量数据的类型
clear hfc*
for s=1:40
    s 
    filename1=[datadir1,filelist1(s).name];
    FLNS1 = ncread(filename1,'FLNS');
    filename5=[datadir5,filelist5(s).name];
    FLNS2 = ncread(filename5,'FLNS');  
    FLNS = cat(3,FLNS1(:,lats_so,41:86),FLNS2(:,lats_so,2:16)); % strange, it should start from 1. But 1 shows abrupt decrease.
    % clear FLNS1 FLNS2
    filename2=[datadir2,filelist2(s).name];
    FSNS1 = ncread(filename2,'FSNS');
    filename6=[datadir6,filelist6(s).name];
    FSNS2 = ncread(filename6,'FSNS');  
    FSNS = cat(3,FSNS1(:,lats_so,41:86),FSNS2(:,lats_so,2:16));
    clear FSNS1 FSNS2
    filename3=[datadir3,filelist3(s).name];
    LHFLX1 = ncread(filename3,'LHFLX');
    filename7=[datadir7,filelist7(s).name];
    LHFLX2 = ncread(filename7,'LHFLX');  
    LHFLX = cat(3,LHFLX1(:,lats_so,41:86),LHFLX2(:,lats_so,2:16));
    clear LHFLX1 LHFLX2
    filename4=[datadir4,filelist4(s).name];
    SHFLX1 = ncread(filename4,'SHFLX');
    filename8=[datadir8,filelist8(s).name];
    SHFLX2 = ncread(filename8,'SHFLX');  
    SHFLX = cat(3,SHFLX1(:,lats_so,41:86),SHFLX2(:,lats_so,2:16));
    clear SHFLX1 SHFLX2
    nhflx = FSNS-SHFLX-LHFLX-FLNS;
    hfc = nhflx.*dx*dy;
    hfc_r = cat(1,hfc(151:360,:,:),hfc(1:150,:,:));
    hfcpac(:,s) = squeeze(sum(sum(hfc_r(1:150,lats_int,:),1),2));
    hfcia(:,s) = squeeze(sum(sum(hfc_r(151:360,lats_int,:),1),2));
    hfcat(:,s) = squeeze(sum(sum(hfc_r(151:230,lats_int,:),1),2));
    hfcio(:,s) = squeeze(sum(sum(hfc_r(231:360,lats_int,:),1),2));
    hfcall(:,s) = squeeze(sum(sum(hfc_r(:,lats_int,:),1),2));
end
clear FLNS FSNS LHFLX SHFLX 
for s = 1:40
    focpac(:,s) = cumsum(hfcpac(:,s))*60*60*24*365/10^21;
    focia(:,s) = cumsum(hfcia(:,s))*60*60*24*365/10^21;
    focat(:,s) = cumsum(hfcat(:,s))*60*60*24*365/10^21;
    focio(:,s) = cumsum(hfcio(:,s))*60*60*24*365/10^21;
    focall(:,s) = cumsum(hfcall(:,s))*60*60*24*365/10^21;
end
%%
close all
figure(1)
hold on
for s = 1:40
    plot(focio(:,s));
end
%%
close all
Fig = figure('position',[700 100 800 400]);
plot(tdio,'r')
hold on
plot(focio,'b')
set(gca,'XLim',[1,61],'ylim',[-100,100],'ytick',[-100:50:100]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
ylabel('ZJ')
legend('tendency','forcing','location','north')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240508_SO_Warming/Tendency_IO_',depthstr,'_35S_55S.png'],'-dpng','-r300')

%% advection
TEMPdir1='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
TEMPlist1=dir([TEMPdir1,'*192001-200512.nc']); %指定批量数据的类型
TEMPdir2='/Volumes/CESM-post2/CESM-post/LE/TEMP/yearly/'; %指定批量数据所在的文件夹
TEMPlist2=dir([TEMPdir2,'*200601*.nc']); %指定批量数据的类型
[UVELdir1, UVELlist1, UISOPdir1, UISOPlist1, USUBMdir1, USUBMlist1]=velhisdir('U');
[UVELdir2, UVELlist2, UISOPdir2, UISOPlist2, USUBMdir2, USUBMlist2]=velrcpdir('U');
[VVELdir1, VVELlist1, VISOPdir1, VISOPlist1, VSUBMdir1, VSUBMlist1]=velhisdir('V');
[VVELdir2, VVELlist2, VISOPdir2, VISOPlist2, VSUBMdir2, VSUBMlist2]=velrcpdir('V');
[WVELdir1, WVELlist1, WISOPdir1, WISOPlist1, WSUBMdir1, WSUBMlist1]=velhisdir('W');
[WVELdir2, WVELlist2, WISOPdir2, WISOPlist2, WSUBMdir2, WSUBMlist2]=velrcpdir('W');
splon = 151; splat = 60; 
lonData_r = cat(1,lonData(splon:360),lonData(1:splon-1));
for s = 1:40
    s
    TEMP1 = ncread([TEMPdir1,TEMPlist1(s).name],'TEMP'); 
    TEMP2 = ncread([TEMPdir2,TEMPlist2(s).name],'TEMP'); 
    TEMP = cat(4,TEMP1(:,:,1:intlev,41:86),TEMP2(:,:,1:intlev,1:16)); 
    clear TEMP1 TEMP2
    Tempint = TEMP(:,:,intlev-1,:)+(TEMP(:,:,intlev,:)-TEMP(:,:,intlev-1,:))*(intdep-depthData(intlev-1))/(depthData(intlev)-depthData(intlev-1));
    Temp0_700 = cat(3,TEMP(:,:,1:intlev-1,:),Tempint);
    clear TEMP
    [dT1 dT2 dT3] = gradient(Temp0_700(:,1:splat,:,:));
    dTdy = dT1./dy;
    dTdx = dT2./dx;
    dTdz = dT3./permute(-dz,[3 2 1]);
    clear Temp0_700 dT1 dT2 dT3
    UVEL1 = ncread([UVELdir1,UVELlist1(s).name],'UVEL'); 
    UVEL2 = ncread([UVELdir2,UVELlist2(s).name],'UVEL'); 
    UVEL = 0.01*cat(4,UVEL1(:,:,1:intlev,41:86),UVEL2(:,:,1:intlev,1:16)); % cm -> m
    clear UVEL1 UVEL2
    UISOP1 = ncread([UISOPdir1,UISOPlist1(s).name],'UISOP'); 
    UISOP2 = ncread([UISOPdir2,UISOPlist2(s).name],'UISOP'); 
    UISOP = 0.01*cat(4,UISOP1(:,:,1:intlev,41:86),UISOP2(:,:,1:intlev,1:16)); % cm -> m
    clear UISOP1 UISOP2
    USUBM1 = ncread([USUBMdir1,USUBMlist1(s).name],'USUBM'); 
    USUBM2 = ncread([USUBMdir2,USUBMlist2(s).name],'USUBM'); 
    USUBM = 0.01*cat(4,USUBM1(:,:,1:intlev,41:86),USUBM2(:,:,1:intlev,1:16)); % cm -> m
    clear USUBM1 USUBM2
    URES = UVEL+UISOP+USUBM; 
    clear UVEL UISOP USUBM

    VVEL1 = ncread([VVELdir1,VVELlist1(s).name],'VVEL'); 
    VVEL2 = ncread([VVELdir2,VVELlist2(s).name],'VVEL'); 
    VVEL = 0.01*cat(4,VVEL1(:,:,1:intlev,41:86),VVEL2(:,:,1:intlev,1:16)); % cm -> m
    clear VVEL1 VVEL2
    VISOP1 = ncread([VISOPdir1,VISOPlist1(s).name],'VISOP'); 
    VISOP2 = ncread([VISOPdir2,VISOPlist2(s).name],'VISOP'); 
    VISOP = 0.01*cat(4,VISOP1(:,:,1:intlev,41:86),VISOP2(:,:,1:intlev,1:16)); % cm -> m
    clear VISOP1 VISOP2
    VSUBM1 = ncread([VSUBMdir1,VSUBMlist1(s).name],'VSUBM'); 
    VSUBM2 = ncread([VSUBMdir2,VSUBMlist2(s).name],'VSUBM'); 
    VSUBM = 0.01*cat(4,VSUBM1(:,:,1:intlev,41:86),VSUBM2(:,:,1:intlev,1:16)); % cm -> m
    clear VSUBM1 VSUBM2
    VRES = VVEL+VISOP+VSUBM; 
    clear VVEL VISOP VSUBM

    WVEL1 = ncread([WVELdir1,WVELlist1(s).name],'WVEL'); 
    WVEL2 = ncread([WVELdir2,WVELlist2(s).name],'WVEL'); 
    WVEL = 0.01*cat(4,WVEL1(:,:,1:intlev,41:86),WVEL2(:,:,1:intlev,1:16)); % cm -> m
    clear WVEL1 WVEL2
    WISOP1 = ncread([WISOPdir1,WISOPlist1(s).name],'WISOP'); 
    WISOP2 = ncread([WISOPdir2,WISOPlist2(s).name],'WISOP'); 
    WISOP = 0.01*cat(4,WISOP1(:,:,1:intlev,41:86),WISOP2(:,:,1:intlev,1:16)); % cm -> m
    clear WISOP1 WISOP2
    WSUBM1 = ncread([WSUBMdir1,WSUBMlist1(s).name],'WSUBM'); 
    WSUBM2 = ncread([WSUBMdir2,WSUBMlist2(s).name],'WSUBM'); 
    WSUBM = 0.01*cat(4,WSUBM1(:,:,1:intlev,41:86),WSUBM2(:,:,1:intlev,1:16)); % cm -> m
    clear WSUBM1 WSUBM2
    WRES = WVEL+WISOP+WSUBM; 
    clear WVEL WISOP WSUBM
    advT = URES(:,1:splat,:,:).*dTdx + VRES(:,1:splat,:,:).*dTdy + WRES(:,1:splat,:,:).*dTdz;
    clear URES dTdx VRES dTdy WRES dTdz
    advTV = ro*cp*advT .*dx *dy .*permute(dz,[3 2 1]);
    advTVr = cat(1,advTV(splon:360,:,:,:),advTV(1:splon-1,:,:,:));
    adp(:,s) = squeeze(nansum(nansum(nansum(advTVr(1:150,lats_int,:,:),1),2),3)); % 60W
    adia(:,s) = squeeze(nansum(nansum(nansum(advTVr(151:360,lats_int,:,:),1),2),3));
    adat(:,s) = squeeze(nansum(nansum(nansum(advTVr(151:230,lats_int,:,:),1),2),3));
    adio(:,s) = squeeze(nansum(nansum(nansum(advTVr(231:360,lats_int,:,:),1),2),3));
    adall(:,s) = squeeze(nansum(nansum(nansum(advTV(:,lats_int,:,:),1),2),3));
end

%%
for s = 1:35
    advpac(:,s) = 60*60*24*365*cumsum(adp(:,s))/10^21;
    advia(:,s) = 60*60*24*365*cumsum(adia(:,s))/10^21;
    advat(:,s) = 60*60*24*365*cumsum(adat(:,s))/10^21;
    advio(:,s) = 60*60*24*365*cumsum(adio(:,s))/10^21;
    advall(:,s) = 60*60*24*365*cumsum(adall(:,s))/10^21;
end
%%
close all
figure(1)
hold on
for s = 1:35
    plot(adp(:,s))
end
%% tendency 
% regstr = 'Pacific';
% ln1 = mean(tdpac,2); ln2 = mean(focpac,2); ln3 = -mean(advpac,2); ln4 = ln2(1:61)+ln3(1:61);

% regstr = 'Atlantic-Indian Ocean';
% ln1 = mean(tdia,2); ln2 = mean(focia,2); ln3 = -mean(advia,2); ln4 = ln2(1:61)+ln3(1:61);

% regstr = 'Atlantic';
% ln1 = mean(tdat,2); ln2 = mean(focat,2); ln3 = -mean(advat,2); ln4 = ln2(1:61)+ln3(1:61);
% 
% regstr = 'Indian Ocean';
% ln1 = mean(tdio,2); ln2 = mean(focio,2); ln3 = -mean(advio,2); ln4 = ln2(1:61)+ln3(1:61);
% 
regstr = 'Zonalbelt';
ln1 = mean(tdall,2); ln2 = mean(focall,2); ln3 = -mean(advall,2); ln4 = ln2(1:61)+ln3(1:61);

close all
Fig = figure('position',[700 100 800 400]);
plot(ln1,'r')
hold on
plot(ln2,'b')
plot(ln3,'g')
plot(ln4,'r--')
set(gca,'XLim',[1,61],'ylim',[-150,150],'ytick',[-150:30:150]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
ylabel('ZJ')
legend('tendency','forcing','transport','forcing+transport','location','northwest')
legend('boxoff')
print(Fig,['/Users/yysong/Desktop/figures/CESM/Yearly/LE/0_LEM/20240513_SO_warming/tendency_',regstr,'.png'],'-dpng','-r300')

%%
close all

Fig = figure('position',[700 100 800 400]);
plot(tdia,'r')
hold on
plot(focia,'b')
plot(advia,'g')
plot(focia(1:46)+advia','k')
set(gca,'XLim',[1,61],'ylim',[-90,90],'ytick',[-90:30:90]);
set(gca,'XTick',[1:10:61]);
set(gca,'XTickLabel',[1960:10:2020],'FontSize',14);
ylabel('ZJ')
legend('tendency','forcing','transport','forcing+transport','location','north')
legend('boxoff')

function [datadir2, filelist2, datadir3, filelist3, datadir4, filelist4]=velhisdir(velname)
datadir2=['/Volumes/CESM-post2/CESM-post/LE/',velname,'VEL/yearly/']; %指定批量数据所在的文件夹
filelist2=dir([datadir2,'*192001-200512.nc']); %指定批量数据的类型
datadir3=['/Volumes/CESM-post2/CESM-post/LE/',velname,'ISOP/yearly/']; %指定批量数据所在的文件夹
filelist3=dir([datadir3,'*192001-200512.nc']); %指定批量数据的类型
datadir4=['/Volumes/CESM-post2/CESM-post/LE/',velname,'SUBM/yearly/']; %指定批量数据所在的文件夹
filelist4=dir([datadir4,'*192001-200512.nc']); %指定批量数据的类型
end
function [datadir6, filelist6, datadir7, filelist7, datadir8, filelist8]=velrcpdir(velname)
datadir6=['/Volumes/CESM-post2/CESM-post/LE/',velname,'VEL/yearly/']; %指定批量数据所在的文件夹
filelist6=dir([datadir6,'*200601*.nc']); %指定批量数据的类型
datadir7=['/Volumes/CESM-post2/CESM-post/LE/',velname,'ISOP/yearly/']; %指定批量数据所在的文件夹
filelist7=dir([datadir7,'*200601*.nc']); %指定批量数据的类型
datadir8=['/Volumes/CESM-post2/CESM-post/LE/',velname,'SUBM/yearly/']; %指定批量数据所在的文件夹
filelist8=dir([datadir8,'*200601*.nc']); %指定批量数据的类型
end
function varMMM = loadsurface(filename1,filename2,varstr)
ncid=netcdf.open(filename1,'NOWRITE');
ncdisp(filename1);
varMMM1 = ncread(filename1,varstr);
varMMM1(:,:,87) = [];
lonData = ncread(filename1,'lon');
latData = ncread(filename1,'lat');

ncid=netcdf.open(filename2,'NOWRITE');
ncdisp(filename2);
varMMM2 = ncread(filename2,varstr);

varMMM = double(cat(3,varMMM1(:,:,1:86),varMMM2(:,:,1:75)));
end
function trd = trend_cal_3D(var)
% var is time*d1*d2
    x = [1:size(var,1)]';
    clear trd
    for i = 1:size(var,2);
        for j = 1:size(var,3);
            par=polyfit(x,var(:,i,j),1); % regression parameters
            trd(i,j) = par(1);
        end
    end
end
function [ts_zs ts] = areamean(var,lons,lats,latData);
    var1 = var(lons,lats,:); 
    var2 = var(lons,lats,1);
    var2(find(isnan(var2) == 0)) = 1; % weight
    ts = reshape(nansum(nansum((cos(latData(lats)'/180*pi)).*var1(:,:,:),1),2)/nansum(cos(latData(lats)'/180*pi).*var2,'all'),size(var1,3),1);
    ts_zs = zscore(ts);
end