clc,clear,close all;
fname = '/Volumes/CESM-post2/VNTensmean_1920-2005.nc';
ncid=netcdf.open(fname,'NOWRITE');
ncdisp(fname);
TemMMM = ncread(fname,['VNT']); 
