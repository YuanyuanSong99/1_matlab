function [var] = read_1lev_1nc(filepath,varname)
% This function is to read single nc file with single level
    ncid=netcdf.open(filepath,'NOWRITE');
    ncdisp(filepath);
    var = ncread(filepath,varname);
end