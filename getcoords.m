% Get filenames, coordinate information for a given large ensemble 
% Feb 2020
% Sam Stevenson

function [mylat,mylon,lat,lon,runnames]=getcoords(ee,varname,dstr,fstr)
    global ensdir
    global ensnames
    global smbox
    
    % Get ensemble data
%     ttmp=strcat('ls -1',{' '},ensdir,ensnames{ee},fstr);
%     [~,list]=system(ttmp{1})
%     runnames= textscan( list, '%s', 'delimiter', '\n' );
%     runnames=runnames{1};
%     runnames=runnames(3:end);
    tmp=dir(strcat(ensdir,ensnames{ee},fstr));
    runnames=cell(length(tmp),1);
    for rr=1:length(runnames)
        runnames{rr}=tmp(rr).name;
    end

    if strcmp(ensnames{ee},'canesm2_lens')
       runnames=runnames([1,3:end]);
    end
    
    % Coordinates
    smnc=netcdf(strcat(ensdir,ensnames{ee},'/',dstr,'/',varname,'/',runnames{1}));
    lat=smnc{'lat'}(:);
    lon=smnc{'lon'}(:);
    mylat=find(lat >= smbox(1) & lat <= smbox(2));
    if smbox(3) == 0 && smbox(4) == 0
        mylon=1:length(lon);
    else
        mylon=find(lon >= smbox(3) & lon <= smbox(4));
    end
    lon(lon > 180)=lon(lon > 180)-360;
    
end