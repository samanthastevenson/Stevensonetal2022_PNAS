% Take result of megadrought property calculation, make something that can
% be added to a composite
% Feb 2020
% Sam Stevenson

function [dr20th,dr20th_dtr,dr21st,dr21st_dtr]=drepochavg(drrisktmp,drrisktmp_dtr,lon,lat,lonfix,latfix)
    % Interpolate epoch averages, differences to common grid
    dr20th=squeeze(nanmean(drrisktmp(:,1,:,:),1));
    dr20th=griddata(lon,lat',dr20th,lonfix,latfix');
    dr21st=squeeze(nanmean(drrisktmp(:,2,:,:),1));
    dr21st=griddata(lon,lat',dr21st,lonfix,latfix');
    dr20th_dtr=squeeze(nanmean(drrisktmp_dtr(:,1,:,:),1));
    dr20th_dtr=griddata(lon,lat',dr20th_dtr,lonfix,latfix');
    dr21st_dtr=squeeze(nanmean(drrisktmp_dtr(:,2,:,:),1));
    dr21st_dtr=griddata(lon,lat',dr21st_dtr,lonfix,latfix'); 
    
    % Reorder array for correct plotting
    myw=find(lonfix < 0);
    mye=find(lonfix >= 0);
    dr20th=dr20th(:,[myw' mye']);
    dr21st=dr21st(:,[myw' mye']);
    dr20th_dtr=dr20th_dtr(:,[myw' mye']);
    dr21st_dtr=dr21st_dtr(:,[myw' mye']);
end