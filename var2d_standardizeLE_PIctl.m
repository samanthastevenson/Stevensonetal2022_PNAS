% Read in PI control data for a large ensemble, for standardization
% purposes
% Feb 2020
% Sam Stevenson

function [pimn,pimn_dtr,pistd,pistd_dtr]=var2d_standardizeLE_PIctl(ctlnames,smbox,ctrlstrtyr,windlen,varname,fac)

    pictrlsm=[];
    piyr=[];
    for cc=1:length(ctlnames)
        ctlname=ctlnames{cc}
        nc=netcdf(ctlname);
        PItime=nc{'time'}(:);
        lat=nc{'lat'}(:);
        lon=nc{'lon'}(:);
        lon(lon > 180)=lon(lon > 180)-360;
        [PIyr,PImon,~]=datenumnoleap(PItime,[ctrlstrtyr 1 1]);
        uyr=unique(PIyr);
        PIyr(1)
        PImon(1)
        PIyr(end)
        PImon(end)

        if (smbox(3) >= 0 && smbox(4) >= 0) || (smbox(3) <= 0 && smbox(4) <= 0)
            mylat=find(lat >= smbox(1) & lat <= smbox(2));
            if smbox(3) == 0 && smbox(4) == 0
                mylon=1:length(lon);
            else
                mylon=find(lon >= smbox(3) & lon <= smbox(4));
            end
            sm=nc{varname}(:,mylat,mylon);
        else
            mylat=find(lat >= smbox(1) & lat <= smbox(2));
            mylonw=find(lon >= smbox(3) & lon < 0);
            mylone=find(lon >= 0 & lon <= smbox(4));
            sme=nc{varname}(:,mylat,mylone);
            smw=nc{varname}(:,mylat,mylonw);
            sm=cat(3,sme,smw);
        end
        
        sm=sm*fac;  % correction factor for unit conversion if necessary
        sm(abs(sm) > 1e10)=0/0;

        % Annual mean
        piannmn=zeros(length(uyr),size(sm,2),size(sm,3));
        for yy=1:length(uyr)
            myy=find(PIyr == uyr(yy));
            piannmn(yy,:,:)=squeeze(nanmean(sm(myy,:,:),1));
        end

        % Store annual mean soil moisture array
        pictrlsm=cat(1,pictrlsm,piannmn);   
        piyr=cat(1,piyr,uyr');
        
    end
        
    % Make "detrended" version for later use
    pictrlsm_dtr=pictrlsm-repmat(nanmean(pictrlsm,1),[size(pictrlsm,1) 1 1]);

    % Compute running mean of anomalies and adjust times such that the running
    % mean at location x is the mean over the preceding 'windlen' timesteps
    for la=1:size(pictrlsm,2)
        for lo=1:size(pictrlsm,3)
            smts=movmean(pictrlsm(:,la,lo),windlen,'omitnan');
            pictrlsm(:,la,lo)=smts;             

            smts=runmean(pictrlsm_dtr(:,la,lo),floor(windlen/2));
            pictrlsm_dtr(:,la,lo)=smts;     
        end
    end 
    pictrlsm=pictrlsm((floor(windlen/2)+1):(end-floor(windlen/2)),:,:);    % Get rid of end effects
    pictrlsm_dtr=pictrlsm_dtr((floor(windlen/2)+1):(end-floor(windlen/2)),:,:);    % Get rid of end effects

    % Store mean, std. dev. values for standardization
    pimn=nanmean(pictrlsm,1);
    pimn_dtr=nanmean(pictrlsm_dtr,1);
    pistd=nanstd(pictrlsm,1);
    pistd_dtr=nanstd(pictrlsm_dtr,1);
    pimn(abs(pimn) < 1e-10)=0/0;
    pimn(abs(pimn_dtr) < 1e-10)=0/0;
    pistd(abs(pistd) < 1e-10)=0/0;
    pistd_dtr(abs(pistd_dtr) < 1e-10)=0/0;

end