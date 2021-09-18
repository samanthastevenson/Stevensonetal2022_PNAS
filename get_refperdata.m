% Collect data during a reference period, return the running-mean value of
% the data concatenated over the entire ensemble
% Feb 2020
% Sam Stevenson

function [sm_refper,sm_refperdtr]=get_refperdata(runnames,myref,smarr,smarr_dtr)
    global windlen

    sm_refper=[];   % Time series during reference period only
    sm_refperdtr=[];
    for rr=1:length(runnames)        
        sm_refperdtr=cat(1,sm_refperdtr,squeeze(smarr_dtr(rr,myref,:,:)));
        sm_refper=cat(1,sm_refper,squeeze(smarr(rr,myref,:,:)));
    end
    
    % Compute running mean of concatenated reference period data and adjust times such that the running
    % mean at location x is the mean over the preceding 'windlen' timesteps
    for la=1:size(sm_refper,2)
        for lo=1:size(sm_refper,3)                
            smrts=movmean(sm_refper(:,la,lo),windlen,'omitnan');
            sm_refper(:,la,lo)=smrts;       
            
            smrts=movmean(sm_refperdtr(:,la,lo),windlen,'omitnan');
            sm_refperdtr(:,la,lo)=smrts;    
        end
    end 
    
    sm_refper=sm_refper((floor(windlen/2)+1):(end-floor(windlen/2)),:,:);    % Get rid of end effects
    sm_refperdtr=sm_refperdtr((floor(windlen/2)+1):(end-floor(windlen/2)),:,:);    % Get rid of end effects
end