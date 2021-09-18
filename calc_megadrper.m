% Given an input soil moisture time series, calculate the persistence of
% megadrought in that time series
% February 2020
% Sam Stevenson

function [drpertmp,drpertmp_dtr]=calc_megadrper(uyrfilt,smtmp,smtmp_dtr,mylat,mylon)
    global cstrt
    global cend
    global thr
    
    drpertmp=NaN(length(cstrt),length(mylat),length(mylon));
    drpertmp_dtr=NaN(length(cstrt),length(mylat),length(mylon));

    for cc=1:length(cstrt)
        for la=1:length(mylat)
            for lo=1:length(mylon)
                if thr < 0
                    % Drought
                    [~,drlen,~,drloc]=tsper(squeeze(smtmp(:,la,lo)),thr);    
                    [~,drlen_dtr,~,drloc_dtr]=tsper(squeeze(smtmp_dtr(:,la,lo)),thr);    
                else
                    % Pluvial
                    [drlen,~,drloc,~]=tsper(smtmp(:,la,lo),thr);    
                    [drlen_dtr,~,drloc_dtr,~]=tsper(smtmp_dtr(:,la,lo),thr);    
                end
                if min(drloc) > 0
                    timmd=uyrfilt(drloc);
                    drend=timmd+drlen;
                    mydrt=find(timmd <= cend(cc) & drend >= cstrt(cc));
                    if ~isempty(mydrt)
                        drpertmp(cc,la,lo)=nanmean(drlen(mydrt));  
                    end
                end
                if min(drloc_dtr) > 0
                    timmd_dtr=uyrfilt(drloc_dtr);
                    drend_dtr=timmd_dtr+drlen_dtr;
                    mydrt_dtr=find(timmd_dtr <= cend(cc) & drend_dtr >= cstrt(cc));
                    if ~isempty(mydrt_dtr)
                        drpertmp_dtr(cc,la,lo)=nanmean(drlen_dtr(mydrt_dtr));  
                    end
                end

            end
        end
    end
end