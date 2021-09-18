% Given an input time series, calculate the severity of megadroughts in
% that time series
% February 2020
% Sam Stevenson

function [drsevertmp,drsevertmp_dtr]=calc_megadrsever(uyrfilt,smtmp,smtmp_dtr,mylat,mylon)
    global cstrt
    global cend
    global thr
    
    drsevertmp=NaN(length(cstrt),length(mylat),length(mylon));
    drsevertmp_dtr=NaN(length(cstrt),length(mylat),length(mylon));
    
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
                    % Store drought and pluvial properties for 20th, 21st centuries
                    myt=find(timmd <= cend(cc) & drend >= cstrt(cc));
                    if ~isempty(myt)
                        drtmp=[];
                        for dd=1:length(myt)
                            drtmp=cat(1,drtmp,nanmean(smtmp(drloc(myt(dd)):(drloc(myt(dd))+drlen(myt(dd))-1),la,lo)));    % average soil moisture anomaly during drought
                        end
                        drsevertmp(cc,la,lo)=nanmean(drtmp);  
                    end
                end
                if min(drloc_dtr) > 0
                    timmd_dtr=uyrfilt(drloc_dtr);
                    drend_dtr=timmd_dtr+drlen_dtr;
                    myt=find(timmd_dtr <= cend(cc) & drend_dtr >= cstrt(cc));
                    if ~isempty(myt)
                        drtmp=[];
                        for dd=1:length(myt)
                            drtmp=cat(1,drtmp,nanmean(smtmp_dtr(drloc_dtr(myt(dd)):(drloc_dtr(myt(dd))+drlen_dtr(myt(dd))-1),la,lo)));    % average soil moisture anomaly during drought
                        end
                        drsevertmp_dtr(cc,la,lo)=nanmean(drtmp);  
                    end
                end
            end
        end
   end
end