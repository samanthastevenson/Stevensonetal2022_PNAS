% Given an input time series, compute the frequency of megadroughts in that
% time series
% February 2020
% Sam Stevenson

function [drfreqtmp,drfreqtmp_dtr]=calc_megadrfreq(uyrfilt,smtmp,smtmp_dtr,mylat,mylon)
    global cstrt
    global cend
    global thr
    
    drfreqtmp=NaN(length(cstrt),length(mylat),length(mylon));
    drfreqtmp_dtr=NaN(length(cstrt),length(mylat),length(mylon));
    
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
                        drfreqtmp(cc,la,lo)=100*length(mydrt)./(cend(cc)-cstrt(cc)+1);   % droughts per century
                    else
                        drfreqtmp(cc,la,lo)=0;
                    end
                else
                    drfreqtmp(cc,la,lo)=0;
                end
                if min(drloc_dtr) > 0
                    timmd_dtr=uyrfilt(drloc_dtr);
                    drend_dtr=timmd_dtr+drlen_dtr;
                    mydrt_dtr=find(timmd_dtr <= cend(cc) & drend_dtr >= cstrt(cc));
                    if ~isempty(mydrt_dtr)
                        drfreqtmp_dtr(cc,la,lo)=100*length(mydrt_dtr)./(cend(cc)-cstrt(cc)+1);   % droughts per century
                    else
                        drfreqtmp_dtr(cc,la,lo)=0;
                    end
                else
                    drfreqtmp_dtr(cc,la,lo)=0;
                end

            end
        end
   end
end