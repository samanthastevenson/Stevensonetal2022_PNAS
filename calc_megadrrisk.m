% Given an input soil moisture time series, calculate the risk of a
% megadrought in that time series
% February 2020
% Sam Stevenson

function [drrisktmp, drrisktmp_dtr]=calc_megadrrisk(uyrfilt,smtmp,smtmp_dtr,mylat,mylon)
    global cstrt
    global cend
    global thr

    drrisktmp=NaN(length(cstrt),length(mylat),length(mylon));
    drrisktmp_dtr=NaN(length(cstrt),length(mylat),length(mylon));
    
    % Calculate risk of drought 
   for la=1:length(mylat)
        for lo=1:length(mylon)
            for cc=1:length(cstrt)
                myt=find(uyrfilt >= cstrt(cc) & uyrfilt <= cend(cc));
                if thr <= 0
                    drs=find(smtmp(myt,la,lo) <= thr);
                else
                    drs=find(smtmp(myt,la,lo) >= thr);                    
                end
                if ~isempty(myt)
                    drrisktmp(cc,la,lo)=length(drs)./length(myt);
                    if length(drs) == 0
                        drrisktmp(cc,la,lo)=0;
                    end
                end

                if thr <= 0
                    drs_dtr=find(smtmp_dtr(myt,la,lo) <= thr);
                else
                    drs_dtr=find(smtmp_dtr(myt,la,lo) >= thr);
                end
                if ~isempty(myt)
                   drrisktmp_dtr(cc,la,lo)=length(drs_dtr)./length(myt);
                end
            end
        end
   end
end