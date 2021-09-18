% Calculate the time at which a time series falls below a given value and
% does not recover again
% August 2019
% Sam Stevenson
% ts = time series of interest
% thr = threshold
% pos = type of exceedance to look for (positive or negative)
% ynn = position within ts where the threshold is crossed 'for good'

function ynn=threshcross(ts,thr,pos)
    
    bin_ts=ts;
    bin_ts(ts > thr)=1;
    bin_ts(ts <= thr)=0;
    bindiff=diff(bin_ts,1);
    ynn=[];
    
    ynn_p = find(bindiff == 1,1,'last')+1;
    ynn_n = find(bindiff == -1,1,'last');
    
    if ~isempty(ynn_p) && ~isempty(ynn_n)
        if pos > 0 && ynn_p > ynn_n
            ynn = find(bindiff == 1,1,'last')+1;
        end
        if pos < 0 && ynn_n > ynn_p
            ynn = find(bindiff == -1,1,'last');
        end
    else
        if pos > 0 
            ynn = ynn_p;
        end
        if pos < 0 
            ynn = ynn_n;
        end        
    end
    
end