% Code to compute the day of the year assuming no leap years
% Created April 19, 2012
% Samantha Stevenson

function [yr,mon,day,hr,min,s]=datenumnoleap(tvec,tstart)

    hr=zeros(1,length(tvec));
    min=zeros(1,length(tvec));
    s=zeros(1,length(tvec));
    mon=zeros(1,length(tvec));
    day=zeros(1,length(tvec));
    yr=zeros(1,length(tvec));
    
    daysofyr=[31,28,31,30,31,30,31,31,30,31,30,31];
    cumday=cumsum(daysofyr);
    
    yrstart=tstart(1);
    monstart=tstart(2);
    
    % Find number of days into the year at starting point
    if monstart > 1
        daystart=tstart(3) + cumday(monstart-1);
    else
        daystart=tstart(3);
    end
           
    for tt=1:length(tvec)
        yr(tt)=floor(tvec(tt)/sum(daysofyr)) + yrstart;        
%         daytmp=tvec(tt)-(yr(tt)-yrstart)*sum(daysofyr)+1;
        daytmp=tvec(tt)-((yr(tt)-yrstart)*365)+1;

        ttmp=find(daytmp > cumday, 1, 'last');
        if ~isempty(ttmp)
            mon(tt)=ttmp+1;
        else
            mon(tt)=1;
        end
        
        if isempty(mon(tt)) mon(tt) = 1; end

        if mon(tt) > 1
            day(tt)=daytmp-cumday(mon(tt)-1);
        else
            day(tt)=daytmp;
        end        
    end
end