% Get time-related information for a given ensemble member
% Feb 2020
% Sam Stevenson

function [ctltime,ctlyr,ctlmon,uyrfilt,myref,myyrs]=gettimeinfo(ee,nc)
    global strtyr
    global minyr
    global maxyr
    global windlen
    global refper
    
    ctltime=nc{'time'}(:);
    [ctlyr,ctlmon,~]=datenumnoleap(ctltime,[strtyr(ee) 1 1]);

    myyrs=find(ctlyr >= minyr(ee) & ctlyr <= maxyr(ee));
    ctlyr=ctlyr(myyrs);
    ctltime=ctltime(myyrs);
    uyr=unique(ctlyr);
    ctlyr(1)
    ctlmon(1)
    ctlyr(end)
    ctlmon(end) 

    % Restrict time periods
    uyrfilt=uyr((floor(windlen/2)+1):(end-floor(windlen/2)));
    % Get reference period locations
    myref=find(uyr >= refper(1) & uyr <= refper(2)); 
end