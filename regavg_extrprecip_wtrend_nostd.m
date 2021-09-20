% Make a time series of running occurrence of precipitation extremes over a
% specified region, in the various multiple large ensembles
% December 2020
% Sam Stevenson


global ensnames
global ensncaps
global ensdir
global windlen
global thr
global minyr
global maxyr
global refper
global smboxarr
global regnames

% Set up environment with global variables
megadrcalcsetup
varname='pr'
% rmpath /glade/u/home/samantha/codes/toolbox/matutil/

% Define other parameters
sgn=-1         % sign of exceedance to consider (1 = wet, -1 = dry)

if sgn == -1
    ndysum=1095    % number of days to total for precip (DRY EXTREMES)
    pthr=1         % percentile of distribution to use as threshold (DRY EXTREMES)
else
    ndysum=90      % number of days to total for precip (WET EXTREMES)
    pthr=99        % percentile of distribution to use as threshold (WET EXTREMES)
end

nyr=30;         % number of years for running window
rnum=7;
smbox=smboxarr(rnum,:);
strtyr=[1850,1850,1860,1850];

for ee=1:length(ensnames)
    
    % Load list of ensemble members, coordinate information
    [~,~,lat,lon,prunnames]=getcoords(ee,varname,'day',strcat('/day/',varname,'/',varname,'_day_',ensncaps{ee},'_historical_rcp85*'));
    
    thr

    for rr=1:length(prunnames)
        prunname=prunnames{rr}
        
        if rr==1            
            % Daily-averaged time, coordinate info
            nc=netcdf(strcat(ensdir,ensnames{ee},'/day/pr/',prunname));

            plat=nc{'lat'}(:);
            plon=nc{'lon'}(:);
            plon(plon > 180)=plon(plon > 180)-360;
            myplat=find(plat >= smbox(1) & plat <= smbox(2));
            if smbox(3) == 0 && smbox(4) == 0
                myplon=1:length(plon);
            else
                if (smbox(3) >= 0 && smbox(4) >= 0) || (smbox(3) <= 0 && smbox(4) <= 0)
                    myplon=find(plon >= smbox(3) & plon <= smbox(4));
                    myplone=[];
                else
                    myplonw=find(plon >= smbox(3) & plon < 0);
                    myplone=find(plon >= 0 & plon <= smbox(4));
                    myplon=[myplonw' myplone'];
                end                
            end
            
            dytime=nc{'time'}(:);     % days since strtyr-1-1
            [dyr,dmon,ddy]=datenumnoleap(dytime,[strtyr(ee) 1 1]);
            min(dyr)
            max(dyr)
            min(dmon)
            max(dmon)
            mydyrs=find(dyr >= minyr(ee) & dyr <= maxyr(ee));
            dytime=datenum([dyr(mydyrs)' dmon(mydyrs)' ddy(mydyrs)']);
            dyr=dyr(mydyrs);
            duyr=unique(dyr);
            
            % Restrict time periods
            duyrfilt=duyr((floor(nyr/2)+1):(end-floor(nyr/2)));
            % Get reference period locations
            myrefd=find(dyr >= refper(1) & dyr <= refper(2)); 
            
            % Output variables
            prarr=zeros(length(prunnames),length(mydyrs),length(myplat),length(myplon));
            extrprfreqtmp=NaN(length(prunnames),length(dytime(floor(ndysum/2):ndysum:end-floor(ndysum/2))));    
            pr_refper=zeros(length(prunnames)*length(myrefd),length(myplat),length(myplon));
            pr_refperdtr=zeros(length(prunnames)*length(myrefd),length(myplat),length(myplon));
            
            if ee == 1     
                extrprfreq=cell(length(ensnames),1);
                timearr=cell(length(ensnames),1);
             end
        end

        
        % Read in daily precip information: reference period only
        nc=netcdf(strcat(ensdir,ensnames{ee},'/day/pr/',prunname));
        if (smbox(3) >= 0 && smbox(4) >= 0) || (smbox(3) <= 0 && smbox(4) <= 0)
            pr=nc{'pr'}(mydyrs,myplat,myplon)*86400;
        else
            precw=nc{'pr'}(mydyrs,myplat,myplonw);
            prece=nc{'pr'}(mydyrs,myplat,myplone);
            pr=cat(3,precw,prece)*86400;
        end
            
        % Save precip over the reference period,
        % for use in computing percentiles
        pr_refper((rr-1)*length(myrefd)+1:rr*length(myrefd),:,:)=squeeze(pr(myrefd,:,:)); 
        prarr(rr,:,:,:)=pr;
    end
    
    'Start running sum of precipitation'
    % Compute running mean of concatenated reference period data and adjust times such that the running
    % mean at location x is the mean over the preceding 'windlen' timesteps
    for la=1:size(pr_refper,2)
        for lo=1:size(pr_refper,3)                
            prrts=movsum(squeeze(pr_refper(:,la,lo)),ndysum,'omitnan');
            pr_refper(:,la,lo)=prrts;         
        end
    end     
    pr_refper=pr_refper((floor(ndysum/2)+1):ndysum:(end-floor(ndysum/2)),:,:);        % Get rid of end effects, overlapping periods
    
    
    % Figure out precipitation thresholds based on some percentile over
    % the reference period
    prthr=prctile(pr_refper,pthr,1);
        
    % Loop through ensemble members again, standardize and compute drought
    % properties for detrended and non-detrended values
    'Begin loop 2: compute occurrence of extremes'
 
    for rr=1:length(prunnames)    
        prunname=prunnames{rr}
        ptmp=zeros(length(dytime(floor(ndysum/2):ndysum:end-floor(ndysum/2))),length(myplat),length(myplon));
        dytmp=dyr(floor(ndysum/2):ndysum:end-floor(ndysum/2));
        
        for la=1:length(myplat)
            for lo=1:length(myplon)                
                prrts=movsum(squeeze(prarr(rr,:,la,lo)),ndysum,'omitnan');
                ptmp(:,la,lo)=prrts(floor(ndysum/2):ndysum:end-floor(ndysum/2));         
            end
        end
    
       % Calculate frequencies of extreme precipitation overall
       thisp=ptmp-prthr;
       if sgn == 1
           % Wet periods: find all places where precip is greater than
           % threshold
           thisp(thisp < 0)=0;
           thisp(thisp > 0)=1;
       else
           % Dry periods: find all places where precip is less than
           % threshold
           thisp(thisp > 0)=0;
           thisp(thisp < 0)=1;
       end
       
       nt=nyr*365./ndysum;  % number of 'ndysum' periods in the window
       extrtmp=zeros(size(thisp));
       for la=1:length(myplat)
           for lo=1:length(myplon) 
               extrtmp(:,la,lo)=movsum(thisp(:,la,lo),floor(nt),'omitnan');
           end
       end
       extrprfreqtmp(rr,:)=squeeze(nanmean(nanmean(extrtmp,2),3));
    end
    
    extrprfreq{ee}=extrprfreqtmp;
    timearr{ee}=dytime(floor(ndysum/2):ndysum:end-floor(ndysum/2));
end

% Plot results
figure(1)
clf
hens=[];
cols=[0 0 1;
    0 0 0;
    1 0 0;
    0 1 0;
    0.5 0.5 0.5;
    ];

for ee=1:length(ensnames)
    etmp=extrprfreq{ee};
    ensmn=squeeze(prctile(etmp,50,1));
    ensmax=squeeze(prctile(etmp,75,1));
    ensmin=squeeze(prctile(etmp,25,1));
    ttmp=timearr{ee};
    [yr,mon,dy]=datevec(ttmp);
    myt=find(yr >= min(yr)+floor(nyr/2) & yr <= max(yr)-floor(nyr/2));
    
    h=area(ttmp(myt),[ensmin(myt)' ensmax(myt)'-ensmin(myt)']);
    set(h(2),'FaceColor',cols(ee,:),'FaceAlpha',0.5)
    set(h(1),'FaceColor','none')
    hold all
    hens=[hens,plot(ttmp(myt),ensmn(myt),'Color',cols(ee,:),'LineWidth',2)];
end
set(gca,'FontSize',24)
datetick
legend(hens,ensnames(1:end))
title(strcat(regnames{rnum},': ',num2str(pthr),'% percentile precip occurrence'))
ylabel(strcat('#/',num2str(nyr),' yrs'))

if sgn == 1
    strcat('/glade/scratch/samantha/plots/precipextreme_',regnames{rnum},'_refyrs',num2str(refper(1)),'-',num2str(refper(2)),'_',num2str(pthr),'pctile_ndysum',num2str(ndysum),'_windlen',num2str(nyr),'yrs.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/precipextreme_',regnames{rnum},'_refyrs',num2str(refper(1)),'-',num2str(refper(2)),'_',num2str(pthr),'pctile_ndysum',num2str(ndysum),'_windlen',num2str(nyr),'yrs.fig'),'fig')
else
    strcat('/glade/scratch/samantha/plots/precipdryextreme_',regnames{rnum},'_refyrs',num2str(refper(1)),'-',num2str(refper(2)),'_',num2str(pthr),'pctile_ndysum',num2str(ndysum),'_windlen',num2str(nyr),'yrs.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/precipdryextreme_',regnames{rnum},'_refyrs',num2str(refper(1)),'-',num2str(refper(2)),'_',num2str(pthr),'pctile_ndysum',num2str(ndysum),'_windlen',num2str(nyr),'yrs.fig'),'fig')
end
