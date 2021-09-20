% Calculate changes to the frequency of extreme precipitation in the future
% relative to a reference period (no detrending or standardization)
% November 2020
% Sam Stevenson


global ensnames
global ensncaps
global ensncaps_ctrl
global ctrlstrtyr
global ensdir
global ctrldir
global windlen
global cstrt
global cend
global smbox
global thr
global minyr
global maxyr
global refper
global smboxarr

% Set up environment with global variables
megadrcalcsetup
varname='pr'

% Define other parameters
sgn=1         % sign of exceedance to consider (1 = wet, -1 = dry)

if sgn == -1
    ndysum=1095    % number of days to total for precip (DRY EXTREMES)
    pthr=1         % percentile of distribution to use as threshold (DRY EXTREMES)
else
    ndysum=90      % number of days to total for precip (WET EXTREMES)
    pthr=99        % percentile of distribution to use as threshold (WET EXTREMES)
end

strtyr=[1850,1850,1860,1850];

for ee=1:length(ensnames)
    
    % Load list of ensemble members, coordinate information
    [mylat,mylon,lat,lon,prunnames]=getcoords(ee,varname,'day',strcat('/day/',varname,'/',varname,'_day_',ensncaps{ee},'_historical_rcp85*'));
    [Xg,Yg]=meshgrid(lon,lat); 
    if ee == 1
        latfix=lat(mylat);
        lonfix=lon(mylon);
        lonplot=lonfix;
    end
    
    thr
    
%     
%     prunnames=strrep(runnames,strcat(varname,'_L'),'pr_A');
%     prunnames=strrep(prunnames,'Amon','day');
%     prunnames=strrep(prunnames,'01-','0101-');
%     prunnames=strrep(prunnames,'12.nc','1231.nc');
        
    for rr=1:length(prunnames)
        prunname=prunnames{rr}
        
        if rr==1            
            % Daily-averaged time, coordinate info
            nc=netcdf(strcat(ensdir,ensnames{ee},'/day/pr/',prunname));

            plat=nc{'lat'}(:);
            plon=nc{'lon'}(:);
            myplat=find(plat >= smbox(1) & plat <= smbox(2));
            if smbox(3) == 0 && smbox(4) == 0
                myplon=1:length(plon);
            else
                myplon=find(plon >= smbox(3) & plon <= smbox(4));
            end
            plon(plon > 180)=plon(plon > 180)-360;
            
            dytime=nc{'time'}(:);     % days since strtyr-1-1
            [dyr,dmon,ddy]=datenumnoleap(dytime+15,[strtyr(ee) 1 1]);
            min(dyr)
            max(dyr)
            min(dmon)
            max(dmon)
            mydyrs=find(dyr >= minyr(ee) & dyr <= maxyr(ee));
            dyr=dyr(mydyrs);
            dytime=dytime(mydyrs);
            duyr=unique(dyr);
            
            % Restrict time periods
            duyrfilt=duyr((floor(windlen/2)+1):(end-floor(windlen/2)));
            % Get reference period locations
            myrefd=find(dyr >= refper(1) & dyr <= refper(2)); 
            
            % Output variables
            extrprfreqtmp=NaN(length(prunnames),length(cstrt),length(mylat),length(mylon));    
            pr_refper=zeros(length(prunnames)*length(myrefd),length(mylat),length(mylon));
            pr_refperdtr=zeros(length(prunnames)*length(myrefd),length(mylat),length(mylon));
            
            if ee == 1                
                extrpr20th=zeros(length(ensnames),length(mylat),length(mylon));
                extrpr21st=zeros(length(ensnames),length(mylat),length(mylon));
             end
        end

        
        % Read in daily precip information: reference period only
        nc=netcdf(strcat(ensdir,ensnames{ee},'/day/pr/',prunname));
        pr=nc{'pr'}(myrefd,mylat,mylon)*86400;
    
        % Save precip over the reference period,
        % for use in computing percentiles
        pr_refper((rr-1)*length(myrefd)+1:rr*length(myrefd),:,:)=squeeze(pr); 
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
        
        % Read in daily precip information
        nc=netcdf(strcat(ensdir,ensnames{ee},'/day/pr/',prunname));
        pr=nc{'pr'}(mydyrs,mylat,mylon)*86400;
        
        for la=1:length(mylat)
            for lo=1:length(mylon)                
                prrts=movsum(squeeze(pr(:,la,lo)),ndysum,'omitnan');
                ptmp(:,la,lo)=prrts(1:ndysum:end-floor(ndysum/2));         
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
       
        for cc=1:length(cstrt)
           extrprfreqtmp(rr,cc,:,:)=nansum(thisp(find(dytmp >= cstrt(cc) & dytmp <= cend(cc)),:,:),1)./(cend(cc)-cstrt(cc)+1);
        end    
    end
    
    % Get epoch averages, put on common grid; save epoch averages,
    % differences
    % Interpolate epoch averages, differences to common grid
    dr20th=squeeze(nanmean(extrprfreqtmp(:,1,:,:),1));
    dr20th=griddata(lon(mylon),lat(mylat)',dr20th,lonfix,latfix');
    dr21st=squeeze(nanmean(extrprfreqtmp(:,2,:,:),1));
    dr21st=griddata(lon(mylon),lat(mylat)',dr21st,lonfix,latfix');
    
    % Reorder array for correct plotting
    myw=find(lonfix < 0);
    mye=find(lonfix >= 0);
    dr20th=dr20th(:,[myw' mye']);
    dr21st=dr21st(:,[myw' mye']);

    extrpr20th(ee,:,:)=dr20th;
    extrpr21st(ee,:,:)=dr21st;  
end

% Redo land masking
[Xg,Yg]=meshgrid(lonplot,latfix);
myw=find(lonplot < 0);
mye=find(lonplot >= 0);
Xg=Xg(:,[myw' mye']);
Yg=Yg(:,[myw' mye']);

lmsk=landmask(Yg,Xg);
extrpr20th=extrpr20th.*repmat(permute(lmsk,[3 1 2]),[size(extrpr20th,1) 1 1]);
extrpr21st=extrpr21st.*repmat(permute(lmsk,[3 1 2]),[size(extrpr21st,1) 1 1]);

ensmnp=extrpr21st-extrpr20th;
ensmnp(ensmnp == 0)=0;
ensmnp(ensmnp <= 0)=-1;
ensmnp(ensmnp > 0)=1;
ensmnp=squeeze(nansum(ensmnp,1));
ensmnp(abs(ensmnp) < 3)=0;
ensmnp(abs(ensmnp) >= 3)=1;

% Plot ensemble mean
figure(1)
clf

% 20th c.
subplot('Position',[0.1 0.7 0.8 0.25])
cla
m_proj('robinson','lon',[-180 180],'lat',[smbox(1) smbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(extrpr20th,1)));
shading flat
hold all
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title(strcat('20th c.'),'FontSize',20)
colorbar

% 21st c.
subplot('Position',[0.1 0.4 0.8 0.25])
cla
m_proj('robinson','lon',[-180 180],'lat',[smbox(1) smbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(extrpr21st,1)));
shading flat
hold all
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title(strcat('21st c.'),'FontSize',20)
colorbar

% Difference
subplot('Position',[0.1 0.1 0.8 0.25])
cla
m_proj('robinson','lon',[-180 180],'lat',[smbox(1) smbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,(squeeze(nanmean(extrpr21st,1))-squeeze(nanmean(extrpr20th,1)))./squeeze(nanmean(extrpr20th,1)));
shading flat
hold all
m_coast('color',[0 0 0]);
m_plot(Xg(ensmnp==0),Yg(ensmnp==0),'.k','MarkerSize',3)
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title(strcat('21st c. - 20th c. (%)'),'FontSize',20)
colorbar
for bb=1:size(smboxarr,1)
    btmp=smboxarr(bb,:);
    m_line([btmp(3) btmp(4)],[btmp(1) btmp(1)],'Color','k','LineWidth',2)
    m_line([btmp(3) btmp(4)],[btmp(2) btmp(2)],'Color','k','LineWidth',2)
    m_line([btmp(3) btmp(3)],[btmp(1) btmp(2)],'Color','k','LineWidth',2)    
    m_line([btmp(4) btmp(4)],[btmp(1) btmp(2)],'Color','k','LineWidth',2)    
end

if sgn == 1
    strcat('/glade/scratch/samantha/plots/',varname,'precipextrememap_ensmn_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/',varname,'precipextrememap_ensmn_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig'),'fig')
else
    strcat('/glade/scratch/samantha/plots/',varname,'precipdryextrememap_ensmn_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/',varname,'precipdryextrememap_ensmn_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig'),'fig')
end
