% Redo megadrought property calculations with a code that's easier to
% figure out what the hell it's actually doing
% February 2020
% Sam Stevenson

global ensnames
global ensncaps
global ensncaps_ctrl
global ctrlstrtyr
global ensdir
global ctrldir
global windlen
global cstrt
global smbox
global thr
global smboxarr
global regnames

% Set up environment with global variables
megadrcalcsetup

% Set up output variables
% Risk
drrisk20th=[];
drrisk21st=[];
drrisk_dtr20th=[];
drrisk_dtr21st=[];
% Persistence
drper20th=[];
drper21st=[];
drper_dtr20th=[];
drper_dtr21st=[];
% Frequency
drfreq20th=[];
drfreq21st=[];
drfreq_dtr20th=[];
drfreq_dtr21st=[];
% Severity
drsever20th=[];
drsever21st=[];
drsever_dtr20th=[];
drsever_dtr21st=[];

varname='mrso'
exreg=1;
exrun=1;

for ee=1:length(ensnames)
    % Get mean, SD values from the PI control for use in standardization
    ttmp=strcat('ls -1',{' '},ctrldir,ensncaps_ctrl{ee},'/',varname,'_*');
    [~,list]=system(ttmp{1})
    ctlnames= textscan( list, '%s', 'delimiter', '\n' );
    ctlnames=ctlnames{1};
    [pimn,pimn_dtr,pistd,pistd_dtr]=var2d_standardizeLE_PIctl(ctlnames,smbox,ctrlstrtyr(ee),windlen,varname,1);

    thr
    
    % Load list of ensemble members, coordinate information
    [mylat,mylon,lat,lon,runnames]=getcoords(ee,varname,strcat('/Lmon/',varname,'/',varname,'_Lmon_',ensncaps{ee},'_historical_rcp85*'));
    [Xg,Yg]=meshgrid(lon,lat); 
    if ee == 1
        latfix=lat(mylat);
        lonfix=lon(mylon);
        lonplot=lonfix;
    end
    
    thr
    
    % Get data for all ensemble members
    for rr=1:length(runnames)
        runname=runnames{rr}

        % Read in soil moisture information
        nc=netcdf(strcat(ensdir,ensnames{ee},'/Lmon/',varname,'/',runname));
        sm=nc{varname}(:,mylat,mylon);
        sm(abs(sm) > 1e10)=0/0;
        sm(abs(sm) < 1e-10)=0/0;
        
        if rr==1
            [ctltime,ctlyr,ctlmon,uyrfilt,myref,myyrs]=gettimeinfo(ee,nc);
            uyr=unique(ctlyr);
            % Output variables
            smarr=NaN(length(runnames),length(uyr),length(mylat),length(mylon));
            % Land mask for this ensemble
            indivmsk=squeeze(smarr(1,1,:,:));
            indivmsk(indivmsk == 0)=0/0;
            indivmsk(~isnan(indivmsk))=1;
            
            if ee == 1
                msk=squeeze(sm(1,:,:));
                msk(~isnan(msk))=1;
                myw=find(lonfix < 0);
                mye=find(lonfix >= 0);
                msk=msk(:,[myw' mye']);
            end
        end
        
        % Filter soil moisture to appropriate interval
        sm=sm(myyrs,:,:);

        % Annual average for non-detrended data
        smmn=zeros(length(uyr),size(smarr,3),size(smarr,4));        
        for yy=1:length(uyr)
            myy=find(ctlyr == uyr(yy));
            smmn(yy,:,:)=squeeze(nanmean(sm(myy,:,:),1));
        end
        smtmp=smmn;
        
        % Save annually averaged data to output variables
        smarr(rr,:,:,:)=smtmp;
    end
    
    thr 
    
    % Detrend by removing ensemble mean
    smmnarr=repmat(nanmean(smarr,1),[size(smarr,1) 1 1 1]);
    smarr_dtr=smarr-smmnarr;
    
    % Collect ref. period data for standardizing 
    [sm_refper,sm_refperdtr]=get_refperdata(runnames,myref,smarr,smarr_dtr);
    
    thr
    
    % Loop through ensemble members again, standardize and compute drought
    % properties for detrended and non-detrended values
    drrisktmp=NaN(length(runnames),length(cstrt),length(mylat),length(mylon));
    drrisktmp_dtr=NaN(length(runnames),length(cstrt),length(mylat),length(mylon));
    drpertmp=NaN(length(runnames),length(cstrt),length(mylat),length(mylon));
    drpertmp_dtr=NaN(length(runnames),length(cstrt),length(mylat),length(mylon));
    drsevertmp=NaN(length(runnames),length(cstrt),length(mylat),length(mylon));
    drsevertmp_dtr=NaN(length(runnames),length(cstrt),length(mylat),length(mylon));
    drfreqtmp=NaN(length(runnames),length(cstrt),length(mylat),length(mylon));
    drfreqtmp_dtr=NaN(length(runnames),length(cstrt),length(mylat),length(mylon));

    for rr=1:length(runnames)    
        smtmp=zeros(length(uyr),size(smarr,3),size(smarr,4));        
        smtmp_dtr=zeros(length(uyr),size(smarr,3),size(smarr,4));        
        
        smvartmp=zeros(length(uyr),size(smarr,3),size(smarr,4));        
        smvartmp_dtr=zeros(length(uyr),size(smarr,3),size(smarr,4));   
        
        % Calculate running mean, variance of soil moisture
        for la=1:length(mylat)
            for lo=1:length(mylon)
                smts=movmean(squeeze(smarr(rr,:,la,lo)),windlen,'omitnan');
%                 smts=circshift(smts,floor(windlen/2));
                smtmp(:,la,lo)=smts;

                smts_dtr=movmean(squeeze(smarr_dtr(rr,:,la,lo)),windlen,'omitnan');
%                 smts_dtr=circshift(smts_dtr,floor(windlen/2));
                smtmp_dtr(:,la,lo)=smts_dtr;
                
                smvarts=movvar(smts,windlen,'omitnan');
%                 smvarts=circshift(smvarts,floor(windlen/2));
                smvartmp(:,la,lo)=smvarts;

                smvarts=movvar(smts_dtr,windlen,'omitnan');
%                 smvarts=circshift(smvarts,floor(windlen/2));
                smvartmp_dtr(:,la,lo)=smvarts;
            end
        end
        smtmp=smtmp((floor(windlen/2)+1):(end-floor(windlen/2)),:,:);    % Get rid of end effects
        smtmp_dtr=smtmp_dtr((floor(windlen/2)+1):(end-floor(windlen/2)),:,:);    % Get rid of end effects

        smvartmp=smvartmp((floor(windlen/2)+1):(end-floor(windlen/2)),:,:);    % Get rid of end effects
        smvartmp_dtr=smvartmp_dtr((floor(windlen/2)+1):(end-floor(windlen/2)),:,:);    % Get rid of end effects
        
        % Standardization ("Method 1") Jan 2020: use PI control
        smmn=repmat(nanmean(sm_refper,1),[length(uyrfilt),1,1]);
        smmn_dtr=repmat(nanmean(sm_refperdtr,1),[length(uyrfilt),1,1]);
%         ctrlmn_dtr=repmat(pimn_dtr,[length(uyrfilt),1,1]);
        ctrlsd=repmat(pistd_dtr,[length(uyrfilt),1,1]);
        smtmp=(smtmp-smmn)./ctrlsd;           
        smtmp_dtr=(smtmp_dtr-smmn_dtr)./ctrlsd;           
%         smtmp_dtr=(smtmp_dtr-ctrlmn_dtr)./ctrlsd;        

        % Calculate drought properties!
        % Risk
        [drrisktmp(rr,:,:,:),drrisktmp_dtr(rr,:,:,:)]=calc_megadrrisk(uyrfilt,smtmp,smtmp_dtr,mylat,mylon);
        
        'Risk done'
        thr 
        
        % Persistence
        [drpertmp(rr,:,:,:),drpertmp_dtr(rr,:,:,:)]=calc_megadrper(uyrfilt,smtmp,smtmp_dtr,mylat,mylon);
        
        'persistence done'
        thr
        
        % Severity
        [drsevertmp(rr,:,:,:),drsevertmp_dtr(rr,:,:,:)]=calc_megadrsever(uyrfilt,smtmp,smtmp_dtr,mylat,mylon);
        
        'severity done'
        thr
        
        % Frequency
        [drfreqtmp(rr,:,:,:),drfreqtmp_dtr(rr,:,:,:)]=calc_megadrfreq(uyrfilt,smtmp,smtmp_dtr,mylat,mylon);
        
        'frequency done'
        thr
        
        if rr == exrun
            % Do a test for the specified example region, figure out whether
            % everything makes sense
            xlon=find(lon(mylon) >= smboxarr(exreg,3) & lon(mylon) <= smboxarr(exreg,4));
            xlat=find(lat(mylat) >= smboxarr(exreg,1) & lat(mylat) <= smboxarr(exreg,2));
            smtmp=squeeze(nanmean(nanmean(smtmp(:,xlat,xlon),2),3));
            smtmp_dtr=squeeze(nanmean(nanmean(smtmp_dtr(:,xlat,xlon),2),3));            
            smvartmp=squeeze(nanmean(nanmean(smvartmp(:,xlat,xlon),2),3));
            smvartmp_dtr=squeeze(nanmean(nanmean(smvartmp_dtr(:,xlat,xlon),2),3));            

            % Calculate megadrought properties
            % Risk
            [exdrrisk,exdrrisk_dtr]=calc_megadrrisk(uyrfilt,smtmp,smtmp_dtr,1,1)

            % Persistence
            [exdrper,exdrper_dtr]=calc_megadrper(uyrfilt,smtmp,smtmp_dtr,1,1)

            % Severity
            [exdrsever,exdrsever_dtr]=calc_megadrsever(uyrfilt,smtmp,smtmp_dtr,1,1)

            % Frequency
            [exdrfreq,exdrfreq_dtr]=calc_megadrfreq(uyrfilt,smtmp,smtmp_dtr,1,1)
            
            % Make an example plot
            figure(1)
            clf
            plot(uyrfilt,squeeze(smtmp),'Color','k','LineWidth',2)
            hold all
            plot(uyrfilt,squeeze(smtmp_dtr),'--','Color','r','LineWidth',2)
            strcat('/glade/scratch/samantha/test_',regnames{exreg},'_',ensnames{ee},'_run',num2str(exrun),'_',num2str(windlen),'yrrunmn.fig')
            saveas(gcf,strcat('/glade/scratch/samantha/test_',regnames{exreg},'_',ensnames{ee},'_run',num2str(exrun),'_',num2str(windlen),'yrrunmn.fig'),'fig')

            figure(1)
            clf
            plot(uyrfilt,squeeze(smvartmp),'Color','k','LineWidth',2)
            hold all
            plot(uyrfilt,squeeze(smvartmp_dtr),'--','Color','r','LineWidth',2)
            strcat('/glade/scratch/samantha/test_',regnames{exreg},'_',ensnames{ee},'_run',num2str(exrun),'_',num2str(windlen),'yrrunvar.fig')
            saveas(gcf,strcat('/glade/scratch/samantha/test_',regnames{exreg},'_',ensnames{ee},'_run',num2str(exrun),'_',num2str(windlen),'yrrunvar.fig'),'fig')
            
        end
    
    end
    


    
    % Get epoch averages, put on common grid; save epoch averages,
    % differences
    % Risk
    [dr20th,dr20th_dtr,dr21st,dr21st_dtr]=drepochavg(drrisktmp,drrisktmp_dtr,lon(mylon),lat(mylat),lonfix,latfix);
    drrisk20th=cat(1,drrisk20th,permute(dr20th,[3 1 2]));
    drrisk21st=cat(1,drrisk21st,permute(dr21st,[3 1 2]));
    drrisk_dtr20th=cat(1,drrisk_dtr20th,permute(dr20th_dtr,[3 1 2]));
    drrisk_dtr21st=cat(1,drrisk_dtr21st,permute(dr21st_dtr,[3 1 2]));
    % Persistence
    [dr20th,dr20th_dtr,dr21st,dr21st_dtr]=drepochavg(drpertmp,drpertmp_dtr,lon(mylon),lat(mylat),lonfix,latfix);
    drper20th=cat(1,drper20th,permute(dr20th,[3 1 2]));
    drper21st=cat(1,drper21st,permute(dr21st,[3 1 2]));
    drper_dtr20th=cat(1,drper_dtr20th,permute(dr20th_dtr,[3 1 2]));
    drper_dtr21st=cat(1,drper_dtr21st,permute(dr21st_dtr,[3 1 2]));
    % Severity
    [dr20th,dr20th_dtr,dr21st,dr21st_dtr]=drepochavg(drsevertmp,drsevertmp_dtr,lon(mylon),lat(mylat),lonfix,latfix);
    drsever20th=cat(1,drsever20th,permute(dr20th,[3 1 2]));
    drsever21st=cat(1,drsever21st,permute(dr21st,[3 1 2]));
    drsever_dtr20th=cat(1,drsever_dtr20th,permute(dr20th_dtr,[3 1 2]));
    drsever_dtr21st=cat(1,drsever_dtr21st,permute(dr21st_dtr,[3 1 2]));    
    % Frequency
    [dr20th,dr20th_dtr,dr21st,dr21st_dtr]=drepochavg(drfreqtmp,drfreqtmp_dtr,lon(mylon),lat(mylat),lonfix,latfix);
    drfreq20th=cat(1,drfreq20th,permute(dr20th,[3 1 2]));
    drfreq21st=cat(1,drfreq21st,permute(dr21st,[3 1 2]));
    drfreq_dtr20th=cat(1,drfreq_dtr20th,permute(dr20th_dtr,[3 1 2]));
    drfreq_dtr21st=cat(1,drfreq_dtr21st,permute(dr21st_dtr,[3 1 2]));
    
    'drepochavg done'
    thr
end

% Plot composite results
plot_megadrpropmap('risk',drrisk20th,drrisk21st,drrisk_dtr20th,drrisk_dtr21st,lonplot,latfix,varname);
plot_megadrpropmap('freq',drfreq20th,drfreq21st,drfreq_dtr20th,drfreq_dtr21st,lonplot,latfix,varname);
plot_megadrpropmap('per',drper20th,drper21st,drper_dtr20th,drper_dtr21st,lonplot,latfix,varname);
plot_megadrpropmap('sever',drsever20th,drsever21st,drsever_dtr20th,drsever_dtr21st,lonplot,latfix,varname);

thr