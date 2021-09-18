% Use changes to mean soil moisture in CLIVAR LENS archives to calcuate the
% "year of new normal" for future projections
% August 2019
% Sam Stevenson
% Updated Jan 2020 to standardize using the pre-industrial control method

% Ensembles to use
ensnames={'cesm_lens','canesm2_lens','gfdl_cm3_lens','csiro_mk36_lens'};
ensncaps={'CESM-LENS','CanESM2','GFDL-CM3','CSIRO-Mk3-6-0'};
strtyr=[0,1850,1860,1850];
ctrlstrtyr=[1,1850,1,1];
enssize=[40,50,20,30];
minyr=1950;
maxyr=2100;
ensdir='/gpfs/fs1/collections/cdg/data/CLIVAR_LE/';
ctrldir='/glade/scratch/samantha/CMIP5/piControl/';

% Drought parameters
windlen=15    % Length of window for running mean
shape=ones(windlen,1)/(windlen);
thr=0.5
if thr < 0
    pos=-1;
else
    pos=1;
end
ynnthr=thr;
refper=[1960,1990];
smbox=[-50,89,0,360];
ensyonn=[];
smdrarr=[];
    smboxarr=[23,38,250,265; ...        % SWUSMEX
            -35,-15,15,30; ...          % SAFR
            -30,-20,115,150; ...        % AUS
            -15,5,285,300; ...          % WAMAZ
            0,30,70,90; ...             % INDIA
            0,20,30,60; ...             % EAFR
            35,50,350,20];              % EUROPE
    smboxarr(smboxarr > 180)=smboxarr(smboxarr > 180)-360;
    regnames={'SWUSMEX','SAFR','AUS','WAMAZ','INDIA','EAFR','EUROPE'};

for ee=1:length(ensnames)
% for ee=1:2
    ensnames{ee}
    
    % Get PI control data
    ttmp=strcat('ls -1',{' '},ctrldir,ensncaps{ee},'/mrso_*');
    [~,list]=system(ttmp{1})
    ctlnames= textscan( list, '%s', 'delimiter', '\n' );
    ctlnames=ctlnames{1};

    pictrlsm=[];
    piyr=[];
    for cc=1:length(ctlnames)
        ctlname=ctlnames{cc}
        nc=netcdf(ctlname);
        PItime=nc{'time'}(:);
        lat=nc{'lat'}(:);
        lon=nc{'lon'}(:);
        [PIyr,PImon,~]=datenumnoleap(PItime,[ctrlstrtyr(ee) 1 1]);
        uyr=unique(PIyr);
        PIyr(1)
        PImon(1)
        PIyr(end)
        PImon(end)

        mylat=find(lat >= smbox(1) & lat <= smbox(2));
        mylon=find(lon >= smbox(3) & lon <= smbox(4));
    
        sm=nc{'mrso'}(:,mylat,mylon);
        sm(abs(sm) > 1e10)=0/0;

        % Annual mean
        piannmn=zeros(length(uyr),size(sm,2),size(sm,3));
        for yy=1:length(uyr)
            myy=find(PIyr == uyr(yy));
            piannmn(yy,:,:)=squeeze(nanmean(sm(myy,:,:),1));
        end

        % Store annual mean soil moisture array
        pictrlsm=cat(1,pictrlsm,piannmn);   
        piyr=cat(1,piyr,uyr');
        
    end
    % Make "detrended" version for later use
    pictrlsm_dtr=pictrlsm-repmat(nanmean(pictrlsm,1),[size(pictrlsm,1) 1 1]);

    % Compute running mean of anomalies and adjust times such that the running
    % mean at location x is the mean over the preceding 'windlen' timesteps
    for la=1:size(pictrlsm,2)
        for lo=1:size(pictrlsm,3)
            smts=movmean(pictrlsm(:,la,lo),windlen,'omitnan');
            smts=circshift(smts,floor(windlen/2));
            pictrlsm(:,la,lo)=smts;             

            smts=runmean(pictrlsm_dtr(:,la,lo),floor(windlen/2));
            smts=circshift(smts,floor(windlen/2));
            pictrlsm_dtr(:,la,lo)=smts;     
        end
    end 

    % Store mean, std. dev. values for standardization
    pimn=nanmean(pictrlsm,1);
    pimn_dtr=nanmean(pictrlsm_dtr,1);
    pistd=nanstd(pictrlsm,1);
    pimn(abs(pimn) < 1e-10)=0/0;
    pimn(abs(pimn_dtr) < 1e-10)=0/0;
    pistd(abs(pistd) < 1e-10)=0/0;
    
%     ttmp=strcat('ls -1',{' '},ensdir,ensnames{ee},'/Lmon/mrso/*rcp85*');
%     [~,list]=system(ttmp{1})
%     runnames= textscan( list, '%s', 'delimiter', '\n' );
%     runnames=runnames{1};
    tmp=dir(strcat(ensdir,ensnames{ee},'/Lmon/mrso/*rcp85*'));
    runnames=cell(length(tmp),1);
    for rr=1:length(runnames)
        runnames{rr}=tmp(rr).name;
    end


    % Coordinates
    smnc=netcdf(strcat(ensdir,ensnames{ee},'/Lmon/mrso/',runnames{1}));
    lat=smnc{'lat'}(:);
    lon=smnc{'lon'}(:);
    mylat=find(lat >= smbox(1) & lat <= smbox(2));
    mylon=find(lon >= smbox(3) & lon <= smbox(4));
    lon(lon > 180)=lon(lon > 180)-360;

    if ee == 1
        % Output variables
        latfix=lat(mylat);
        lonfix=lon(mylon);
    end
    [Xg,Yg]=meshgrid(lon,lat); 
    
    yonnarr=zeros(length(mylat),length(mylon));
    smdrref=zeros(length(runnames),length(mylat),length(mylon));
    sm_refper=[];   % Time series during reference period only
    sm_refperdtr=[];
    
    for rr=1:length(runnames)
        runname=runnames{rr}

        % Read in soil moisture information
        nc=netcdf(strcat(ensdir,ensnames{ee},'/Lmon/mrso/',runname));
        
        if rr == 1
            ctltime=nc{'time'}(:);
            [ctlyr,ctlmon,~]=datenumnoleap(ctltime,[strtyr(ee) 1 1]);
            ctlyr(1)
            ctlmon(1)
            ctlyr(end)
            ctlmon(end)
            myt=find(ctlyr >= minyr & ctlyr <= maxyr);
            ctltime=ctltime(myt);
            [ctlyr,ctlmon,~]=datenumnoleap(ctltime,[strtyr(ee) 1 1]);
            uyr=unique(ctlyr);
            myref=find(uyr >= refper(1) & uyr <= refper(2));       
            uyrfilt=uyr((floor(windlen/2)+1):end);
            smarr=zeros(length(runnames),length(uyrfilt),length(mylat),length(mylon));
        end
        sm=nc{'mrso'}(myt,mylat,mylon);
        sm(abs(sm) > 1e10)=0/0;
        
        % Calculate annual mean
        smann=zeros(length(uyr),length(mylat),length(mylon));
        for yy=1:length(uyr)
            myy=find(ctlyr == uyr(yy));
            smann(yy,:,:)=squeeze(nanmean(sm(myy,:,:),1));
        end
        sm=smann;
                
        % Compute running mean of anomalies and adjust times such that the running
        % mean at location x is the mean over the preceding 'windlen' timesteps
        for la=1:length(mylat)
            for lo=1:length(mylon)
                smts=runmean(sm(:,la,lo),floor(windlen/2));
                smts=circshift(smts,floor(windlen/2));
                sm(:,la,lo)=smts;
            end
        end
        smarr(rr,:,:,:)=sm((floor(windlen/2)+1):end,:,:);

        % Collect ref. period data for standardizing 
        sm_refper=cat(1,sm_refper,squeeze(smarr(rr,myref,:,:)));

    end
    
    % Compute running mean of concatenated reference period data and adjust times such that the running
    % mean at location x is the mean over the preceding 'windlen' timesteps
    for la=1:size(sm_refper,2)
        for lo=1:size(sm_refper,3)                
            smrts=movmean(sm_refper(:,la,lo),windlen,'omitnan');
            smrts=circshift(smrts,floor(windlen/2));
            sm_refper(:,la,lo)=smrts;          
        end
    end
    
%     % Calculate mean, SD over all reference period data
%     smmnref=nanmean(sm_refper,1);
%     smdrref=squeeze(smmnref+ynnthr*pistd);   % drought threshold

    % Standardize data
    smmnref=repmat(permute(nanmean(sm_refper,1),[4 1 2 3]),[size(smarr,1),size(smarr,2),1,1]);
    ctrlsd=repmat(permute(nanstd(sm_refper,1),[4 1 2 3]),[size(smarr,1),size(smarr,2),1,1]);
    smarr=(smarr-smmnref)./ctrlsd;           
        
    % Figure out at what point the forced trend in the 15-year running 
    % mean becomes equal to what is considered a "drought" in the reference
    % period: "year of new normal"
    smarravg=squeeze(nanmean(smarr,1));
    for la=1:length(mylat)
        for lo=1:length(mylon)
            tmp=threshcross(smarravg(:,la,lo),ynnthr,pos);
            if tmp ~= 0
                yonnarr(la,lo)=uyrfilt(tmp);
            end
        end
    end

    lon=lon(mylon);
    lat=lat(mylat);
    yonnarr(yonnarr == 0)=0/0;
    
    % Interpolate epoch averages, differences to common grid
    yonnarr=griddata(lon,lat',yonnarr,lonfix,latfix');
%     smdrref=griddata(lon,lat',smdrref,lonfix,latfix');
%     smmnref=griddata(lon,lat',squeeze(smmnref),lonfix,latfix');
    pistdref=griddata(lon,lat',squeeze(pistd),lonfix,latfix');
    
    lonplot=lonfix;    
    % Reorder array for correct plotting
    myw=find(lonplot < 0);
    mye=find(lonplot >= 0);
    yonnarr=yonnarr(:,[myw' mye']);
%     smdrref=smdrref(:,[myw' mye']);
%     smmnref=smmnref(:,[myw' mye']);

    % Save epoch averages, differences
    ensyonn=cat(1,ensyonn,permute(yonnarr,[3 1 2]));
    size(ensyonn)
        
    if ee == 1
        msk=squeeze(pistdref);
        msk(~isnan(msk))=1;
    end

    % Plot results
    yonnarr(abs(yonnarr) > 1e10)=0/0;
    
    figure(1)
    clf
    m_proj('miller','lon',[0],'lat',[smbox(1) smbox(2)])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lonplot([myw' mye']),latfix,yonnarr.*msk(:,[myw' mye']));
    shading flat
    hold all
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    if ynnthr > 0
        title('Time of Emergence: Pluvial','FontSize',20)
    else
        title('Time of Emergence: Drought','FontSize',20)
    end
    colorbar
    if ynnthr > 0
        strcat('/glade/scratch/samantha/plots/mrsoyrnewnorm_',num2str(ynnthr),'sig_',num2str(windlen),'yrmegapl_',ensnames{ee},'_',num2str(refper(1)),'-',num2str(refper(2)),'ref.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/mrsoyrnewnorm_',num2str(ynnthr),'sig_',num2str(windlen),'yrmegapl_',ensnames{ee},'_',num2str(refper(1)),'-',num2str(refper(2)),'ref.fig'),'fig')
    else
        strcat('/glade/scratch/samantha/plots/mrsoyrnewnorm_',num2str(ynnthr),'sig_',num2str(windlen),'yrmegadr_',ensnames{ee},'_',num2str(refper(1)),'-',num2str(refper(2)),'ref.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/mrsoyrnewnorm_',num2str(ynnthr),'sig_',num2str(windlen),'yrmegadr_',ensnames{ee},'_',num2str(refper(1)),'-',num2str(refper(2)),'ref.fig'),'fig')
    end
end

[Xg,Yg]=meshgrid(lonplot([myw' mye']),latfix); 

% Calculate multi-ensemble mean: require that signal emerges in at least 2
% out of 4 ensembles
emsk=ensyonn;
emsk(~isnan(emsk))=1;
emsk(isnan(emsk))=0;
emsk=squeeze(sum(emsk,1));
emsk(emsk < 2)=0/0;
esum=emsk;
emsk(~isnan(emsk))=1;

ensyonnavg=squeeze(nanmean(ensyonn,1)).*emsk;
ensyonnavg(ensyonnavg >= 2080 | isnan(ensyonnavg) | ensyonnavg == 0)=0/0;

% Save an arry for later use in calculating fraction of points with/without
% emergence
lonfix=lonplot([myw' mye']);
save(strcat('/glade/scratch/samantha/mrsoTOE_ensmndata_',num2str(ynnthr),'sig.mat'),'ensyonnavg','latfix','lonfix')

% Plot ensemble mean
figure(1)
clf
m_proj('miller','lon',[0],'lat',[smbox(1) smbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,ensyonnavg);
shading flat
hold all
m_plot(Xg(esum >= 3),Yg(esum >= 3),'.k','MarkerSize',3)
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
if ynnthr > 0
    title('Time of Emergence: Pluvial','FontSize',20)
else
    title('Time of Emergence: Drought','FontSize',20)
end
colorbar
for bb=1:size(smboxarr,1)
    btmp=smboxarr(bb,:);
    m_line([btmp(3) btmp(4)],[btmp(1) btmp(1)],'Color','k','LineWidth',2)
    m_line([btmp(3) btmp(4)],[btmp(2) btmp(2)],'Color','k','LineWidth',2)
    m_line([btmp(3) btmp(3)],[btmp(1) btmp(2)],'Color','k','LineWidth',2)    
    m_line([btmp(4) btmp(4)],[btmp(1) btmp(2)],'Color','k','LineWidth',2)    
end
if ynnthr > 0
    strcat('/glade/scratch/samantha/plots/mrsoTOE_',num2str(ynnthr),'sig_',num2str(windlen),'yrmegapl_ensmn_',num2str(refper(1)),'-',num2str(refper(2)),'ref_PIctlstd.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/mrsoTOE_',num2str(ynnthr),'sig_',num2str(windlen),'yrmegapl_ensmn_',num2str(refper(1)),'-',num2str(refper(2)),'ref_PIctlstd.fig'),'fig')
else
    strcat('/glade/scratch/samantha/plots/mrsoTOE_',num2str(ynnthr),'sig_',num2str(windlen),'yrmegadr_ensmn_',num2str(refper(1)),'-',num2str(refper(2)),'ref_PIctlstd.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/mrsoTOE_',num2str(ynnthr),'sig_',num2str(windlen),'yrmegadr_ensmn_',num2str(refper(1)),'-',num2str(refper(2)),'ref_PIctlstd.fig'),'fig')
end
