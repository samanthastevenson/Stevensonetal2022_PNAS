% Make maps of surface temperature change for each available large ensemble, both
% with and without background trends
% June 2019
% Sam Stevenson
% Adapting Sept 2021 to make seasonal means

% Ensembles to use
ensnames={'cesm_lens','canesm2_lens','gfdl_cm3_lens','csiro_mk36_lens'};
ensncaps={'CESM-LENS','CanESM2','GFDL-CM3','CSIRO-Mk3-6-0'};
enssize=[40,50,20,30];
strtyr=[0,1850,1860,1850];
minyr=[1920,1850,1860,1850];
maxyr=[2100,2100,2100,2100];
ctrlstrtyr=[1,1850,1,1];


ensdir='/gpfs/fs1/collections/cdg/data/CLIVAR_LE/';
ctrldir='/glade/scratch/samantha/CMIP5/piControl/';

% Parameters
var='mrso'
varcaps='mrso'
mcomp='L';
refper=[1960,1990];
cstrt=[1950,2040];
cend=[2005,2080];
tsbox=[-50,89,0,360]
dostd=1;
maskland=0;
cm=1-[[1*(1:-.1:0)' 1*(1:-.1:0)' 0*(1:-.1:0)'];[0*(0:.1:1)' 1*(0:.1:1)' 1*(0:.1:1)']];  % colorbar

for ee=1:length(ensnames)
    ensnames{ee}
    
    if dostd == 1
    % Get PI control data
    ttmp=strcat('ls -1',{' '},ctrldir,ensncaps{ee},'/',var,'_*');
    [~,list]=system(ttmp{1})
    ctlnames= textscan( list, '%s', 'delimiter', '\n' );
    ctlnames=ctlnames{1};

    pictrlsm_amjjas=[];
    pictrlsm_ondjfm=[];
    piyr=[];
    for cc=1:length(ctlnames)
        ctlname=ctlnames{cc}
        nc=netcdf(ctlname);
        PItime=nc{'time'}(:);
        lat=nc{'lat'}(:);
        lon=nc{'lon'}(:);
        [PIyr,PImon,~]=datenumnoleap(PItime-15,[ctrlstrtyr(ee) 1 1]);
        uyr=unique(PIyr);
        PIyr(1)
        PImon(1)
        PIyr(end)
        PImon(end)

        mylat=find(lat >= tsbox(1) & lat <= tsbox(2));
        mylon=find(lon >= tsbox(3) & lon <= tsbox(4));
        lon(lon > 180)=lon(lon > 180)-360;
    
        if strcmp(var,'pr')
            sm=nc{var}(:,mylat,mylon)*86400;
        else
            sm=nc{var}(:,mylat,mylon);
        end
        sm(abs(sm) > 1e10)=0/0;

        % Annual mean: AMJJAS, ONDJFM
        piannmn_amjjas=zeros(length(uyr),size(sm,2),size(sm,3));
        piannmn_ondjfm=zeros(length(uyr)-1,size(sm,2),size(sm,3));
        for yy=1:length(uyr)
            myy=find(PIyr == uyr(yy) & PImon == 3);
            piannmn_amjjas(yy,:,:)=squeeze(nanmean(sm(myy:myy+2,:,:),1));
            if yy < length(uyr)
                myy=find(PIyr == uyr(yy) & PImon == 9);
                piannmn_ondjfm(yy,:,:)=squeeze(nanmean(sm(myy:myy+2,:,:),1));
            end
        end

        % Store annual mean soil moisture array
        pictrlsm_amjjas=cat(1,pictrlsm_amjjas,piannmn_amjjas);   
        pictrlsm_ondjfm=cat(1,pictrlsm_ondjfm,piannmn_ondjfm);   
        piyr=cat(1,piyr,uyr');
    end
    % Store mean, std. dev. values for standardization
    pimn_amjjas=nanmean(pictrlsm_amjjas,1);
    pistd_amjjas=nanstd(pictrlsm_amjjas,1);    
    pimn_ondjfm=nanmean(pictrlsm_ondjfm,1);
    pistd_ondjfm=nanstd(pictrlsm_ondjfm,1);  
    end
    
    ttmp=strcat('ls -1',{' '},ensdir,ensnames{ee},'/',mcomp,'mon/',var,'/*historical_rcp85*');
    [~,list]=system(ttmp{1})
    runnames= textscan( list, '%s', 'delimiter', '\n' );
    runnames=runnames{1};

    % Coordinates
    tsnc=netcdf(runnames{1});
    lat=tsnc{'lat'}(:);
    lon=tsnc{'lon'}(:);
    min(lon)
    max(lon)
    mylat=find(lat >= tsbox(1) & lat <= tsbox(2));
    mylon=find(lon >= tsbox(3) & lon <= tsbox(4));
    lon(lon > 180)=lon(lon > 180)-360;

    if ee == 1
        % Output variables
        latfix=lat(mylat);
        lonfix=lon(mylon);
        [Xg,Yg]=meshgrid(lonfix,latfix);    
        tscomp20th_amjjas=zeros(length(ensnames),max(enssize),length(mylat),length(mylon));
        tscomp21st_amjjas=zeros(length(ensnames),max(enssize),length(mylat),length(mylon));
        tscomp20th_ondjfm=zeros(length(ensnames),max(enssize),length(mylat),length(mylon));
        tscomp21st_ondjfm=zeros(length(ensnames),max(enssize),length(mylat),length(mylon));
    end
        
    for rr=1:length(runnames)
        runname=runnames{rr}

        % Read in temperature information
        nc=netcdf(runname);
        if strcmp(var,'pr')
            ts=nc{var}(:,mylat,mylon)*86400;
        else
            ts=nc{var}(:,mylat,mylon);
        end
        ts(abs(ts) > 1e10)=0/0;
        
        if rr == 1
            ctltime=nc{'time'}(:);
            [ctlyr,ctlmon,~]=datenumnoleap(ctltime,[strtyr(ee) 1 1]);
            uyr=unique(ctlyr);
            ctlyr(1)
            ctlmon(1)
            ctlyr(end)
            ctlmon(end)
            avgtime_amjjas=zeros(length(uyr),1);
            avgtime_ondjfm=zeros(length(uyr)-1,1);
        
            for yy=1:length(uyr)
                myy=find(ctlyr == uyr(yy) & ctlmon == 3);
                avgtime_amjjas(yy)=nanmean(ctltime(myy:myy+2),1);
                if yy < length(uyr)
                    myy=find(ctlyr == uyr(yy) & ctlmon == 9);
                    avgtime_ondjfm(yy)=nanmean(ctltime(myy:myy+2),1);
                end
            end
            ctltime_amjjas=avgtime_amjjas;
            ctltime_ondjfm=avgtime_ondjfm;
            myt=find(uyr >= minyr(ee) & uyr <= maxyr(ee));
            uyr=uyr(myt);
        
            tsmnarr_amjjas=zeros(length(uyr),length(mylat),length(mylon));
            tsmnarr_ondjfm=zeros(length(uyr)-1,length(mylat),length(mylon));
            tsarr_amjjas=zeros(length(runnames),length(uyr),length(mylat),length(mylon));
            tsarr_ondjfm=zeros(length(runnames),length(uyr)-1,length(mylat),length(mylon));
            tscomparr_amjjas=zeros(length(cstrt),length(runnames),length(mylat),length(mylon));
            tscomparr_ondjfm=zeros(length(cstrt)-1,length(runnames),length(mylat),length(mylon));
        end
        
        % Compute annual mean temperature
        tsmn_amjjas=zeros(length(uyr),size(ts,2),size(ts,3));
        tsmn_ondjfm=zeros(length(uyr)-1,size(ts,2),size(ts,3));

        for yy=1:length(uyr)
            myy=find(ctlyr == uyr(yy) & ctlmon == 3);
            tsmn_amjjas(yy,:,:)=nanmean(ts(myy:myy+2,:,:),1);
            if yy < length(uyr)
                myy=find(ctlyr == uyr(yy) & ctlmon == 9);
                tsmn_ondjfm(yy,:,:)=nanmean(ts(myy:myy+2,:,:),1);
            end
        end
        ts_amjjas=tsmn_amjjas;
        ts_ondjfm=tsmn_ondjfm;
        
        % Restrict time periods
        myt=find(uyr >= minyr(ee) & uyr <= maxyr(ee));
        ts_amjjas=ts_amjjas(myt,:,:);
        ts_ondjfm=ts_ondjfm(myt(1:(end-1)),:,:);
        
        % Save 'raw' annual data to output variables
        tsarr_amjjas(rr,:,:,:)=ts_amjjas;
        tsarr_ondjfm(rr,:,:,:)=ts_ondjfm;
        tsmnarr_amjjas=tsmnarr_amjjas+ts_amjjas./length(runnames);      
        tsmnarr_ondjfm=tsmnarr_ondjfm+ts_ondjfm./length(runnames);      
        
        if dostd == 1
            % Old method, using individual ens. member
%             myref=find(uyr >= refper(1) & uyr <= refper(2));
%             % Calculate mean and standard deviation, normalize
%             tsmn=repmat(nanmean(ts(myref,:,:),1),[size(ts,1),1,1]);
%             tsstd=repmat(nanstd(ts(myref,:,:),1),[size(ts,1),1,1]);
%             ts=(ts-tsmn)./tsstd;   

            % New method ("Method 1") Jan 2020: use PI control
            myref=find(uyr >= refper(1) & uyr <= refper(2));
            tsmn_amjjas=repmat(nanmean(ts_amjjas(myref,:,:),1),[size(ts_amjjas,1),1,1]);
            tsmn_ondjfm=repmat(nanmean(ts_ondjfm(myref,:,:),1),[size(ts_ondjfm,1),1,1]);
            ctrlmn_amjjas=repmat(pimn_amjjas,[length(uyr),1,1]);
            ctrlmn_ondjfm=repmat(pimn_ondjfm,[length(uyr)-1,1,1]);
            ctrlsd_amjjas=repmat(pistd_amjjas,[length(uyr),1,1]);
            ctrlsd_ondjfm=repmat(pistd_ondjfm,[length(uyr)-1,1,1]);
            ts_amjjas=(ts_amjjas-tsmn_amjjas)./ctrlsd_amjjas;             
            ts_ondjfm=(ts_ondjfm-tsmn_ondjfm)./ctrlsd_ondjfm;             
        end
        
        for cc=1:length(cstrt)
            myt=find(uyr >= cstrt(cc) & uyr <= cend(cc));
            tscomparr_amjjas(cc,rr,:,:)=squeeze(nanmean(ts_amjjas(myt,:,:),1));   
            tscomparr_ondjfm(cc,rr,:,:)=squeeze(nanmean(ts_ondjfm(myt(1:(end-1)),:,:),1));   
        end
        
    end
    tscomparr_amjjas(isinf(tscomparr_amjjas))=0/0;
    tscomparr_ondjfm(isinf(tscomparr_ondjfm))=0/0;

    tscomparr_amjjas(tscomparr_amjjas == 0)=0/0;
    tscomparr_ondjfm(tscomparr_ondjfm == 0)=0/0;

    % Composite, compute statistical significance of differences in drought properties 
    % between 20th and 21st centuries
    tscompp_amjjas=zeros(size(ts,2),size(ts,3));
    tscompp_ondjfm=zeros(size(ts,2),size(ts,3));

    % p-values for differences within this individual ensemble
    for la=1:size(ts,2)
        for lo=1:size(ts,3)
            tmp20th_amjjas=squeeze(tscomparr_amjjas(1,:,la,lo));
            tmp21st_amjjas=squeeze(tscomparr_amjjas(2,:,la,lo));
            tmp20th_ondjfm=squeeze(tscomparr_ondjfm(1,:,la,lo));
            tmp21st_ondjfm=squeeze(tscomparr_ondjfm(2,:,la,lo));

            if ~isempty(find(~isnan(tmp20th_amjjas))) && ~isempty(find(~isnan(tmp21st_amjjas)))
                [tscompp_amjjas(la,lo),~]=ranksum(tmp20th_amjjas,tmp21st_amjjas);
            end
            
            if ~isempty(find(~isnan(tmp20th_ondjfm))) && ~isempty(find(~isnan(tmp21st_ondjfm)))
                [tscompp_ondjfm(la,lo),~]=ranksum(tmp20th_ondjfm,tmp21st_ondjfm);
            end
        end
    end
    tscompp_amjjas(tscompp_amjjas < 0.05)=0;  % significant
    tscompp_amjjas(tscompp_amjjas >= 0.05)=1; % insignificant
    tscompp_ondjfm(tscompp_ondjfm < 0.05)=0;  % significant
    tscompp_ondjfm(tscompp_ondjfm >= 0.05)=1; % insignificant

    lon=lon(mylon);
    lat=lat(mylat);
    
    % Interpolate epoch averages, differences to common grid
    tscompp_amjjas=griddata(lon,lat',tscompp_amjjas,lonfix,latfix');
    tscompp_ondjfm=griddata(lon,lat',tscompp_ondjfm,lonfix,latfix');
    
    ts20th_amjjas=zeros(length(runnames),length(latfix),length(lonfix));
    ts20th_ondjfm=zeros(length(runnames),length(latfix),length(lonfix));
    ts21st_amjjas=zeros(length(runnames),length(latfix),length(lonfix));
    ts21st_ondjfm=zeros(length(runnames),length(latfix),length(lonfix));
    
    for rr=1:length(runnames)
        ts20th_amjjas(rr,:,:)=griddata(lon,lat',squeeze(tscomparr_amjjas(1,rr,:,:)),lonfix,latfix');
        ts20th_ondjfm(rr,:,:)=griddata(lon,lat',squeeze(tscomparr_ondjfm(1,rr,:,:)),lonfix,latfix');
        ts21st_amjjas(rr,:,:)=griddata(lon,lat',squeeze(tscomparr_amjjas(2,rr,:,:)),lonfix,latfix');
        ts21st_ondjfm(rr,:,:)=griddata(lon,lat',squeeze(tscomparr_ondjfm(2,rr,:,:)),lonfix,latfix');       
    end
    
    lonplot=lonfix;
    
    % Reorder array for correct plotting
    myw=find(lonplot < 0);
    mye=find(lonplot >= 0);
    ts20th_amjjas=ts20th_amjjas(:,:,[myw' mye']);
    ts20th_ondjfm=ts20th_ondjfm(:,:,[myw' mye']);
    ts21st_amjjas=ts21st_amjjas(:,:,[myw' mye']);
    ts21st_ondjfm=ts21st_ondjfm(:,:,[myw' mye']);
    tscompp_amjjas=tscompp_amjjas(:,[myw' mye']);
    tscompp_ondjfm=tscompp_ondjfm(:,[myw' mye']);
    if ee == 1
        msk=squeeze(tscomparr_ondjfm(1,1,:,:));
        msk=griddata(lon,lat',msk,lonfix,latfix');
        msk(~isnan(msk))=1;
        msk=msk(:,[myw' mye']);
    end    

    [Xg,Yg]=meshgrid(lonplot([myw' mye']),latfix); 
    
    if maskland == 1
        lmsk=double(landmask(Yg,Xg));
        lmsk(lmsk == 0)=0/0;
        ts20th_amjjas=ts20th_amjjas.*repmat(permute(lmsk,[3 1 2]),[size(ts20th_amjjas,1) 1 1]);
        ts20th_ondjfm=ts20th_ondjfm.*repmat(permute(lmsk,[3 1 2]),[size(ts20th_ondjfm,1) 1 1]);
        ts21st_amjjas=ts21st_amjjas.*repmat(permute(lmsk,[3 1 2]),[size(ts21st_amjjas,1) 1 1]);
        ts21st_ondjfm=ts21st_ondjfm.*repmat(permute(lmsk,[3 1 2]),[size(ts21st_ondjfm,1) 1 1]);
    end
    
    % Save epoch averages, differences
    tscomp20th_amjjas(ee,1:length(runnames),:,:)=ts20th_amjjas;
    tscomp20th_ondjfm(ee,1:length(runnames),:,:)=ts20th_ondjfm;
    tscomp21st_amjjas(ee,1:length(runnames),:,:)=ts21st_amjjas;
    tscomp21st_ondjfm(ee,1:length(runnames),:,:)=ts21st_ondjfm;
        
    
    % Plot results
    figure(1)
    clf

    % 20th c.
    subplot('Position',[0.1 0.7 0.4 0.25])
    cla
    m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(ts20th_amjjas,1)));
    shading flat
    hold all
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title('20th c.','FontSize',20)
    colorbar

    % 21st c.
    subplot('Position',[0.1 0.4 0.4 0.25])
    cla
    m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(ts21st_amjjas,1)));
    shading flat
    hold all
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title('21st c.','FontSize',20)
    colorbar

    % Difference
    tmp=(ts21st_amjjas-ts20th_amjjas);
    tmp(abs(tmp) > 1e10)=0/0;
    subplot('Position',[0.1 0.1 0.4 0.25])
    cla
    m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(tmp,1)));
    shading flat
    hold all
    m_plot(Xg(tscompp_amjjas==1),Yg(tscompp_amjjas==1),'.k','MarkerSize',3)
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title('\Delta: 21st c. - 20th c.','FontSize',20)
    colorbar


    % 20th c.
    subplot('Position',[0.55 0.7 0.4 0.25])
    cla
    m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(ts20th_ondjfm,1)));
    shading flat
    hold all
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title('20th c.','FontSize',20)
    colorbar

    % 21st c.
    subplot('Position',[0.55 0.4 0.4 0.25])
    cla
    m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(ts21st_ondjfm,1)));
    shading flat
    hold all
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title('21st c.','FontSize',20)
    colorbar

    % Difference
    subplot('Position',[0.55 0.1 0.4 0.25])
    cla
    m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(ts21st_ondjfm-ts20th_ondjfm,1)));
    shading flat
    hold all
    m_plot(Xg(tscompp_ondjfm==1),Yg(tscompp_ondjfm==1),'.k','MarkerSize',3)
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title('\Delta: 21st c. - 20th c.','FontSize',20)
    colorbar
    strcat('/glade/scratch/samantha/plots/',varcaps,'sestrend_',ensnames{ee},'_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/',varcaps,'sestrend_',ensnames{ee},'_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig'),'fig')
end

tscomp20th_amjjas(abs(tscomp20th_amjjas) > 1e10)=0/0;
tscomp20th_ondjfm(abs(tscomp20th_ondjfm) > 1e10)=0/0;
tscomp21st_amjjas(abs(tscomp21st_amjjas) > 1e10)=0/0;
tscomp21st_ondjfm(abs(tscomp21st_ondjfm) > 1e10)=0/0;
tscomp20th_amjjas(tscomp20th_amjjas == 0)=0/0;
tscomp20th_ondjfm(tscomp20th_ondjfm == 0)=0/0;
tscomp21st_amjjas(tscomp21st_amjjas == 0)=0/0;
tscomp21st_ondjfm(tscomp21st_ondjfm == 0)=0/0;

% OLD METHOD
% Calculate statistical significance of ensemble mean
% follow IPCC-like conventions: 3 or more models agree on sign of change
ensmnp_amjjas=squeeze(nanmean(tscomp21st_amjjas-tscomp20th_amjjas,2));
ensmnp_amjjas(ensmnp_amjjas <= 0)=-1;
ensmnp_amjjas(ensmnp_amjjas > 0)=1;
ensmnp_amjjas=squeeze(nansum(ensmnp_amjjas,1)).*msk;
ensmnp_amjjas(abs(ensmnp_amjjas) < 3)=0;
ensmnp_amjjas(abs(ensmnp_amjjas) > 0)=1;

ensmnp_ondjfm=squeeze(nanmean(tscomp21st_ondjfm-tscomp20th_ondjfm,2));
ensmnp_ondjfm(ensmnp_ondjfm <= 0)=-1;
ensmnp_ondjfm(ensmnp_ondjfm > 0)=1;
ensmnp_ondjfm=squeeze(nansum(ensmnp_ondjfm,1)).*msk;
ensmnp_ondjfm(abs(ensmnp_ondjfm) < 3)=0;
ensmnp_ondjfm(abs(ensmnp_ondjfm) > 0)=1;

[Xg,Yg]=meshgrid(lonplot([myw' mye']),latfix); 

    
% Plot ensemble mean
figure(1)
clf

% 20th c.
subplot('Position',[0.1 0.7 0.4 0.25])
cla
m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(nanmean(tscomp20th_amjjas,2),1)).*msk);
shading flat
hold all
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title('20th c.','FontSize',20)
colorbar

% 21st c.
subplot('Position',[0.1 0.4 0.4 0.25])
cla
m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(nanmean(tscomp21st_amjjas,2),1)).*msk);
shading flat
hold all
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title('21st c.','FontSize',20)
colorbar

% Difference
subplot('Position',[0.1 0.1 0.4 0.25])
cla
m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,(squeeze(nanmean(nanmean(tscomp21st_amjjas,2),1))-squeeze(nanmean(nanmean(tscomp20th_amjjas,2),1))).*msk);
shading flat
hold all
m_coast('color',[0 0 0]);
m_plot(Xg(ensmnp_amjjas==0),Yg(ensmnp_amjjas==0),'.k','MarkerSize',3)
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title('\Delta: 21st c. - 20th c.','FontSize',20)
colorbar


% 20th c.
subplot('Position',[0.55 0.7 0.4 0.25])
cla
m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(nanmean(tscomp20th_ondjfm,2),1)).*msk);
shading flat
hold all
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title('20th c.','FontSize',20)
colorbar

% 21st c.
subplot('Position',[0.55 0.4 0.4 0.25])
cla
m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,squeeze(nanmean(nanmean(tscomp21st_ondjfm,2),1)).*msk);
shading flat
hold all
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title('21st c.','FontSize',20)
colorbar

% Difference
subplot('Position',[0.55 0.1 0.4 0.25])
cla
m_proj('miller','lon',[0],'lat',[tsbox(1) tsbox(2)])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lonplot([myw' mye']),latfix,(squeeze(nanmean(nanmean(tscomp21st_ondjfm,2),1))-squeeze(nanmean(nanmean(tscomp20th_ondjfm,2),1))).*msk);
shading flat
hold all
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
m_plot(Xg(ensmnp_ondjfm==0),Yg(ensmnp_ondjfm==0),'.k','MarkerSize',3)
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title('\Delta: 21st c. - 20th c.','FontSize',20)
colorbar

if dostd == 1
    strcat('/glade/scratch/samantha/plots/',varcaps,'sestrend_std_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/',varcaps,'sestrend_std_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig'),'fig')
else
    strcat('/glade/scratch/samantha/plots/',varcaps,'sestrend_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/',varcaps,'sestrend_',num2str(cstrt(1)),'-',num2str(cend(1)),'vs',num2str(cstrt(2)),'-',num2str(cend(2)),'.fig'),'fig')
end
