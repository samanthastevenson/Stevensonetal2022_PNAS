% Load parameters for megadrought calculations
% February 2020
% Sam Stevenson

function megadrcalcsetup
    global ensnames
    global ensncaps
    global ensncaps_ctrl
    global strtyr
    global minyr
    global maxyr
    global ctrlstrtyr
    global enssize
    global ensdir
    global ctrldir
    global windlen
    global thr
    global refper
    global cstrt
    global cend
    global smbox
    global smboxarr
    global regnames
    
    % Ensembles to use
    ensnames={'cesm_lens','canesm2_lens','gfdl_cm3_lens','csiro_mk36_lens'};
    ensncaps={'CESM1-CAM5','CanESM2','GFDL-CM3','CSIRO-Mk3-6-0'};
    ensncaps_ctrl={'CESM-LENS','CanESM2','GFDL-CM3','CSIRO-Mk3-6-0'};
    strtyr=[0,1850,1860,1850];
    minyr=[1950,1950,1950,1950];
    maxyr=[2100,2100,2100,2100];
    ctrlstrtyr=[1,1850,1,1];
    enssize=[40,50,20,30];
    
    ensdir='/gpfs/fs1/collections/cdg/data/CLIVAR_LE/';
    ctrldir='/glade/scratch/samantha/CMIP5/piControl/';

    % Drought parameters
    windlen=15    % Length of window for running mean
    thr=0.5
    refper=[1960,1990];
    cstrt=[1950,2040];
    cend=[2005,2080];
    smbox=[-50,89,0,360];
    smbox(smbox > 180)=smbox(smbox > 180)-360;
    smboxarr=[23,38,250,265; ...        % SWUSMEX
            -35,-15,15,30; ...          % SAFR
            -30,-20,115,150; ...        % AUS
            -15,5,285,300; ...          % WAMAZ
            0,30,70,90; ...             % INDIA
            0,20,30,60; ...             % EAFR
            35,50,350,20];              % EUROPE
    smboxarr(smboxarr > 180)=smboxarr(smboxarr > 180)-360;
    regnames={'SWUSMEX','SAFR','AUS','WAMAZ','INDIA','EAFR','EUROPE'};
end