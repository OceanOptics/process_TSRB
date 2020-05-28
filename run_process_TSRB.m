function TSRB = run_process_TSRB(filedir,tiltlimit)
%
% Code to process radiometry data from the Satlantic HyperPro in buoy mode,
% a.k.a. TSRB (Tethered Spectral Radiometer Buoy)
%
% This function requires *.dat files that have been processed from *.raw files
% using SatCon; comma delimited, apply time tages, and apply appropriate
% immersion coefficients. Name processed data by appending the cal file name
% to the original file name.
%
% For each *.raw file, five *.dat files should be
% generated. The order of the five files will not matter, but they must be
% appear in groups of five in the folder that contains the *.dat files.
%
%   Inputs:
%
%       filedir = directory of *.dat files
%       tiltlimit (optional; default = 5) = tilt value in degrees above which all data will be removed
%
%   Output:
%
%       TSRB  = structure of processed radiometry data, or an array of
%               structures for more than one .raw sample file. Structure contents:
%
%               TSRB.Ed = downwelling irradiance spectra
%               TSRB.Lu = upwelling radiance spectra
%               TSRB.Lw = water-leaving upwelling radiance spectra
%               TSRB.Rrs = Remote-sensing reflectance spectra
%               TSRB.wl = wavelength
%               TSRB.datetime = date and time of each spectrum in matlab datenum format
%               TSRB.ChlOC = chlorophyll concentration estimated using OC4 band ratio algorithm
%               TSRB.units = 'Ed: uW cm^-2 nm^-1   Lu & Lw: uW cm^-2 nm^-1 sr^-1   Rrs: sr^-1   wl: nm  ChlOC: mg m^-3';
%
%
%  Naming convention from cal files/SatCon:
%
%       PED - Ed dark file (or HED for Es dark file)
%       HPE - Ed light file (or Hse for Es light file)
%       PLD - Lu dark file
%       HPL - Lu light file
%       MPR - Tilt etc. file
%
%  The following functions are called:
%
%       read_datfiles.m
%       read_tiltfile.m
%       calculate_Rrs.m
%       get_ap_oc_from_Rrs.m
%       get_water_iops.m
%
%  The following constants are used in calculate_Rrs.m:
%
%   n_wat = 1.34;   % refractive index of seawater
%   t = 0.98;       % radiance transmittance of the surface (Mobley, 1994)
%   z = 0.20;       % meters below surface of Lu sensor
%   mu_u = 0.5;     % average upwelling cosine
%
% A. Chase 2017, University of Maine, Orono ME USA
% alison.p.chase@maine.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    tiltlimit = 5;
end

files = dir([filedir,'*.dat']);

% Loop through number of samples (number of .dat files/5)
n = 0;
for imeas = 1:length(files)/5
    
    disp(['Looping through ',num2str(length(files)/5),' samples'])
    
    % Read in and dark-correct Ed data
    for i=1:5
        currfile = files(n+i).name;
        if contains(currfile,'PED') == 1 || contains(currfile,'HED') == 1
            Edark = currfile; disp(['Ed dark file: ',Edark])
        end
        if contains(currfile,'HPE') == 1 || contains(currfile,'Hse') == 1
            Elight = currfile; disp(['Ed light file: ',Elight])
        end
    end
    
    [EDdata,wlEd,tt2,SN,INTtt2,filetype] = read_datfiles(Edark,Elight,filedir);
    
    % Read in and dark-correct Lu data
    for i=1:5
        currfile = files(n+i).name;
        if contains(currfile,'PLD') == 1
            Ldark = currfile;disp(['Lu dark file: ',Ldark])            
        end
        if contains(currfile,'HPL') == 1
            Llight = currfile;disp(['Lu light file: ',Llight])
        end
    end
    [LUdata,wlLu,tt22,SN2,INTtt22,filetype2] = read_datfiles(Ldark,Llight,filedir);
    
    % Read in the tilt etc. file
    for i=1:5
        currfile = files(n+i).name;
        if contains(currfile,'MPR') == 1
            tiltfile = currfile;disp(['Tilt file: ',tiltfile])
        end
    end
    [tiltx,tilty,tilt_tt2] = read_tiltfile(tiltfile,filedir,tiltlimit);
    
    % interpolate the larger array (either Ed or Lu) to the smaller array
    % (in time; go to lowest time resolution)
    if max(size(tt2))>max(size(tt22))
        EDdata=interp1(tt2,EDdata,tt22);
        tt2=tt22;
    else
        LUdata=interp1(tt22,LUdata,tt2);
        tt22=tt2;
    end
        
    %find match-ups between the radiometric times and the "good" tilt times
    tt2_keep=[];
    tilt_keep=[];
    
    for itilt=1:max(size(tt2))
        match=find(tt2(itilt)==tilt_tt2);
        if isempty(match)==0
            tt2_keep=[tt2_keep,itilt];
            tilt_keep=[tilt_keep,match(1)];
        end
    end
    
    Ed=EDdata(tt2_keep,:);
    Lu=LUdata(tt2_keep,:);
    time_keep = tt2(tt2_keep);
    
    % interpolate the Lu to the Ed (in spectral space)
    for iii=1:max(size(Lu(:,1)))
        Lu(iii,:)=interp1(wlLu,Lu(iii,:),wlEd,'linear','extrap');
    end
    
    % Only keep the Ed and Lu data that falls between the 25th and 75th
    % percentile of Ed data. Not repeated for Lu spectra since there is less
    % variation and keeping Lu IQR spectra results in removal of most spectra
    % due to low variability in the red.
    Ed_prc = prctile(Ed,[25,50,75]);
    Ed_keep=[];
    for jjj=1:max(size(Ed(:,1)))
        test_Ed_high=Ed(jjj,:)>Ed_prc(3,:);
        test_Ed_low=Ed(jjj,:)<Ed_prc(1,:);
        if max(test_Ed_high)+max(test_Ed_low)>0
            Ed_keep=[Ed_keep,jjj];
        end
    end
    full_jjj=1:jjj;
    good_Ed=setxor(full_jjj,Ed_keep);
    
    time_mid = time_keep(good_Ed);
    Ed_mid = Ed(good_Ed,:);
    Lu_mid = Lu(good_Ed,:); % select the Lu data that matches the middle 50th percentile of Ed data
    
    % Calculate the Ed and Lu medians and std devs
    for istd=1:max(size(Ed_mid(1,:)))
        Ed_std(istd) = nanstd(Ed_mid(:,istd));
        Lu_std(istd) = nanstd(Lu_mid(:,istd));
        Ed_avg(istd) = nanmean(Ed_mid(:,istd));
        Lu_avg(istd) = nanmean(Lu_mid(:,istd));
    end
    
    [Rrs,Lw,ChlOC] = calculate_Rrs(Ed_mid,Lu_mid,wlEd);
    
    % store a structure with all Ed, Lu, Lw, Rrs values
    TSRB(imeas).Ed = Ed_mid;
    TSRB(imeas).Lu = Lu_mid;
    TSRB(imeas).Lw = Lw;
    TSRB(imeas).Rrs = Rrs';
    TSRB(imeas).wl = wlEd;
    TSRB(imeas).datetime = time_mid;
    TSRB(imeas).ChlOC = ChlOC;
    TSRB(imeas).units = 'Ed: uW cm^-2 nm^-1   Lu & Lw: uW cm^-2 nm^-1 sr^-1   Rrs: sr^-1   wl: nm  ChlOC: mg m^-3';
    
    n = n+5;
    
    plot(wlEd,Rrs')
    pause
    
end
