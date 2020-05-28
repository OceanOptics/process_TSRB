function [Rrs,Lw,chlOCnew] = calculate_Rrs(Ed_mid,Lu_mid,wl)

% To calculate Rrs, we need water-leaving radiance (Lw).
% Extrapolate Lu to Lu(0-) and then get Lw using:
%   K_Lu ~(a_w+a_p)/0.5, <mu_u>~0.5
%   Lu(lambda, 0-) = Lu(lambda, z)*exp(-KLu*z)
%   Lw(lambda)=t* Lu(0-,lambda)/(n^2)

mu_u = 0.5;     % average upwelling cosine
z = 0.20;       % meters below surface of Lu sensor
n_wat = 1.34;   % refractive index of seawater
t = 0.98;       % radiance transmittance of the surface (Mobley, 1994)

% For a_p in the above equation, use an a_p
% spectrum calcuated from an estimated Chl value.
% Estimate Chl from Rrs using band ratio, then use A and B coeffs from
% Chase et al. 2017 JGR paper to get an a_p spectrum, then iterate if necessary.

% Calculate the Ed and Lu medians and std devs
for istd=1:max(size(Ed_mid(1,:)))
    Ed_std(istd) = nanstd(Ed_mid(:,istd));
    Lu_std(istd) = nanstd(Lu_mid(:,istd));
    Ed_avg(istd) = nanmean(Ed_mid(:,istd));
    Lu_avg(istd) = nanmean(Lu_mid(:,istd));
end

% First guess at Rrs for input to function to get a_p
Rrs_guess = Lu_avg./Ed_avg;

[ap,apwl,chlOC] = get_ap_oc_from_Rrs(Rrs_guess,wl);

% get a_w and bb_w values
temp = NaN; sal = NaN;
[a_sw,bb_sw] = get_water_iops(wl,temp,sal);

ap_int = interp1(apwl,ap,wl,'linear','extrap');

%Determine diffuse attenuation coefficient (K_Lu)
K_Lu=(a_sw+ap_int)./mu_u;

% --> add self shading here if desired and compute the 'true' 
% radiance using measured radiance:
%delta = shading_corr(K_Lu,r,z,lat,lon,tt2);
% for iii=1:max(size(Lu(:,1)))
%     Lu_true(iii,:)=Lu(iii,:).*delta;
% end
% Lu = Lu_true;

%extrapolate Lu to Lu(0-) to then calculate Lw 
for ilu=1:max(size(Lu_mid(:,1)))
    Lu_0(ilu,:)=Lu_mid(ilu,:).*exp(-K_Lu.*z);
end

Lw = t*Lu_0./(n_wat^2);

%Calculate: Rrs=Lw/Ed and then avg all Rrs spectra
% also calculate rrs (just below the surface)
Rrs = [];
rrs = [];
for irrs=1:max(size(Ed_mid(:,1)))
    Rrs(irrs,:)=Lw(irrs,:)./Ed_mid(irrs,:);
    rrs(irrs,:)=Lu_0(irrs,:)./Ed_mid(irrs,:);
end
Rrs = Rrs';
rrs = rrs';

for kk=1:max(size(Rrs(:,1)))
    Rrs_std(kk)=nanstd(Rrs(kk,:));
    Rrs_avg(kk)=nanmean(Rrs(kk,:));
    Rrs_med(kk)=nanmedian(Rrs(kk,:));
    
    rrs_std(kk)=nanstd(rrs(kk,:));
    rrs_avg(kk)=nanmean(rrs(kk,:));
    rrs_med(kk)=nanmedian(rrs(kk,:));    
end

% Calculate chl from band ratio of the new Rrs spectrum 
[~,~,chlOCnew] = get_ap_oc_from_Rrs(Rrs_avg,wl);

%disp(['chlOC: ',num2str(chlOC)])
disp(['Chl from OC4: ',num2str(chlOCnew)])
