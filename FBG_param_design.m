% function for designing the FBG
% input parameters:   dispersion parameter (ps/nm.km)
%                     neff for 1550nm
% output structure:   FBG parameter struct: length of FBG fiber (in cm)
%                                           chirp variable (in nm/cm)
%                                           del_neff
%                                           epodize function choice
%                                           spectrum
%                                           max length of fiber FBG can
%                                           compensate
function [FBG_param] = FBG_param_design(Dt,neff)
if nargin>2
    error('FBG_desing: too many input arguements')
elseif nargin<2
    error('FBG_design: not enough input arguements')
end

fprintf([repmat('#',1,20) 'FBG Design' repmat('#',1,50) '\n']);
% if input recieve empty value, default value get selected

fprintf('insert design parameters for FBG\n')
Lg = input('length of FBG grating (in cm): '); 
chirp_var = input('chirp variable (nm/cm): ');
del_neff = input('index change: ');
ep_fxn_choice = input('epodize function choice: ');


% setting up defalult values
if isempty(Lg)
    Lg = 1;
end
if isempty(chirp_var)
    chirp_var = -0.4;
end
if isempty(del_neff)
    del_neff = 0.0004;
end
if isempty(ep_fxn_choice)
    ep_fxn_choice = 2;
end

% generating optical spectrum for FBG
fbg_spectrum = FBG(Lg,neff,del_neff,chirp_var,ep_fxn_choice);
plot(fbg_spectrum.w,fbg_spectrum.r);
% FWHM in nm
FWHM = pulsewidth(fbg_spectrum.r,fbg_spectrum.w);
if isempty(FWHM)
    error('FBG_design: FWHM cannot be calculated using pulsewidth')
end

%calculating Q factor of reflectance spectrum
mx = max(fbg_spectrum.r);
index_hm = find(abs(fbg_spectrum.r - mx/2) < 5e-2);
if isempty(index_hm)
    error('fbg_param_design: cannot find index of half maximum -> decrease precision')
end

if length(index_hm) == 1
    error('fbg_param_design: decrease precision')
end

index_l = index_hm(1);
index_r = index_hm(end);
index_c = round((index_l+ index_r)/2);
wavelength_c = fbg_spectrum.w(index_c);
Q = wavelength_c/FWHM;

s = ['Lg = ' num2str(Lg) 'cm, n_{eff} = ' num2str(neff) ', \deltan = ' ...
    num2str(del_neff) ', \delta\lambda/\deltaz = ' num2str(chirp_var) 'nm/cm'];

%plotting 
subplot(211)
plot(fbg_spectrum.w,fbg_spectrum.r)
xlabel('wavelength(nm)')
ylabel('reflectence coeff')
title(['Q factor: ' num2str(Q) ', \lambda_{c}: ' num2str(wavelength_c) 'nm'])
hold on
patch(fbg_spectrum.w,(fbg_spectrum.w > fbg_spectrum.w(index_l)).*(fbg_spectrum.w < fbg_spectrum.w(index_r)),...
    'blue','FaceAlpha',0.2,'EdgeColor','none') 
text(1545,0.9,s,'color','red');


subplot(212)
plot(fbg_spectrum.w,fbg_spectrum.tau)
xlabel('wavelength(nm)')
ylabel('time delay (ps)')
title(['\Delta\tau_{max}: ' num2str((2*neff*Lg)/3e-2) 'ps'])
hold on
len = max(fbg_spectrum.tau) - min(fbg_spectrum.tau);
upb = max(fbg_spectrum.tau);
lowb = min(fbg_spectrum.tau);
patch(fbg_spectrum.w,len*(fbg_spectrum.w > fbg_spectrum.w(index_l)).*(fbg_spectrum.w < fbg_spectrum.w(index_r)) + lowb,...
    'blue','FaceAlpha',0.2,'EdgeColor','none')
ylim([lowb upb])

%approx max length of fiber fbg can compensate (km)
Lf_max = (2*neff*Lg)/(3e-2*FWHM*Dt);
suptitle(['Lf_{max} : ' num2str(Lf_max) 'km'])

disp([index_l index_c index_r])

FBG_param = struct('Lg',Lg,...
                        'chirp_var',chirp_var,...
                        'del_neff',del_neff,...
                        'fxn_choice',ep_fxn_choice,...
                        'spectrum',fbg_spectrum,...
                        'Lf_max',Lf_max);


end








                    
