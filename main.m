close all; clear; clc;
% OC MTE Project
% Main script depticting fiber dispersion and attenuation loss and dispersoin
% conepnsation with the use of FBG grating fiber. SMF-28 fiber is used for default parameters
% and gaussian profile is used for optical laser pulse
% the following opitcal system is maintaining at least 1GHz bandwidth
% wavelength span for whole simulaton is 1545nm to 1555nm

%% Making structure for SMF-28 fiber
smf28 = optic_fiber_parameter();


%% Generating optical laser pulse
% time axis in ps
t = linspace(0,1e4,1e4);
%pulse of 1ns width
input_pulse = (t > 4.5e3).*(t < 5.5e3);
%generating optical specturm for pulse
wavelengths = linspace(1550-5,1550+5,100); %wavelength array in nm
%FWHM 4nm, normalized intensity in dB
sigma = 4/(sqrt(8*log(2)));
spectrumi = 10*log(gaussmf(wavelengths,[sigma 1550]));
%pulse matrix
spc = (10.^(spectrumi/10))/(sigma*sqrt(2*pi))*(wavelengths(2)-wavelengths(1));
pulse = zeros(length(wavelengths),length(t));
pulse = (pulse + input_pulse).*spc';


%input pulse structure
in_pulse = struct('t',t,'pulse',pulse,'spectrum',spectrumi,'wavelengths',wavelengths);

pulse1 = sum(pulse);
plot(t,pulse1);
hold on
plot(t,pulse(40,:))
plot(t,pulse(50,:))
plot(t,pulse(70,:))
xlabel('time in ps')
ylabel('normalized power')
title('input pulse')
legend('combined wave packet',['w = ' num2str(wavelengths(40))],...
    ['w = ' num2str(wavelengths(50))],['w = ' num2str(wavelengths(70))]);

%% optic fiber 
%bandwidth-length product with no spectral modification from FBG
FWHM = 4;
disp_per_km = smf28.Dt*FWHM; %(ps/km)
BWL = 1e3/(2*disp_per_km); %(GHz.km)
Lf = floor(BWL); %in km
% Lf = 20;

%calling smf fxn for asserting dispersion and attenuation in pulse
out_pulse1 = smfloss(in_pulse,Lf);

%ploting
figure
subplot(211)
%ploting input vs output wave of smf
res_out_pulse = sum(out_pulse1.pulse);
plot(t,res_out_pulse)
hold on
plot(t,pulse1);
legend(['L =' num2str(Lf)],'input wave');
xlabel('time in ps')
ylabel('normalized power')
title('input pulse vs output pulse of smf')

subplot(212)
%ploting optical spectrum of input and output wave
plot(wavelengths,in_pulse.spectrum);
hold on
plot(wavelengths,out_pulse1.spectrum);
xlabel('wavelength in nm')
ylabel('normalized power in dB')
legend('in wave','out wave');
title('spectrum of input and output wave')


figure
plot(t,out_pulse1.pulse(10:10:end-10,:));
legend(num2str(wavelengths(10:10:end-10)'))
title('pulse components')

%% FBG design
figure
FBG_use = FBG_param_design(smf28.Dt,smf28.neff);
FBG_spectrum = FBG_use.spectrum;
% FBG_spectrum = FBG(1,smf28.neff,0.0004,-0.4,2);
% figure
% subplot(211)
% plot(FBG_spectrum.w,FBG_spectrum.r)
% xlabel('wavelength nm')
% ylabel('reflectivity')
% title('reflectivity')
% subplot(212)
% plot(FBG_spectrum.w,FBG_spectrum.tau)
% xlabel('wavelength nm')
% ylabel('time delay')
% title('time delay')

out_pulse2 = compensate(out_pulse1,FBG_spectrum);

res_fbg_out = sum(out_pulse2.pulse);

figure
subplot(211)
plot(t,res_out_pulse);
hold on
plot(t,res_fbg_out);
xlabel('time in ps')
ylabel('normalized power')
legend('FBG in','FBG out')
title('FBG in vs FBG out pulse')

subplot(212)
plot(wavelengths,out_pulse1.spectrum);
hold on
plot(wavelengths,out_pulse2.spectrum);
xlabel('wavelength in nm')
ylabel('normalized power in dB')
legend('FBG in','FBG out')
title('optical spectrum') 

figure
subplot(131)
patch(t,input_pulse,'blue','FaceAlpha',0.4)
xlabel('time in ps'); ylabel('amp'); ylim([0 1.2]); 
title('input signal')

subplot(132)
patch(t,res_out_pulse/max(res_out_pulse),'red','FaceAlpha',0.4)
xlabel('time in ps'); ylabel('amp'); ylim([0 1.2]); 
title(['after smf of L_{f}: ' num2str(Lf) 'km'])

subplot(133)
patch(t,res_fbg_out/max(res_fbg_out),'green','FaceAlpha',0.4)
xlabel('time in ps'); ylabel('amp'); ylim([0 1.2]); 
title('after FBG')

suptitle('FBG compensation result')






