close all; clear; clc;
% OC MTE Project
% Main script depticting fiber dispersion and attenuation loss and dispersoin
% conepnsation with the use of FBG grating fiber. SMF-28 fiber is used for default parameters
% and gaussian profile is used for optical laser pulse
% the following opitcal system is maintaining at least 1GHz bandwidth

%% Making structure for SMF-28 fiber
smf28 = optic_fiber_parameter();


%% Generating optical laser pulse
% time axis in ps
t = linspace(0,1e4,1e4);
%pulse of 1ns width
pulsei = @(t)gaussmf(t,[500 5e3]);
%generating optical specturm for pulse
wavelengths = linspace(1550-20,1550+20,1e3); %wavelength array in nm
%FWHM 10nm, normalized intensity in dB
spectrumi = 10*log(gaussmf(wavelengths,[10/(sqrt(8*log(2))) 1550]));

%input pulse structure
in_pulse = struct('pulse',pulsei,'spectrum',spectrumi,'wavelengths',wavelengths);

%% optic fiber 
%bandwidth-length product with no spectral modification from FBG
FWHM = 10;
disp_per_km = smf28.Dt*FWHM; %(ps/km)
BWL = 1e3/(2*disp_per_km); %(GHz.km)
Lf = floor(BWL); %in km

%calling smf fxn for asserting dispersion and attenuation in pulse
out_pulse1 = smfloss(in_pulse,Lf);

%ploting
plot(t,pulsei(t));
hold on
plot(t,out_pulse1.pulse(t));
legend('input wave','L = 10km');
xlabel('time in ps')
ylabel('normalized intensity')


%% FBG design
[r,tau,w] = FBG();
out_pulse2 = compensate(out_pulse1,50,r,w);

figure
subplot(211)
plot(t,out_pulse1.pulse(t));
hold on
plot(t,out_pulse2.pulse(t));
xlabel('time in ps')
ylabel('normalized intensity')
legend('FBG in','FBG out')

subplot(212)
plot(wavelengths,out_pulse1.spectrum);
hold on
plot(wavelengths,out_pulse2.spectrum);
xlabel('wavelength in nm')
ylabel('Intensity')
legend('FBG in','FBG out')







