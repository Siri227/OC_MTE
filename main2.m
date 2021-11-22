close all; clear; clc;
% script to show dispersion in SMF-28 w.r.t length and compensation of
% pulse signals depicting ISI with the help of FBG
%% Making structure of SMF with very low attenuatoin loss (for better plot view) 
smf = optic_fiber_parameter(18,0.092,0.0018,0,1.4682);

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
ylabel('amplitude')
title('input pulse')
legend('combined wave packet',['w = ' num2str(wavelengths(40))],...
    ['w = ' num2str(wavelengths(50))],['w = ' num2str(wavelengths(70))]);

%% optic fiber 
% length of optical fibe
Lf = [5,10,50];

%calling smf fxn for asserting dispersion and attenuation in pulse
out_pulse1 = smfloss(in_pulse,Lf(1),smf.alpha,smf.So);
pulse2 = sum(out_pulse1.pulse);
out_pulse2 = smfloss(in_pulse,Lf(2),smf.alpha,smf.So);
pulse3 = sum(out_pulse2.pulse);
out_pulse3 = smfloss(in_pulse,Lf(3),smf.alpha,smf.So);
pulse4 = sum(out_pulse3.pulse);

%ploting
figure
%ploting input vs output wave of smf

%for Lf = 0km
patch(t,pulse1,'blue','FaceAlpha',0.5)
hold on

%for Lf = 5km
patch(t,pulse2,'green','FaceAlpha',0.5)

%for Lf = 10km
patch(t,pulse3,'red','FaceAlpha',0.5)

%for Lf = 50km
patch(t,pulse4,'yellow','FaceAlpha',0.5)

xlabel('time in ps')
ylabel('power')
title('dispersoin comparision at different L_{f}')
legend([repmat('L = ',4,1) num2str([0 Lf]')]);


%% checking for ISI with bandwidth of 1GHz
% making pulse signals with each pulse of 1ns width 0010100100
s = [0;0;1;0;1;0;0;1;0;0];
y = circshift(input_pulse,-2.5e3) + circshift(input_pulse,2.5e3) + circshift(input_pulse,-0.5e3);
% making pusle matrix with same optical spectrum
pulse = zeros(length(wavelengths),length(t));
pulse = (pulse + y).*spc';

pulse1 = sum(pulse);

in_pulse = struct('t',t,'pulse',pulse,'spectrum',spectrumi,'wavelengths',wavelengths);

% max optica fiber length to support 1GHz bandwidth (km)
Lfmax = (1e12/(2*smf.Dt*4*1e9));

% output at Lfmax
out_pulse1 = smfloss(in_pulse,Lfmax,smf.alpha,smf.So);
pulse2 = sum(out_pulse1.pulse);

% outputat 3*Lfmax
out_pulse2 = smfloss(in_pulse,3*Lfmax,smf.alpha,smf.So);
pulse3 = sum(out_pulse2.pulse);

%plotting
figure
subplot(311) %input
plot(t,pulse1)
xlabel('ps'); ylabel('amp'); ylim([0 1.5]);
text(5e2:1e3:1e4,repmat(1.1,1,10),num2str(s),'color','red')
title('input signal')

subplot(312) %at Lf = Lfmax
plot(t,pulse2)
xlabel('ps'); ylabel('amp'); ylim([0 1.5]);
text(5e2:1e3:1e4,repmat(1.1,1,10),num2str(s),'color','red')
title(['at L = ' num2str(Lfmax) 'km: max length for 1GHz'])

subplot(313) %at Lf = 3*Lfmax
plot(t,pulse3)
xlabel('ps'); ylabel('amp'); ylim([0 1.5]);
text(5e2:1e3:1e4,repmat(1.1,1,10),num2str(s),'color','red')
hold on
patch(t,pulse3.*((t>3e3).*(t<4e3) + (t>5e3).*(t<7e3)),'red','FaceAlpha',0.5)
text(3e3,0.7,'indistinguishable from 1')
title(['at L = ' num2str(3*Lfmax) 'km'])

suptitle('ISI due to dispersion')

%% applying FBG to the signal
% making fbg spectrum
% parameters to choose lg = 1, chirp_var = -0.4, index change = 0.0004,
% epodize function = 2
figure
FBG_use = FBG_param_design(smf.Dt,smf.neff);
FBG_spectrum = FBG_use.spectrum;

fbg_out1 = compensate(out_pulse1,FBG_spectrum);
fbg_pulse1 = sum(fbg_out1.pulse);

fbg_out2 = compensate(out_pulse2,FBG_spectrum);
fbg_pulse2 = sum(fbg_out2.pulse);

%plotting
figure
subplot(211) %at Lf = Lfmax
plot(t,pulse2)
xlabel('ps'); ylabel('amp'); ylim([0 1.5]);
text(5e2:1e3:1e4,repmat(1.1,1,10),num2str(s),'color','red')
title(['at L = ' num2str(Lfmax) 'km: max length for 1GHz'])
hold on
patch(t,pulse2.*((t>3e3).*(t<4e3) + (t>5e3).*(t<7e3)),'red','FaceAlpha',0.5)
text(3.5e3,0.7,'ISI')

subplot(212) %at Lf = Lfmax
plot(t,fbg_pulse1)
xlabel('ps'); ylabel('amp'); ylim([0 1.5]);
text(5e2:1e3:1e4,repmat(1.1,1,10),num2str(s),'color','red')
title('compensated pulse')

figure
subplot(211) %at Lf = 3*Lfmax
plot(t,pulse3)
xlabel('ps'); ylabel('amp'); ylim([0 1.5]);
text(5e2:1e3:1e4,repmat(1.1,1,10),num2str(s),'color','red')
title(['at L = ' num2str(3*Lfmax) 'km'])
hold on
patch(t,pulse3.*((t>3e3).*(t<4e3) + (t>5e3).*(t<7e3)),'red','FaceAlpha',0.5)
text(3.5e3,0.7,'ISI')

subplot(212) %at Lf = Lfmax
plot(t,fbg_pulse2)
xlabel('ps'); ylabel('amp'); ylim([0 1.5]);
text(5e2:1e3:1e4,repmat(1.1,1,10),num2str(s),'color','red')
title('compensated pulse')


