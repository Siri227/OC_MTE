% script for pulse dispersion and attenuatoin 
% input parameters: optical laser pulse struct
%                   length of fiber in km
%                   attenuation parameter (alpha) db/km
%                   second order disperson parameter (So) in ps/nm2.km
% output = modified optic pusle struct
% optic struct field:   pulse in time domain (ps): function handle 
%                       optical spectrum: array
%                       wavelengths axis: array
function [out_pulse] = smfloss(in_pulse,L,varargin)
nopin = length(varargin);
if nopin > 2
    error('Too many input arguements in smfloss')
end

defaults = {0.14 + 0.04*rand(1),0.092};
defaults(1:nopin) = varargin;

[alpha,So] = defaults{:};

%% dispersion calculation
%rms spectrum width of input pulse
sigmaw = pulsewidth(in_pulse.spectrum,in_pulse.wavelengths);
%total dispersoin in ps/nm.km
Dt = (So*1550/4)*(1 - (1310/1550)^4);
%pulse broadening calculation in ps
pulse_inc = sigmaw*L*Dt;

%% signal attenuatoin in db
dbloss = alpha*L;
af = 10^(-dbloss/10);
%% generating output pulse
t = linspace(0,1e4,1e4);
pulsei = in_pulse.pulse(t);
%max value index
[~,t0] = max(pulsei);
%scaling factor
sf = pulsewidth(pulsei,t)/(pulsewidth(pulsei,t) + pulse_inc);
%shifting to center 
fc = @(t)in_pulse.pulse(t+t0);
%scaling
fs = @(t)fc(t*sf);
%shifting back
pulsef = @(t) af*fs(t-t0);

spectrumf = af*in_pulse.spectrum;
out_pulse = struct('pulse',pulsef,'spectrum',spectrumf,'wavelengths',in_pulse.wavelengths);

fprintf('rms spectrum width:'); disp(sigmaw);
fprintf('Dispersion parameter:'); disp(Dt);
fprintf('dispersion inc/km:'); disp(pulse_inc/L);
end


