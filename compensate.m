%pulse compensation function
% input parameters:   optical laser pulse structure
%                     pulse compensation 
%                     reflectance coeff array 
%                     wavelength for reflectance coeff array
% ouput: modified optical pulse structure
function [out_pulse] = compensate(in_pulse,pulse_comp,r,wavelengths)
if nargin>4
    error('Too many arguement in compesate fxn')
elseif nargin<4
    error('Not enough arguement in compensate fxn')
end

%% generating output pulse
t = linspace(0,1e4,1e4);
pi = in_pulse.pulse(t);
% plot(t,pi);
% hold on
if isempty(pulsewidth(pi,t))
    error('compensate fxn: unable to find pulse width')
end
%max value index
[~,t0] = max(pi);
%scaling factor
sf = (pulsewidth(pi,t) + pulse_comp)/pulsewidth(pi,t);
%shifting to center 
fc = @(t)in_pulse.pulse(t+t0);
%scaling
fs = @(t)fc(t*sf);
%amp scaling
af = max(r);
%shifting back
pf = @(t) af*fs(t-t0);
% plot(t,pf(t))

%% generating output spectrum
w = in_pulse.wavelengths;
edge_db_left = 10*log10(r(1));
edge_db_right = 10*log10(r(end));
spectrumf = zeros(1,length(w));
for i = 1:length(w)
    index = find(abs(wavelengths - w(i))<=1e-2);
    if ~isempty(index)
        if(length(index)>1)
            spectrumf(i) = in_pulse.spectrum(i) + 10*log10(r(index(1)));
        else
            spectrumf(i) = in_pulse.spectrum(i) + 10*log10(r(index));
        end
    elseif i < length(w)/2
        spectrumf(i) = in_pulse.spectrum(i) + edge_db_left;
    else
        spectrumf(i) = in_pulse.spectrum(i) + edge_db_right;
    end
end   


%% generating output structure
out_pulse = struct('pulse',pf,'spectrum',spectrumf,'wavelengths',w);

end


