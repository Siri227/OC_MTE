%test smfloss 
close all
clear

t = linspace(0,1e3,1e3);
pulsei = @(t)gaussmf(t,[20 500]);
wavelengths = linspace(1550-50,1550+50,1e3);
spectrumi = gaussmf(wavelengths,[0.05 1550]);

%input pulse structure
in_pulse = struct('pulse',pulsei,'spectrum',spectrumi,'wavelengths',wavelengths);

%calling smf loss function for diff lengths
out_pulse1 = smfloss(in_pulse,10);
out_pulse2 = smfloss(in_pulse,50);
out_pulse3 = smfloss(in_pulse,100);

%ploting
plot(t,pulsei(t));
hold on
plot(t,out_pulse1.pulse(t));
plot(t,out_pulse2.pulse(t));
plot(t,out_pulse3.pulse(t));
legend('input wave','L = 10','L = 50','L = 100');
xlabel('time in ps')
ylabel('normalized intensity')
