%clear;
% r -> reflection spectrum of FBG
% tau -> time delay of FBG
% lamda -> wavelength on which calculation has been made
function [r,tau,lamda] = FBG(varargin)
%% input grating parameter
% setting up default values
if length(varargin)>5
    error('Too Many Arguements in FBG function')
end

defaults = {2,0.04,1.447,10*1e-9,100};
defaults(1:nargin) = varargin;

[choice,L,n_eff,C,N] = defaults{:};


delta_neff = 0.00005;    % index change
bragg_period = 1e-4;    % bragg period
del_n = 0.0004;         % depth of modulation
lambdad = 1.55*1e-6;     % design grating wavelength

bragg_p = lambdad/(2*n_eff);   % bragg period 

%% refractive index variation along z
syms x y; % x as z bcs of fcontour
% fprintf('choose apodizatoin function:\n')
% fprintf('1: uniform grating \n2: gaussian profile \n3: raused cosine \n');
% choice = input('choose from above(frm 1 to 3): ');
% while ((choice > 3) || (choice < 1))
%     choice = input('invalid choice\nenter again: ');
% end
%apodization function
gz1 = 1;                     % uniform grating

a = 64;                      % gaussian profile
gz2 = exp(-a*((x - L/2)/L).^4);
 
gz3 = 0.5*(1 + cos(pi*(x - L/2)/L)); % raised cosine
gzs = [gz1 gz2 gz3];

gz = gzs(choice);

%% steps for simulatoin
%N;              %Segmentation of FBG
M=1501;            %Total points of computation

%% adopization function
ap_fn = double(subs(gz,linspace(0,L,N)));
delta_n_eff = delta_neff*ap_fn;

%% input wavelength
lamda1=lambdad*1e9 - 1;
lamda2=lambdad*1e9 + 1;
lamda=linspace(lamda1,lamda2,M)*1e-9; %Wavelength span
delta_lamda=(lamda2-lamda1)/M*1e-9;

%% initializing reflectence, phase, delay vector
r=zeros(1,length(lamda));
phi=zeros(1,length(lamda));
tau=zeros(1,length(lamda));

%% transmisson matrix method
for k=1:M
    F=[1,0;0,1];  %Initialize transmission matrix
    for i=1:N
        %£¨1£©Uniform
        %delta_n_eff=0.00005;
        %delta_n_eff=0.00005*exp((-64*(-L/2+i*L/N)^4)/L^4);
        lamda_D=(1550 + C*L*(i/N - 0.5))*1e-9;  %The center wavelength of each fiber grating segment
        sigma=2*pi*n_eff*(1/lamda(k)-1/lamda_D)+2*pi*delta_n_eff(i)/lamda(k)+(4*pi*n_eff)*C*(-L/2+i*L/N)/lamda_D^2;  %Self coupling coefficient
        kappa=pi*delta_n_eff(i)/lamda(k);  %AC coupling coefficient
        omega=sqrt(kappa^2-sigma^2);
        
        f11=cosh(omega*L/N)-1i*(sigma/omega)*sinh(omega*L/N);
        f12=1i*(kappa/omega)*sinh(omega*L/N);
        f21=-1i*(kappa/omega)*sinh(omega*L/N);
        f22=cosh(omega*L/N)+1i*(sigma/omega)*sinh(omega*L/N);
        F=F*[f11,f12;f21,f22];  %Transmission matrix
    end
    r(k)=(abs(F(3)/F(1)))^2;  %Magnitude of refractive rate
    phi(k)=phase(F(3)/F(1));  %Phase of refractive rate
end

dtg = [phi(1) diff(phi)];
s = std(dtg);
plot(diff(dtg))
dtg(abs([0 diff(dtg)])>s) = nan;

tau = -((lamda.^2)/(2*pi*3e8)).*(dtg/delta_lamda)*1e12;


end %function end


%% display graph
% subplot(2,1,1)
% plot(lamda*1e9,r,'b');
% grid on
% xlabel('Wavelength/nm');
% ylabel('Reflection');
% title('Reflection Spectrum of Chirped Grating');
% 
% subplot(2,1,2)
% plot(lamda*1e9,tau,'r');
% grid on
% xlabel('Wavelength/nm');
% ylabel('Time Delay');
% title('Time Delay of Chirped Grating');