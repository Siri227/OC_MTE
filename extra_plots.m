close all; clear; clc;
% scripts for extra plots
% 1. depicting importance of adopization function

%% apodization function
fprintf([repmat('#',1,100) '\n'])
% optical fiber params
Dt = 18;
neff = 1.4682;
% FBG params
Lg = 0.5;
chirp_var = -2.5;
del_neff = 0.0005;

uniform_fbg = FBG(Lg,neff,del_neff,chirp_var,1);
gaussian_fbg = FBG(Lg,neff,del_neff,chirp_var,2);
fprintf('calcultion done for spectral response for different adopization function\n');

w = uniform_fbg.w;

s = ['Lg = ' num2str(Lg) 'cm, n_{eff} = ' num2str(neff) ', \deltan = ' ...
    num2str(del_neff) ', \delta\lambda/\deltaz = ' num2str(chirp_var) 'nm/cm'];

subplot(121)
plot(w,uniform_fbg.r,'r','LineWidth',2)
xlabel('\lambda(nm)'); ylabel('Reflectence'); ylim([0 1]);
title('uniform grating')
text(1545,0.9,s,'color','red');

subplot(122)
plot(w,gaussian_fbg.r,'LineWidth',2)
xlabel('\lambda(nm)'); ylabel('Reflectence'); ylim([0 1]);
title('gaussian profile')
text(1545,0.9,s,'color','red');

suptitle('Comparision of uniform and gaussian profile grating reflectence spectrum')

%% variation of length
fprintf([repmat('#',1,100) '\n'])
%length in cms
Lg = [0.5 0.7 1];
chirp_var = -1; neff = 1.4682; del_neff = 0.0004; 
% calculating response for each Lg for uniform apodize grating
u1 = FBG(Lg(1),neff,del_neff,chirp_var,1);
u2 = FBG(Lg(2),neff,del_neff,chirp_var,1);
u3 = FBG(Lg(3),neff,del_neff,chirp_var,1);
fprintf('spectrum calculatoin done for uniform adopize grating for different length\n')
%calculating response for each Lg for gaussian apodize 
g1 = FBG(Lg(1),neff,del_neff,chirp_var,2);
g2 = FBG(Lg(2),neff,del_neff,chirp_var,2);
g3 = FBG(Lg(3),neff,del_neff,chirp_var,2);
fprintf('spectrum calculatoin done for gaussian adopize grating for different length\n')
%plotting
figure
subplot(121) %for uniform chirped grating
plot(w,u1.r,'b','LineWidth',2);
hold on
plot(w,u2.r,'g','LineWidth',2);
plot(w,u3.r,'r','LineWidth',2);
xlabel('\lambda(nm)'); ylabel('Reflectence'); ylim([0 1]); grid on;
legend([repmat('L_{g} = ',3,1) num2str(Lg') repmat('cm',3,1)])

subplot(122) %for gaussian grating
plot(w,g1.r,'b','LineWidth',2);
hold on
plot(w,g2.r,'g','LineWidth',2);
plot(w,g3.r,'r','LineWidth',2);
xlabel('\lambda(nm)'); ylabel('Reflectence'); ylim([0 1]); grid on;
legend([repmat('L_{g} = ',3,1) num2str(Lg') repmat('cm',3,1)])

s = ['n_{eff} = ' num2str(neff) ', \deltan = ' ...
    num2str(del_neff) ', \delta\lambda/\deltaz = ' num2str(chirp_var) 'nm/cm'];
suptitle(s)


%% variation of chirp variable
fprintf([repmat('#',1,100) '\n'])
%chirp variable in nm/cm
chirp_var = [-2.5 -1.5 -0.5];
Lg = 0.5; neff = 1.4682; del_neff = 0.0004; 
% calculating response for each Lg for uniform apodize grating
u1 = FBG(Lg,neff,del_neff,chirp_var(1),1);
u2 = FBG(Lg,neff,del_neff,chirp_var(2),1);
u3 = FBG(Lg,neff,del_neff,chirp_var(3),1);
fprintf('spectrum calculatoin done for uniform adopize grating for different chirp variable\n')
%calculating response for each Lg for gaussian apodize 
g1 = FBG(Lg,neff,del_neff,chirp_var(1),2);
g2 = FBG(Lg,neff,del_neff,chirp_var(2),2);
g3 = FBG(Lg,neff,del_neff,chirp_var(3),2);
fprintf('spectrum calculatoin done for gaussian adopize grating for different chirp vairable\n')
%plotting
figure
subplot(121) %for uniform chirped grating
plot(w,u1.r,'b','LineWidth',2);
hold on
plot(w,u2.r,'g','LineWidth',2);
plot(w,u3.r,'r','LineWidth',2);
xlabel('\lambda(nm)'); ylabel('Reflectence'); ylim([0 1]); grid on;
legend([repmat('\delta\lambda/\deltaz = ',3,1) num2str(chirp_var') repmat('nm/cm',3,1)])

subplot(122) %for gaussian grating
plot(w,g1.r,'b','LineWidth',2);
hold on
plot(w,g2.r,'g','LineWidth',2);
plot(w,g3.r,'r','LineWidth',2);
xlabel('\lambda(nm)'); ylabel('Reflectence'); ylim([0 1]); grid on;
legend([repmat('\delta\lambda/\deltaz = ',3,1) num2str(chirp_var') repmat('nm/cm',3,1)])

s = ['Lg = ' num2str(Lg) 'cm, n_{eff} = ' num2str(neff) ', \deltan = ' ...
    num2str(del_neff)];
suptitle(s)


%% variation of index change
fprintf([repmat('#',1,100) '\n'])
% index change
del_neff = [0.0002 0.0004 0.0006];
Lg = 0.5; neff = 1.4682; chirp_var = -1; 
% calculating response for each Lg for uniform apodize grating
u1 = FBG(Lg,neff,del_neff(1),chirp_var,1);
u2 = FBG(Lg,neff,del_neff(2),chirp_var,1);
u3 = FBG(Lg,neff,del_neff(3),chirp_var,1);
fprintf('spectrum calculatoin done for uniform adopize grating for different index change\n')
%calculating response for each Lg for gaussian apodize 
g1 = FBG(Lg,neff,del_neff(1),chirp_var,2);
g2 = FBG(Lg,neff,del_neff(2),chirp_var,2);
g3 = FBG(Lg,neff,del_neff(3),chirp_var,2);
fprintf('spectrum calculatoin done for gaussian adopize grating for different index change\n')
%plotting
figure
subplot(121) %for uniform chirped grating
plot(w,u1.r,'b','LineWidth',2);
hold on
plot(w,u2.r,'g','LineWidth',2);
plot(w,u3.r,'r','LineWidth',2);
xlabel('\lambda(nm)'); ylabel('Reflectence'); ylim([0 1]); grid on;
legend([repmat('\deltan = ',3,1) num2str(del_neff')])

subplot(122) %for gaussian grating
plot(w,g1.r,'b','LineWidth',2);
hold on
plot(w,g2.r,'g','LineWidth',2);
plot(w,g3.r,'r','LineWidth',2);
xlabel('\lambda(nm)'); ylabel('Reflectence'); ylim([0 1]); grid on;
legend([repmat('\deltan = ',3,1) num2str(del_neff')])

s = ['Lg = ' num2str(Lg) 'cm, n_{eff} = ' num2str(neff) ', \delta\lambda/\deltaz = ' num2str(chirp_var) 'nm/cm'];
suptitle(s)







