% directional emissivity simulation for human skin ( using the parameter of
% "Theoretical Modeling of Skin Emissivity, 1992- by Shahram Hejazi and
%  Robert Spangler
clear all; close all;

n1=1.33;
n2=1;
alpha=10^5; % m^(-1)   
lamda=10*10^(-6);% m

k1= alpha*lamda/4;
E0= 1-((n1-1)/(n1+1))^2;
[phi,rho_phi,emis]=Direct_Emiss(n1,n2,alpha,lamda);

% 
figure;
plot(phi*180/pi, emis);
xlabel('Viewing Angle (deg)');
ylabel('1-Emissivity(theta)') ;
grid on;
