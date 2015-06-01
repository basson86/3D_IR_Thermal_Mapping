%clear all; %close all;
function [phi,rho_phi,emis]=Direct_Emiss(n1,n2,alpha,lamda);


% optical properties of skin & water
% n1=3.1;
% n2=1;
% alpha= 10*10^5; % m^(-1)   
% lamda=4*10^(-6) ;% m

%% incident angle:
phi=0:0.01:pi/2;

%% complex refrative index : n
%k1= alpha.*lamda./(4.*phi);
k1= alpha*lamda/4;
n= n1+ i.*k1;
Sin_X=n2./(n1.^2+ k1.^2).*(n1-i*k1).*sin(phi);

%%
p_sq= 1/2*(-n1.^2+ k1.^2+n2.^2*sin(phi).^2)+ 1/2*sqrt(4*n1.^2*k1.^2+ (n1.^2-k1.^2-n2.^2*sin(phi).^2).^2);
q_sq= 1/2*(n1.^2- k1.^2-n2.^2*sin(phi).^2)+ 1/2*sqrt(4*n1.^2*k1.^2+ (n1.^2-k1.^2-n2.^2*sin(phi).^2).^2);

p=sqrt(p_sq);
q=sqrt(q_sq);

rho_n=((q-n2.*cos(phi)).^2 + p.^2)./((q+n2.*cos(phi)).^2 + p.^2);
rho_p= (((n1.^2 - k1.^2).*cos(phi)-n2.*q).^2 + (2*n1.*k1.*cos(phi)-n2.*p).^2)./...
       (((n1.^2 - k1.^2).*cos(phi)+n2.*q).^2 + (2*n1.*k1.*cos(phi)+n2.*p).^2);
   
rho_phi= (rho_n+rho_p)./2;
emis= 1-rho_phi;



