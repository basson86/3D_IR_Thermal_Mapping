function [phi,rho_phi,emis]=Direct_Emiss_dielectric(n1)
% dieletric model for directional emissivity simulation( using the parameter of
% "Theoretical Modeling of Skin Emissivity, 1992- by Shahram Hejazi and
%  Robert Spangler
%% incident angle:
phi=0:0.01:pi/2;

%%
beta=sqrt(n1^2-sin(phi).^2);
rho_phi=1/2.*(((beta-cos(phi))./(beta+cos(phi))).^2).*(1+((beta.*cos(phi)- sin(phi).^2)./(beta.*cos(phi)+sin(phi).^2)).^2);
emis= 1-rho_phi;

