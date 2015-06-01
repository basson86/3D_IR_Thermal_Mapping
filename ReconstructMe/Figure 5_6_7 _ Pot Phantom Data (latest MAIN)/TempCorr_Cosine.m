function [Tc]= TempCorr_Cosine(viewing_angle,Tm,Ta)

% convert the viewing angle to radius unit


radius=viewing_angle/90*pi/2;

% Emissivity for Anodized Aluminum

E0=0.7704;

Emissivity= 0.7704 *cos(radius);
rho_phi= 1- Emissivity;

%% computing corrected temeprature 
Tc = (Tm-(rho_phi*Ta))/(1-rho_phi);
