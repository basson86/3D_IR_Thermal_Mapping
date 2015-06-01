function [Tc]= deltaT_by_skin_viewing_angle_Cosine(viewing_angle,Tm,Ta)

phi=viewing_angle/90*pi/2;

emis= cos(phi);

rho_phi=1-emis;
%% computing corrected temeprature 
Tc = (Tm-(rho_phi*Ta))/(1-rho_phi);
