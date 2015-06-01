function [Tc]= deltaT_by_skin_viewing_angle_D(viewing_angle,Tm,Ta)

phi=viewing_angle/90*pi/2;

% optical properties of skin & water
n1=1.33;
%n2=1;
alpha=10^5; % m^(-1)   
lamda=10*10^(-6);% m

k1= alpha*lamda/4;
E0= 1-((n1-1)/(n1+1))^2;

beta=sqrt(n1^2-sin(phi).^2);
rho_phi=1/2.*(((beta-cos(phi))./(beta+cos(phi))).^2).*(1+((beta.*cos(phi)- sin(phi).^2)./(beta.*cos(phi)+sin(phi).^2)).^2);
emis= 1-rho_phi;

%% computing corrected temeprature 
Tc = (Tm-(rho_phi*Ta))/(1-rho_phi);
