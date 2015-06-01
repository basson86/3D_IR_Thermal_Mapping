function [Tc]= TempCorr_dielectric_AnodizedAL(viewing_angle,Tm,Ta)
% dieletric directional emissivity simulation for Anodized AL ( using the parameter of
% "Theoretical Modeling of Skin Emissivity, 1992- by Shahram Hejazi and
%  Robert Spangler

% convert the viewing angle to radius unit
phi=viewing_angle/90*pi/2;

% optical properties of Anodized Aluminum
n1=2.84;
n2=1;
alpha=3*10^(-3)*10^5; % m^(-1)   
lamda=4*10^(-6);% m

k1= alpha*lamda/4;
E0= 1-((n1-1)/(n1+1))^2;

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

%% computing corrected temeprature 
Tc = (Tm-(rho_phi*Ta))/(1-rho_phi);
