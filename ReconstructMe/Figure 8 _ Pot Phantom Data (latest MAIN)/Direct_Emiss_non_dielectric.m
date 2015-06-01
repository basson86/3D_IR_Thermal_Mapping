function [phi,rho_phi,emis]=Direct_Emiss_non_dielectric(n1,n2,alpha,lamda)
% non-dieletric model for directional emissivity simulation( using the parameter of
% "Theoretical Modeling of Skin Emissivity, 1992- by Shahram Hejazi and
%  Robert Spangler

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


