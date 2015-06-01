
clear all; close all;

% optical properties of skin & water
n1=1.33;
n2=1;
alpha=10^5; % m^(-1)   
lamda=10*10^(-6);% m

% directional emissivity simulation for human skin ( using the parameter of
% "Theoretical Modeling of Skin Emissivity, 1992- by Shahram Hejazi and
%  Robert Spangler
k1= alpha*lamda/4;
E0= 1-((n1-1)/(n1+1))^2;
[phi,rho_phi,emis_ND]=Direct_Emiss_non_dielectric(n1,n2,alpha,lamda);
[phi,rho_phi,emis_D]=Direct_Emiss_dielectric(n1);

%% 
figure(1);
polar(phi,emis_ND, 'r');
hold on;
polar(phi,emis_D, 'b');
hold on;
polar(phi,E0*cos(phi), 'g');
hold on;
polar(phi,E0*ones(1,158),'m');
view([90 -90]);


xl = get(gca,'XLim'); yl = get(gca,'YLim');
set(gca,'XLim', [0 xl(2)], 'YLim', [0 yl(2)]);
legend('Non-dielectric model','Dielectric model','Cosine function');


figure(2);
plot(phi*180/pi, emis_ND,'r');
hold on;
plot(phi*180/pi, emis_D,'b');
hold on;
plot(phi*180/pi, E0*cos(phi),'g');
hold on;
plot(phi*180/pi, E0*ones(1,158),'m');

xlabel('Viewing Angle (deg)');
ylabel('(Reflectivity): 1-Emissivity(theta)') ;
legend('Non-dielectric model','Dielectric model','Cosine function');

grid on;

