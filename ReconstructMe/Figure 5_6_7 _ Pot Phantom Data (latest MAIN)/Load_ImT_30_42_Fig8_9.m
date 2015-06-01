clear all; close all; clc;
% load the averaged IR image over time at different temperature
initial_T=34;
final_T=44;
int_T=2;

% the vertical range for averaging
upper_B= 75;
lower_B= 175;
sample_row=(upper_B+lower_B)/2;

% the herizontal range of the curve surface
left_B=8;
right_B=312;
center=round((left_B+right_B)/2);
radius=(right_B-left_B)/2;

marked_deg_ind=[33,40,59,84,113,142,169,199,226,253,274,290,300];

%figure;
%subplot(1,2,2);
% for k=1:length(marked_deg_ind);
% plot(marked_deg_ind(1,k),22:40,'r.');
% hold on;
% end
% hold on;

for i= initial_T:int_T:final_T
%%
    % load Averaged Image at temperature i
    load(strcat('Aver_ImT_',num2str(i)));
    %[Aver_ImT,error,f]=TempConvertCorret_0306(Aver_ImT);
    
     % take the mean value of temeprature over the vertical direction
     MeanT_Vertical((i-initial_T)/int_T+1,:)=mean(Aver_ImT(upper_B:lower_B,:));
     Norm_MeanT_Vertical((i-initial_T)/int_T+1,:)= MeanT_Vertical((i-initial_T)/int_T+1,:)./max(MeanT_Vertical((i-initial_T)/int_T+1,:));
     
     % assuming the maximum temperature is at the center of the pot
     MaxT= MeanT_Vertical((i-initial_T)/int_T+1,center);
     
     B_MaxT= i; % assuming the water temperature is the temperature from blackbody
     
     % approximate the emissivity at each temperature level
     Emi_appr((i-initial_T)/int_T+1,:)= MeanT_Vertical((i-initial_T)/int_T+1,center)/B_MaxT;
     
     Del_T((i-initial_T)/int_T+1,:)= -(MeanT_Vertical((i-initial_T)/int_T+1,:)-MaxT);
     Del_T_Noraml((i-initial_T)/int_T+1,:)= (MeanT_Vertical((i-initial_T)/int_T+1,:)-MaxT)./(22.63-MaxT);

     Del_T_Emis((i-initial_T)/int_T+1,:)= (MeanT_Vertical((i-initial_T)/int_T+1,:)-B_MaxT)./(22.63-B_MaxT);

     % compute theta according to projected position
     theta= asind((left_B-center:right_B-center)./radius);
     
% Show the movie frame
%     subplot(1,2,1);
%     imagesc(Aver_ImT,[21 40]); 
%     hold on;
%     plot(1:320,upper_B,'-');
%     hold on;
%     plot(1:320,lower_B,'-');
%     hold on;
%     plot(left_B,1:256,'-');
%     hold on;
%     plot(right_B,1:256,'-');


%     subplot(1,2,2);
%     plot(MeanT_Vertical((i-initial_T)/int_T+1,:));
%     hold on;
%     xlim([1 320]);
%     ylim([22 40]);
%     grid on;
%     ylabel('degC');

% temperature distribution along the curved surface
 %    subplot(1,2,2);    
%present in positive angle only (folding together)
%     theta= asind((0:right_B-center)./radius);
%     plot(theta,  MeanT_Vertical((i-initial_T)/int_T+1,center:1:right_B),'b');
%     hold on;
%     plot(theta,  MeanT_Vertical((i-initial_T)/int_T+1,center:-1:left_B),'b--');
%     grid on;
% symetrical view

    plot(theta,  MeanT_Vertical((i-initial_T)/int_T+1,left_B:1:right_B),'b');
    hold on;
    xlim([-90 90]); 
    ylim([20 40]);
    grid on;
    xlabel('Viewing Angle (deg)');
    ylabel('Temperature (^oC)');
    hold on;
    
% get frame for replay the movie    
%F((i-initial_T)/int_T+1)=getframe(gcf);   
end
%%
% figure;
% theta= asind((left_B-center:1:right_B-center)./radius);
% plot(theta, MeanT_Vertical(1,left_B:1:right_B),'b');
% ylim([27 40]);

%figure;
% plot(theta, Norm_MeanT_Vertical(1:7,left_B:1:right_B));
% xlim([-90 90]); 
% ylim([0 1.5]);
% grid on;
% xlabel('incident angle (deg)');
% ylabel('T/Tmax (degC)');
% legend('30','32','34','36','38','40','42','44');
%%
figure;
%theta= asind((left_B-center:1:right_B-center)./radius);
%plot(theta, Del_T_Noraml(:,left_B:1:right_B));
theta= asind((0:1:right_B-center)./radius);
plot(theta, Del_T_Noraml(:,center:1:right_B),'.-');
xlim([0 90]); 
ylim([0 1]);
grid on;
xlabel('Viewing Angle (deg)');
ylabel('Tn-T(\theta)/Tn-Ta');
legend('34 ^oC','36 ^oC','38 ^oC','40 ^oC','42 ^oC','44 ^oC');
grid on;


figure;
theta= asind((0:1:right_B-center)./radius);
plot(theta, Del_T_Emis(:,center:1:right_B),'.-');
xlim([0 90]); 
ylim([0 1]);
grid on;
xlabel('Viewing Angle (deg)');
ylabel('Tw-T(\theta)/Tw-Ta');
legend('34 ^oC','36 ^oC','38 ^oC','40 ^oC','42 ^oC','44 ^oC');
grid on;



%% Theoritical Modelling
n1=2.84;
n2=1;
alpha=3*10^(-3)*10^5; % m^(-1)   
lamda=4*10^(-6);% m

k1= alpha*lamda/4;
E0= 1-((n1-1)/(n1+1))^2;
[phi,rho_phi,emis]=Direct_Emiss(n1,n2,alpha,lamda);


% comparison with experimental data
Avg_DelT_Emis= mean(Del_T_Emis,1);
% 
figure(5);
[AX,H1,H2] = plotyy(theta, Avg_DelT_Emis(:,center:1:right_B),phi*180/pi, rho_phi,'plot');
set(get(AX(1),'Ylabel'),'String','Ts-T(theta)/Ts-Ta') ;
set(get(AX(2),'Ylabel'),'String','1-Emissivity(theta)') ;
set(H1,'LineStyle','*');
xlabel('Viewing Angle (deg)');
legend('mean of experimental data','simulation from non-dielectric model')
grid on;

figure(6);
plot(theta, Avg_DelT_Emis(:,center:1:right_B),'bx');
hold on;
plot(phi*180/pi, 1.1*rho_phi,'g');
legend('mean of experimental data','scaling exproximation from non-dielectric model')
grid on;

%%

