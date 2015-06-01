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

%%
Tamb=23; %(unit:degC) the measured ambient temperature for correction (based on the background of the objecti in the IR image)
         % it was set to 22.63 in the ASME paper : curvature effect quantification

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
     Del_T_Noraml((i-initial_T)/int_T+1,:)= (MeanT_Vertical((i-initial_T)/int_T+1,:)-MaxT)./(Tamb-MaxT);

     Del_T_Emis((i-initial_T)/int_T+1,:)= (MeanT_Vertical((i-initial_T)/int_T+1,:)-B_MaxT)./(Tamb-B_MaxT);

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
theta_half= asind((0:1:right_B-center)./radius);
plot(theta_half, Del_T_Noraml(:,center:1:right_B),'.-');
xlim([0 90]); 
ylim([0 1]);
grid on;
xlabel('Viewing Angle (deg)');
ylabel('Tn-T(\theta)/Tn-Ta');
legend('34 ^oC','36 ^oC','38 ^oC','40 ^oC','42 ^oC','44 ^oC');
grid on;


figure;
plot(theta_half, Del_T_Emis(:,center:1:right_B),'.-');
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
[AX,H1,H2] = plotyy(theta_half, Avg_DelT_Emis(:,center:1:right_B),phi*180/pi, rho_phi,'plot');
set(get(AX(1),'Ylabel'),'String','Ts-T(theta)/Ts-Ta') ;
set(get(AX(2),'Ylabel'),'String','1-Emissivity(theta)') ;
set(H1,'LineStyle','*');
xlabel('Viewing Angle (deg)');
legend('mean of experimental data','simulation from non-dielectric model')
grid on;

figure(6);
plot(theta_half, Avg_DelT_Emis(:,center:1:right_B),'bx');
hold on;
plot(phi*180/pi, 1.1*rho_phi,'g');
legend('mean of experimental data','scaling approximation from non-dielectric model')
grid on;

%% load the "36 degC experimental" image data to see how the formula-based correction works
pdata= 'C:\Users\Heat Transfer Lab\Desktop\ReconstructMe\0909 cyliner pot experiment\Angle correction simulation for pot image (for 3D IR)\36';
path(path, pdata);

stamp='125249';
[ImR36]= Read_IR_raw_data(stamp,0);
[ImT36]= TempConvert(ImR36);

ImT36_corr=ImT36;
%% take the central line temperature as the true temperature reading, to find the scaling factor


Normal_Temp=ImT36(upper_B:lower_B,center);
Normal_Temp_corr=CorrT_by_Viewing_Angle_Anodized_AL(0,ImT36(upper_B:lower_B,center),Tamb);

% take mean value of both "Normal_Temp" and "Normal_Temp_corr" to get a
% correction scaling factor: "Corr_Scale", which is based on their ratio

Mean_Normal_T=mean(Normal_Temp);
Mean_Normal_corr_T=mean(Normal_Temp_corr);
Corr_Scale=Mean_Normal_T/Mean_Normal_corr_T;

% apply the apprximated curvature error correction to the IR image: ImT36 
for i=left_B:right_B;
    
    % compute the x-dependent viwing angle along the cylinder surface, for
    % the approximation of error correction using non-dielectric model

    if i< center
    Viewing_angle(:,i)= -asind((i-center:i-center)./radius);
    else 
    Viewing_angle(:,i)= asind((i-center:i-center)./radius);
    end
    
    % compute the corrected temeprature based on the viewing angle and
    % non-dielectric model approximation, and scale it with "Corr_Scale" to
    % make it more truthful
    ImT36_corr(:,i)=CorrT_by_Viewing_Angle_Anodized_AL(Viewing_angle(:,i),ImT36(:,i),Tamb)*Corr_Scale;
    
end

%% compare the IR image before/after the viewing angle-based correction
figure (7) % before applying the correction
subplot(1,2,1);
imagesc(ImT36);
caxis([25 35]);
colorbar;
hold on;
%mark the range of interset in the IR image
plot(left_B+5:1:right_B-5,upper_B,'b.','MarkerSize',6);
hold on;
plot(left_B+5:1:right_B-5,lower_B,'b.','MarkerSize',6);
hold on;
% mark the central line where we sample the temperature to get the scaling
% factor
plot(center,upper_B:lower_B,'k.','MarkerSize',6);

subplot(1,2,2); % after applying the correction
imagesc(ImT36_corr);
caxis([25 35]);
colorbar;
hold on;
%mark the range of interset in the IR image
plot(left_B+5:1:right_B-5,upper_B,'b.','MarkerSize',6);
hold on;
plot(left_B+5:1:right_B-5,lower_B,'b.','MarkerSize',6);
hold on;
% mark the central line where we sample the temperature to get the scaling
% factor
plot(center,upper_B:lower_B,'k.','MarkerSize',6);

% take the mean temperature between the upper & lowe bound of ImT36
MeanT_Vertical_36=mean(ImT36(upper_B:lower_B,:));
MeanT_Vertical_36_corr=mean(ImT36_corr(upper_B:lower_B,:));

% examine the 1D line along the curved surface befre/after correciton
figure (8);
subplot(1,2,1);
plot(theta(6:length(theta)-5),MeanT_Vertical_36(:,left_B+5:1:right_B-5),'b');
xlim([-70,70]);
ylim([25 35]);
grid on;
ylabel('Temperature (^oC)');
subplot(1,2,2);
plot(theta(6:length(theta)-5),MeanT_Vertical_36_corr(:,left_B+5:1:right_B-5),'b');
xlim([-70,70]);
ylim([25 35]);
grid on;
ylabel('Temperature (^oC)');
