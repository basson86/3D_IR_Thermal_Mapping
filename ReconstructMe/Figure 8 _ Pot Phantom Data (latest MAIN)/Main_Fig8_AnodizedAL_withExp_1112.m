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
figure(1);
theta_half= asind((0:1:right_B-center)./radius);
plot(theta_half, Del_T_Noraml(:,center:1:right_B),'.-');
xlim([0 90]); 
ylim([0 1]);
grid on;
xlabel('Viewing Angle (deg)');
ylabel('Tn-T(\theta)/Tn-Ta');
legend('34 ^oC','36 ^oC','38 ^oC','40 ^oC','42 ^oC','44 ^oC');
grid on;


figure(2);
plot(theta_half, Del_T_Emis(:,center:1:right_B),'.-');
xlim([0 90]); 
ylim([0 1]);
grid on;
xlabel('Viewing Angle (deg)');
ylabel('Tw-T(\theta)/Tw-Ta');
legend('34 ^oC','36 ^oC','38 ^oC','40 ^oC','42 ^oC','44 ^oC');
grid on;



%% Theoritical Modelling
% the optical properties of anodized pot
n1=2.84;
n2=1;
alpha=3*10^(-3)*10^5; % m^(-1)   
lamda=4*10^(-6);% m

k1= alpha*lamda/4;
E0= 1-((n1-1)/(n1+1))^2;

[phi,rho_phi_ND,emis_ND]=Direct_Emiss_non_dielectric(n1,n2,alpha,lamda);
[phi,rho_phi_D,emis_D]=Direct_Emiss_dielectric(n1);

% comparison with experimental data
Avg_DelT_Emis= mean(Del_T_Emis,1);
% 
figure(3);
[AX,H1,H2] = plotyy(theta_half, Avg_DelT_Emis(:,center:1:right_B),phi*180/pi, rho_phi_ND,'plot');
set(get(AX(1),'Ylabel'),'String','Ts-T(theta)/Ts-Ta') ;
set(get(AX(2),'Ylabel'),'String','1-Emissivity(theta)') ;
set(H1,'LineStyle','*');
xlabel('Viewing Angle (deg)');
legend('mean of experimental data','simulation from non-dielectric model')
grid on;


figure(4);
plot(theta_half, Avg_DelT_Emis(:,center:1:right_B),'bx');
hold on;
plot(phi*180/pi, 1.1*rho_phi_ND,'r');
hold on;
plot(phi*180/pi, 1.1*rho_phi_D,'b');
legend('mean of experimental data','scaling exproximation from non-dielectric model')
grid on;
%% directional emissivity polar plot comparison (non-dielectric model, cosine model, and circular model)

figure(5);
polar(phi,emis_ND, 'r');
hold on;
polar(phi,emis_D, 'b');
hold on;
polar(phi,E0*cos(phi), 'g');
view([90 -90])

legend('Non-dielectric model','Dielectric model','cosine function')
xl = get(gca,'XLim'); yl = get(gca,'YLim');
set(gca,'XLim', [0 xl(2)], 'YLim', [0 yl(2)]);


%% load the "36 degC experimental" image data to see how the formula-based correction works
pdata= 'C:\Users\Heat Transfer Lab\Desktop\ReconstructMe\0909 cyliner pot experiment\Angle correction simulation for pot image (for 3D IR)\36';
path(path, pdata);

stamp='125249';
[ImR36]= Read_IR_raw_data(stamp,0);
[ImT36]= TempConvert(ImR36);

ImT36_corr_ND=ImT36;
ImT36_corr_CS=ImT36;
%% take the central line temperature as the true temperature reading, to find the scaling factor
n1=2.84;
n2=1;

Normal_Temp=ImT36(upper_B:lower_B,center);
Normal_Temp_corr_ND=TempCorr_non_dielectric_AnodizedAL(0,ImT36(upper_B:lower_B,center),22.63);
Normal_Temp_corr_CS=TempCorr_Cosine(0,ImT36(upper_B:lower_B,center),22.63);


% take mean value of both "Normal_Temp" and "Normal_Temp_corr" to get a
% correction scaling factor: "Corr_Scale", which is based on their ratio

Mean_Normal_T=mean(Normal_Temp);
Mean_Normal_corr_T_ND=mean(Normal_Temp_corr_ND);
Corr_Scale_ND=Mean_Normal_T/Mean_Normal_corr_T_ND;

Mean_Normal_corr_T_CS=mean(Normal_Temp_corr_CS);
Corr_Scale_CS=Mean_Normal_T/Mean_Normal_corr_T_CS;



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
    ImT36_corr_ND_nS(:,i)=TempCorr_non_dielectric_AnodizedAL(Viewing_angle(:,i),ImT36(:,i),22.63);
    ImT36_corr_ND(:,i)=TempCorr_non_dielectric_AnodizedAL(Viewing_angle(:,i),ImT36(:,i),22.63)*Corr_Scale_ND;
    ImT36_corr_CS(:,i)=TempCorr_Cosine(Viewing_angle(:,i),ImT36(:,i),22.63)*Corr_Scale_CS;
    
end

%% compare the IR image before/after the viewing angle-based correction (non-dielectric model)
figure (7) % before applying the correction
subplot(1,2,1);
imagesc(ImT36);
caxis([30 33.5]);
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
axis off;

subplot(1,2,2); % after applying the correction
imagesc(ImT36_corr_ND);
caxis([30 33.5]);
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
axis off;
% take the mean temperature between the upper & lowe bound of ImT36
MeanT_Vertical_36=mean(ImT36(upper_B:lower_B,:));
MeanT_Vertical_36_corr_ND_ns=mean(ImT36_corr_ND_nS(upper_B:lower_B,:));
MeanT_Vertical_36_corr_ND=mean(ImT36_corr_ND(upper_B:lower_B,:));



% examine the 1D line along the curved surface befre/after correciton
figure (8);
subplot(1,2,1);
plot(theta(6:length(theta)-5),MeanT_Vertical_36(:,left_B+5:1:right_B-5),'b');
xlim([-70,70]);
ylim([25 37]);
grid on;
ylabel('Temperature (^oC)');
subplot(1,2,2);
plot(theta(6:length(theta)-5),MeanT_Vertical_36_corr_ND(:,left_B+5:1:right_B-5),'b');
hold on;
% also plot the corrected temperature before scaling
plot(theta(6:length(theta)-5),MeanT_Vertical_36_corr_ND_ns(:,left_B+5:1:right_B-5),'r-');
hold on;
xlim([-70,70]);
ylim([25 37]);

grid on;
ylabel('Temperature (^oC)');
%
%% compare the IR image before/after the viewing angle-based correction (Cosine function)
figure (9) % before applying the correction
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
imagesc(ImT36_corr_CS);
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
MeanT_Vertical_36_corr_CS=mean(ImT36_corr_CS(upper_B:lower_B,:));

% examine the 1D line along the curved surface befre/after correciton
figure (10);
subplot(1,2,1);
plot(theta(6:length(theta)-5),MeanT_Vertical_36(:,left_B+5:1:right_B-5),'b');
xlim([-70,70]);
ylim([25 35]);
grid on;
ylabel('Temperature (^oC)');
subplot(1,2,2);
plot(theta(6:length(theta)-5),MeanT_Vertical_36_corr_CS(:,left_B+5:1:right_B-5),'b');
xlim([-70,70]);
ylim([25 35]);
grid on;
ylabel('Temperature (^oC)');