% Tze-Yuan Cheng 2013/02/21 
clear all; close all;clc;

%% The range of temperature for calibration
LowestTemp= 100; % (unit: 0.1 deg C) 
HighestTemp= 400;
% The increment temperature 
IncrementTemp = 005 ;% (unit: 0.1 deg C)  

%% The frame number setting for each image data
initial=0000; % the initial frame
final= 0049;  % the final frame
interval= 1;  % the interval number between 2 frames

%% 3D array to store the error(i,j,k) 
ImT_Volume = [];
T_Error_Volume = [];


%% Automatically read the IR image raw data at all temperature levels
for BBTemp= LowestTemp:IncrementTemp:HighestTemp
     BBTemp_Stamp= num2str(BBTemp);
     
     % Get frame-averaged image at each temperature level
     ImR_Avg_i=Average_multiple_frames(BBTemp_Stamp,initial,interval,final);
     
     % Conver the raw data to temeprature using the previous temperature
     % calibration polynomial coefficient
     ImT_Avg_i=TempConvert_NEW(ImR_Avg_i);
     
     % assemble all frame-averaed Temperature image to a 3D volume
     ImT_Volume = cat(3, ImT_Volume, ImT_Avg_i);
     
     % Compute coordinate-dependent eror: E(i,j,k) with respect to the controlled 
     % black body temperature 
     T_Error=BBTemp/10-ImT_Avg_i;
     
    
     % assemble all the error frame to a 3D volume
     T_Error_Volume = cat(3, T_Error_Volume, T_Error);          
end      

% save('0301_ImT_Volume.mat','ImT_Volume');
% save('0301_T_Error_toBB_Volume.mat','T_Error_Volume');


%% Applying the polynomial fit (polyfitn funciton) pixel by pixel 
% for the vignetting error over the temperature range K
load('0301_T_Error_toBB_Volume.mat');
load('0301_ImT_Volume.mat');

K= LowestTemp/10:IncrementTemp/10:HighestTemp/10;


for I=1:size(ImT_Volume,1);
    for J= 1:size(ImT_Volume,2);
            % pull out the sampled vignetting error from the calibration
            % data
            Error_at_IJ_all_T=reshape(T_Error_Volume(I,J,:),1,61);
            
            % pull out the "measured temperature" at pixel (i,j)
            T_measured=reshape(ImT_Volume(I,J,:),1,61);
            
            % polynomial fit w.r.t the measured temperature at pixel(i,j)
            p_IJ_Measured= polyfitn(T_measured,Error_at_IJ_all_T,3);
            
            % polynomial fit w.r.t the black body ground truth temperature K
            p_IJ_Temp= polyfitn(K,Error_at_IJ_all_T,3);
            
            % store the polynomial fit coefficients (error v.s measured
            % temperature) at pixel ij
            Coeff_m{I,J}=p_IJ_Measured.Coefficients; 
            
            % store the polynomial fit coefficients (error v.s "black body"
            % temperature) at pixel ij
            Coeff_bb{I,J}=p_IJ_Temp.Coefficients; 
            
    end
end

% save ('0301_Coeff_IJ_to_BBTemp.mat','Coeff_bb');
% 
% save ('0301_Coeff_IJ_to_MTemp.mat','Coeff_m');


%% Selecting 3 representative temparature images :10, 25, and 40 degc, 
%  test the correting effect by applying the polynomial obtained from the last step to the raw image
% load('0301_T_Error_toBB_Volume.mat');
% load('0301_ImT_Volume.mat');

% load the 2D "coeffiecent array"(each pixel has 4 polynomial coefficients in
% each element
load('0301_Coeff_IJ_to_MTemp.mat');

% select 3 tempreature data from the 3D array ImT_Volume
k1= 1; % 0 deg C
k2= 31; % 25 deg C
k3= 61; % 40 deg C

% get the calibration temperature image data
ImT1=ImT_Volume(:,:,k1);
ImT2=ImT_Volume(:,:,k2);
ImT3=ImT_Volume(:,:,k3);

% Get the 4 polynomial coefficient 2D arrays. z1...z4 
% ( each array has various coefficient at every pixel (I,J)) 
for I=1:size(Coeff_m,1);
    for J= 1:size(Coeff_m,2);
            z4(I,J)=Coeff_m{I,J}(1,1);
            z3(I,J)=Coeff_m{I,J}(1,2);
            z2(I,J)=Coeff_m{I,J}(1,3);
            z1(I,J)=Coeff_m{I,J}(1,4);
    end
end

% save ('0301_z1.mat','z1');
% save ('0301_z2.mat','z2');
% save ('0301_z3.mat','z3');
% save ('0301_z4.mat','z4');

% Based on the 4 coefficient arrays, compute the fitted error using 
% the "measured temperature : ImT1 - ImT3"
Fitted_Error_T1= z1+ ImT1.*z2 + (ImT1.^2).*z3 + (ImT1.^3).*z4 ;
Fitted_Error_T2= z1+ ImT2.*z2 + (ImT2.^2).*z3 + (ImT2.^3).*z4 ;
Fitted_Error_T3= z1+ ImT3.*z2 + (ImT3.^2).*z3 + (ImT3.^3).*z4 ;

% Pull-out the "true error" ( relative to black body temperature ) from the calibration data
True_Error_T1=T_Error_Volume(:,:,k1);
True_Error_T2=T_Error_Volume(:,:,k2);
True_Error_T3=T_Error_Volume(:,:,k3);

% apply the fitted error to the measured temperature image
% : T_corrteed = T_measured + Fitted_error

ImT1_corr= ImT1 + Fitted_Error_T1;
ImT2_corr= ImT2 + Fitted_Error_T2;
ImT3_corr= ImT3 + Fitted_Error_T3;

% plot the "fitted error" and "true error" together for comparison
figure;
subplot(1,2,1);
mesh(True_Error_T1);
title('measure errors from experiment');
subplot(1,2,2);
mesh(Fitted_Error_T1);
title('computed errors from polynomial fit');

%% Visualize the temperature image before/after applying the correction using "fitted-error"
% for the 3 selected temperature 

% plot the temperature image of 10 deg C
figure;
subplot(1,2,1);
imagesc(ImT1_corr,[min(min(ImT1)),max(max(ImT1))]);
colorbar;
subplot(1,2,2);
imagesc(ImT1);
colorbar;

% plot the temperature image of 25 deg C
figure;
subplot(1,2,1);
imagesc(ImT2_corr,[min(min(ImT2)),max(max(ImT2))]);
colorbar;
subplot(1,2,2);
imagesc(ImT2);
colorbar;

% plot the temperature image of 40 deg C
figure;
subplot(1,2,1);
imagesc(ImT3_corr,[min(min(ImT3)),max(max(ImT3))]);
colorbar;
subplot(1,2,2);
imagesc(ImT3);
colorbar;

%% testing the correction effect by converting any random image 
ImR_test=Read_IR_raw_data_NEW('315',10);

% Include the vignetting correction function into the temperature
% conversion function
ImT_test_corr=TempConvert_NEW_Vig_Correct(ImR_test);

% visualize the image for refernce
figure;
imagesc(ImT_test_corr);
colorbar;










%% ---------------------------------------------------------
% %% old code, reserved for future reference
% k=31;
% c1=Coeff_IJK{k}(1,1);
% c2=Coeff_IJK{k}(1,2);
% c3=Coeff_IJK{k}(1,3);
% c4=Coeff_IJK{k}(1,4);
% c5=Coeff_IJK{k}(1,5);
% c6=Coeff_IJK{k}(1,6);
% 
% for i=1:size(Coeff,1);
%     for j= 1:size(Coeff,1);       
%     ImT_corr(i,j)=c1*i^2 + c2*i*j + c3*i+ c4*j^2+ c5*j+ c6;       
%     end
% end
