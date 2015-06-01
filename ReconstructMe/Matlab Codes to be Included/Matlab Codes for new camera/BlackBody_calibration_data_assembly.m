% For calibration purpose, read the black body image data which is named by its controlled
% temperature : Tze-Yuan Cheng 2013/02/21 
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

%% 3D array to store the image data at each temeprature
Temp_Sample_Volume_Whole = [];
Temp_Sample_Volume_Central = [];

%% 2D Array to store the spatial temeprature mean over the central portion of image
Intensity_vs_temp = [];

%% Automatically read the IR image raw data at all temperature levels
for Temp= LowestTemp:IncrementTemp:HighestTemp
    Temp_Stamp= num2str(Temp);
    
     %% Get a single 'averaged image' from the IR raw data sequence at each
     % temperature level
     ImR_Avg_i=Average_multiple_frames(Temp_Stamp,initial,interval,final);

     % Get the average intensity over the central portion: 60 x 60 portion of the image 
     ImR_Avg_Central_i = ImR_Avg_i(251:310,251:310);
        
     %% Append the image to the "Temperature Sample 3D array" (Whole image) 
     Temp_Sample_Volume_Whole = cat(3, Temp_Sample_Volume_Whole, ImR_Avg_i);

     % Append the image to the "Temperature Sample 3D array" (Central portion only)
     Temp_Sample_Volume_Central = cat(3,Temp_Sample_Volume_Central , ImR_Avg_Central_i);

     %% Get the spatial mean of tempeature over the central portion of image
     Spatial_mean_central_i= mean(mean(ImR_Avg_Central_i));
     Intensity_vs_temp= [Intensity_vs_temp;[Spatial_mean_central_i,Temp/10]];   

end

%% Result visualization

% plot the final image frame for reference
figure;
imagesc(Temp_Sample_Volume_Whole(:,:,2));

% plot the calibration curve based on spaitially & temperally-averaged data(intensity v.s bb temperature)
figure;
plot(Intensity_vs_temp(:,1),Intensity_vs_temp(:,2),'-*');
grid on;
xlabel('Raw IR Image Intensity');
ylabel('Controlled BB Temperature (degC)');
%% save the calibration data curve
%save('0301_BB_Calibration_curve_60by60.mat','Intensity_vs_temp');


%% polynomial fit of the calibration curve
%load('0301_BB_Calibration_curve_60by60.mat');

poly = polyfit(Intensity_vs_temp(:,1)',Intensity_vs_temp(:,2)',4);
poly

% residual error
Fitting_curve = polyval(poly,Intensity_vs_temp(:,1)');
Rerr= Intensity_vs_temp(:,2)'-Fitting_curve;

% display the fitting results
figure
subplot(2,1,1)
plot(Intensity_vs_temp(:,1),Intensity_vs_temp(:,2),'o',Intensity_vs_temp(:,1),Fitting_curve,'-')
grid on;
xlabel('Raw IR Image Intensity');
ylabel('Controlled BB Temperature (degC)');

% residual errors
subplot(2,1,2)
bar(Intensity_vs_temp(:,2),Rerr');
xlim([10,40])
grid on;
xlabel('Controlled BB Temperature (degC)');
ylabel('Residual Errors');


%% save the polynimial fitting results
%save ('0301_BB_Calibration_Polynimial.mat','poly')

