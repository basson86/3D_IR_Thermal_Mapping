% Start with clean Matlab Workspace
clear all; close all; clc

%% Image Case Information
% stamp='131830';  % practice trial 2012/02/16 -Hand_2_markers_round_1, gel-pack temperature
% initial= 69;
% % final= 20; 
% final= 70; 
% SSstamp='131453';  % steady state image
% SS_N= 1;
% time_step=1;
% Frame_int= 1;


stamp='135132';  % 2012/02/28 -Copper plate warm 0

SSstamp=stamp;  
initial= 39;
final= 150; 

SS_N= 1;
time_step=1;
Frame_int= 1;
lowT=8;
highT=37;

%% Temperature Conversion
% Convert the pixel intensity to temperature scale
[ImR1]= Read_IR_raw_data(stamp,initial);
[ImT1]= TempConvert(ImR1);

% Get the first movie frame
I= mat2gray(ImT1);

% steady state image for comparison 
[ImRS]= Read_IR_raw_data(SSstamp,SS_N);
[ImTS]= TempConvert(ImRS);

% Get the steady-state frame
I_S= mat2gray(ImTS);

%% Compute the mean value of the temperature in the steady-state image
mean_temp= median(median(ImTS));
Var_ROI= ImTS- mean_temp;

% display the temperature variation map in the steady state
figure;
imagesc(Var_ROI, [-1 max(max(Var_ROI))]);
colorbar;

figure;
%% Frame Reading loop
for i= initial:Frame_int:final
    % Get the a movie frame        
    [ImR_i]= Read_IR_raw_data(stamp,i);
    [ImT_i]= TempConvert(ImR_i);
    
    % Show the movie frame
    imagesc(ImT_i,[lowT highT]);  
    colorbar;
    title(strcat('frame No.=',num2str(i),' time=',num2str((i-initial+1)*time_step),'sec'));
    
% get frame for replay the movie    
F((i-initial)/Frame_int+1)=getframe(gcf);         
end
%% replay the movie
ht=figure;
movie(ht,F,30,5);
