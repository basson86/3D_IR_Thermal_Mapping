%% Start with clear workspace
clear all; close all; clc;
%% Load the IR image raw file
Data_path = 'C:\Users\Heat Transfer Lab\Desktop\IR_Camera_Black_Body_Calibration\data_0612';
path(path,Data_path);

SS_stamp='2014_06_12_15_08_47';
frame_num = 4;
[ImRS] = Read_IR_raw_data_NEW(SS_stamp,frame_num);
[ImTS] = TempConvert_NEW (ImRS);

str = num2str(frame_num);
imagesc(ImTS);
colormap('jet');
colorbar;
axis off;
h = text (10,25,str,'Color','w');


 %% Reading multiple images
% 
% for i = minnum:1:maxnum
%     [ImRS] = Read_IR_raw_data_NEW(SS_stamp,i);
%     [ImTS] = TempConvert_NEW (ImRS);
%     
% end