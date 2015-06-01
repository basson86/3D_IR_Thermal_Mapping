
clear all;close all;

%% path for new IR camera file processing
p= 'C:\Users\Heat Transfer Lab\Documents\MATLAB\kinect_calibration_with_data_2_1\Our Experiment\Matlab Codes for new camera';
path(path,p);

p= 'C:\Users\Heat Transfer Lab\Documents\MATLAB\kinect_calibration_with_data_2_1\Our Experiment\Data_0814\';
path(path,p);

%%

for i=0:47

stamp=num2str(i,'%04.0f');
number=0;
[imR_ir]=Read_IR_raw_data_NEW(stamp,number);

imR_gray=mat2gray(imR_ir);

imwrite(imR_gray,strcat(stamp,'_IR.jpg'),'jpg');

end