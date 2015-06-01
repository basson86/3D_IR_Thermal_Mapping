function [ImT_corr]=TempConvert_NEW_Vig_Correct(ImR_i)

% From the calibration of Test 1(0222)
% m4 = 3.658757446064370e-12;
% m3 = -3.544860807308792e-07;
% m2 = 0.012874702415691;
% m1 = -2.077252604084555e+02;
% m0 = 1.256123903487680e+06;

% From the calibration of Test 2(0301)
m4 = -9.699380162015556e-13;
m3 = 9.470542369512357e-08;
m2 = -0.003469484795606;
m1 = 56.540702500396690;
m0 = -3.459387438171146e+05;

Temp_map=m0 + m1.*ImR_i  + m2 .* ImR_i.^2 + m3 .* ImR_i.^3 + m4 .*ImR_i.^4;


% load the 4 "pixel-based" polynomial coefficient arrays to correct
% Vignetting error
load ('0301_z1.mat');
load ('0301_z2.mat');
load ('0301_z3.mat');
load ('0301_z4.mat');

Fitted_Error= z1+ Temp_map.*z2 + (Temp_map.^2).*z3 + (Temp_map.^3).*z4 ;


ImT_corr= Temp_map + Fitted_Error;