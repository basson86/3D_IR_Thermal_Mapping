% Read multiple frames for a single black body temperature level 
% by Tze-Yuan Cheng 2013/02/18

function [ImR_Avg]=Average_multiple_frames(time_stamp,initial,interval,final)

%Number_of_frame=(final-initial+1)/interval; % total number of frames to be read

Number_of_frame = 0;

ImR_Sum= zeros(476,632);% Initialize the total sum of all images frames for averaging

% reading loop
for i=initial:interval:final

ImR_i= Read_IR_raw_data_NEW(time_stamp,i);

ImR_Sum = ImR_Sum + ImR_i;

Number_of_frame=Number_of_frame+1; 
end

% Compute the resulting 'averaged IR fram':
ImR_Avg= ImR_Sum./Number_of_frame;





