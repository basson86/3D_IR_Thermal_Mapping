% Read IR raw data and reshape it to 632x475 matrix
% by Tze-Yuan Cheng 2013/02/18
function [imR_ir]=Read_IR_raw_data_NEW(stamp,number)

imgnameR= strcat(stamp,'_',num2str(number,'%04.0f'),'.raw'); 
% the format of frame number is 4 digit integrs

fid= fopen(imgnameR,'r','b');
dataR= fread(fid, [300832 1], 'uint16');
fclose(fid);

% the dimension of the image is 476 x 632
imR_ir=reshape(dataR, [632 476]);
imR_ir=imR_ir';