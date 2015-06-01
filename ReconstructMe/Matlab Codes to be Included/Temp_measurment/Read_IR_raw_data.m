function [imR_ir]=Read_IR_raw_data(stamp,number)

imgnameR= strcat(stamp,'-image',num2str(number),'.ir.raw');
fid= fopen(imgnameR,'r','b');
header = fread(fid, [1 240], 'char');
dataR= fread(fid, [81920 1], 'uint16');
fclose(fid);
imR_ir=reshape(dataR, [320 256]);
imR_ir=imR_ir';