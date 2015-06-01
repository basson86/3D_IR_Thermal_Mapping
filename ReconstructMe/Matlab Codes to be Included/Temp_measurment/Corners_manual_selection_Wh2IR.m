function [CW,CIR]=Corners_manual_selection_Wh2IR(WLfilename,imR_ir)

% read the white light image and resize it to 256 x 320
imW_original =im2double(imread(WLfilename,'JPEG'));
imW_resize = imresize(imW_original,[256 320]);
imW_resize = rgb2gray(imW_resize);

figure;
% White light image
imshow(imW_resize);
title('please select the points of 4 corners');
CW =ginput(4);
CW=CW';
% [BW, CWx, CWy] = roipoly(imW_resize);
% CW=[CWx';CWy'];
CW(3,:)=ones(1,length(CW));
close all

figure;
% Reference Frame of IR image
imagesc(imR_ir); 
title('please select the points of 4 corners');
CIR= ginput(4);
CIR=CIR';
CIR(3,:)=ones(1,length(CIR));
close all

