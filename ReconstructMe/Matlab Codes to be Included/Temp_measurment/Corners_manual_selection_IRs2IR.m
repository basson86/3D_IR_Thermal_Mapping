function [CW,CIR]=Corners_manual_selection_IRs2IR(ImRS,imR_ir)

% read the steady-state IR image for corners mapping
imRs_gray =mat2gray(ImRS);

figure;
% White light image
imshow(imRs_gray);
colormap('jet');
title('please select the points of 4 corners');
CW =ginput(4);
CW=CW';
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

