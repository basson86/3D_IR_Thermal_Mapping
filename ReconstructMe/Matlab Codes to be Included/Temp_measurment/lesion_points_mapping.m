function[l_p,h_p,Lesion_B_IR_mapped,CIR_mapped]=lesion_points_mapping(Lesion_B_W,imR_ir,point_num)

%% Transformation using homography estimation funciton:homography_ndlt
% Display the image for manual selection of corresponding point pairs

%[CW,CIR]=Corners_manual_selection(imW_resize,imR_ir);

load ('CW.mat');
load ('CIR.mat');


%% Obtain the homography H from X1 & X2 using homography_ndlt
[H_W2IR]=homography_ndlt(CW,CIR);

% Transform the corners' points & lesion boundary location to IR using
% "H_W2IR"

CIR_mapped= H_W2IR * CW;
Lesion_B_IR_mapped= H_W2IR * Lesion_B_W;

for i=1:length(CIR_mapped)
CIR_mapped(:,i)=CIR_mapped(:,i)/CIR_mapped(3,i);
end

for i=1:length(Lesion_B_IR_mapped)
Lesion_B_IR_mapped(:,i)=Lesion_B_IR_mapped(:,i)/Lesion_B_IR_mapped(3,i);
end

% Display the mapping results of H_W2IR
figure;
imcontour(imR_ir,25); 
colorbar;
hold on;
scatter(CIR_mapped(1,:), CIR_mapped(2,:), 40, 'g', 'filled');
scatter(Lesion_B_IR_mapped(1,:), Lesion_B_IR_mapped(2,:), 8, 'r', 'filled');
title('please select the lesion point');
l_p= ginput(point_num);
title('please select the healthy tissue point');
h_p= ginput(point_num);