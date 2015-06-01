function[l_p,h_p,Lesion_B_IR_mapped,CIR_mapped,point_num]=lesion_points_mapping_auto(Lesion_B_W,imR_ir,CW,CIR,lowT,highT)

%% Transformation using homography estimation funciton:homography_ndlt
% Display the image for manual selection of corresponding point pairs

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
%imcontour(imR_ir,25); 
imshow(imR_ir,[lowT highT]),colormap('jet');
colorbar;
hold on;
scatter(CIR_mapped(1,:), CIR_mapped(2,:), 40, 'g', 'filled');
scatter(Lesion_B_IR_mapped(1,:), Lesion_B_IR_mapped(2,:), 2, 'k', 'filled');
title('please select the lesion point');

% get the centroid of lesion region
l_p_b= Lesion_B_IR_mapped(1:2,:)';
l_centroid= round(mean(l_p_b));
% get the probe region of lesion by selecting the 5x5 pixels mask centered
% at its centroid
[l_x,l_y]= meshgrid(l_centroid(1,1)-1:1:l_centroid(1,1)+1,l_centroid(1,2)-1:1:l_centroid(1,2)+1);
l_p(:,1)= reshape(l_x,9,1);
l_p(:,2)= reshape(l_y,9,1);
point_num=size(l_p,1);

%selection of healthy tissue coordinates
title('please select the "centroid" of healthy tissue point');
% determine the "translated displacement" from lesion to healthy point
h_centroid = ginput(1);
% shift the same point mask to "healthy site"
h_p(:,1)= l_p(:,1)+(h_centroid(:,1)-l_centroid(:,1));
h_p(:,2)= l_p(:,2)+(h_centroid(:,2)-l_centroid(:,2));


% show the selected "lesion and healthy" tissue point
scatter(l_p(:,1), l_p(:,2), 10, 'r');
scatter(h_p(:,1), h_p(:,2), 10, 'm','MarkerFaceColor',[1,1,1]);
colormap('jet');

%% put the row coordinate to the second index
% Lesion_B_temp(1,:)=Lesion_B_IR_mapped(2,:); 
% Lesion_B_temp(2,:)=Lesion_B_IR_mapped(1,:);
% Lesion_B_IR_mapped(1:2,:)=Lesion_B_temp;