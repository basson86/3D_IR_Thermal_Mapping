function[Cut_line,P1,P2,Lesion_B_IR_mapped,CIR_mapped]=lesion_points_mapping_manual_line(Lesion_B_W,imR_ir,CW,CIR,lowT,highT)

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
imagesc(imR_ir,[lowT highT]),colormap('jet');
colorbar;
hold on;
scatter(CIR_mapped(1,:), CIR_mapped(2,:), 40, 'g', 'filled');
scatter(Lesion_B_IR_mapped(1,:), Lesion_B_IR_mapped(2,:), 10, 'r', 'filled');
title('please select the first point P1 of the cut line');

%% get the "first point: P1" of the cut line
P1= ginput(1);

P1_x=round(P1(1,1));
P1_y=round(P1(1,2));
    

%% get the "second point: P2" of the cut line
title('please select the second point of the cut line');

P2= ginput(1);

P2_x=round(P2(1,1));
P2_y=round(P2(1,2));

%% Based on the line determined by P1 and P2, spread a "directional points" along that line 

Cut_line= [];

slop= (P2_y-P1_y) / (P2_x-P1_x);

for x= P1_x:1:P2_x

    y= (slop*(x-P1_x)) + P1_y;
    
    % if y is integer, put them into the Cut_line coordinate array
    %if (y - round(y))==0;
        
    Cut_line=[Cut_line;[x,y]]; 
        
    %end
end   
    
% show the selected "lesion and healthy" tissue point
scatter(Cut_line(:,1), Cut_line(:,2), 10, 'r');
colormap('jet');

