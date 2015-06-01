function[Lesion_B_IR_mapped]=lesion_mapping_SS(Lesion_B_W,imR_ir,CW,CIR)

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
imagesc(imR_ir);
colorbar;
hold on;
scatter(CIR_mapped(1,:), CIR_mapped(2,:), 40, 'g', 'filled');
scatter(Lesion_B_IR_mapped(1,:), Lesion_B_IR_mapped(2,:), 8, 'r', 'filled');

