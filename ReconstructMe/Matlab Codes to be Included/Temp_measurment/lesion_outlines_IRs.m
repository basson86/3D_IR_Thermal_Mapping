function[Lesion_IRs,imRs_gray]=lesion_outlines_IRs(ImRS)

% read the steady state IR image 
imRs_gray = mat2gray(ImRS);

% Segment the lesion outlines using the funciton: "lesion_segmentation"
[Lesion_IRs_o]=lesion_segmentation(imRs_gray);
Lesion_IRs_o=Lesion_IRs_o'
Lesion_IRs(1,:)= Lesion_IRs_o(2,:);
Lesion_IRs(2,:)= Lesion_IRs_o(1,:);
Lesion_IRs(3,:)=ones(1,length(Lesion_IRs));

