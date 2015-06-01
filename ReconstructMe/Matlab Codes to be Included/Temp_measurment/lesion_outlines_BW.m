function[Lesion_B_W,imW_resize]=lesion_outlines_BW(WLfilename)

% read the white light image and resize it to 256 x 320
imW_original =im2double(imread(WLfilename,'JPEG'));
imW_resize = imresize(imW_original,[256 320]);
imW_resize = rgb2gray(imW_resize);

% Segment the lesion outlines using the funciton: "lesion_segmentation"
[Lesion_B_W_o]=lesion_segmentation(imW_resize);
Lesion_B_W_o=Lesion_B_W_o'
Lesion_B_W(1,:)= Lesion_B_W_o(2,:);
Lesion_B_W(2,:)= Lesion_B_W_o(1,:);
Lesion_B_W(3,:)=ones(1,length(Lesion_B_W));
