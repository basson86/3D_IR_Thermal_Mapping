function [coord_rc]=lesion_segmentation(img);

%select 2 seed to apply the random walker segmentation
[X Y]=size(img);
imshow(img);
hold on;
title('select the 2 points, one in lesion region, the other in normal skin');
[sx,sy] = ginput(2);
close all

sx=round(sx);
sy=round(sy);
s1x= sx(1);s2x= sx(2);
s1y= sy(1);s2y= sy(2);


%two seeds after resizing to 256 x 320

%Apply the random walker algorithms
[mask,probabilities] = random_walker(img,[sub2ind([X Y],s1y,s1x), ...
sub2ind([X Y],s2y,s2x)],[1,2]);

% segment the white pixel ( the region of interest)
E_white = regionprops(mask,'pixellist');
White_pxllist = E_white.PixelList;

figure;
[imgMasks,segOutline,imgMarkup]=segoutput(img,mask);
imagesc(imgMarkup);
colormap('gray');
axis equal
axis tight
axis off
hold on
title('Outlined mask');

%Get the coordinates of pixel segmented as lesion
[r, c] = find(segOutline==0);
coord_rc=cat(2,r,c);



