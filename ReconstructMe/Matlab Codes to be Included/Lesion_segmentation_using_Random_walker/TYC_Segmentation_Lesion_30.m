clear all
close all

img_r =im2double(imread('DSC00834','JPEG'));
img = imresize(img_r,[256 320]);

img = rgb2gray(img);
[X Y]=size(img);

%two seeds after resizing to 148 x 128
s1x=141; s1y=138; 
s2x=166; s2y=142; 



%two seeds after resizing to 256 x 320

%Apply the random walker algorithms
[mask,probabilities] = random_walker(img,[sub2ind([X Y],s1y,s1x), ...
    sub2ind([X Y],s2y,s2x)],[1,2]);

% segment the white pixel ( the region of interest)
E_white = regionprops(mask,'pixellist');
White_pxllist = E_white.PixelList;


%Display results
figure
imagesc(img);
colormap('gray')
axis equal
axis tight
axis off
hold on
plot(s1x,s1y,'g.','MarkerSize',10)
plot(s2x,s2y,'b.','MarkerSize',10)
title('Image with foreground (green) and background (blue) seeds')

figure
imagesc(mask)
colormap('gray')
axis equal
%axis tight
axis off
hold on
plot(s1x,s1y,'g.','MarkerSize',10)
plot(s2x,s2y,'b.','MarkerSize',10)
title('Output mask');

figure
[imgMasks,segOutline,imgMarkup]=segoutput(img,mask);
imagesc(imgMarkup);
colormap('gray')
axis equal
axis tight
axis off
hold on
title('Outlined mask');

%Get the coordinates of pixel segmented as lesion
[r, c] = find(segOutline==0);
coord_rc=cat(2,r,c);
