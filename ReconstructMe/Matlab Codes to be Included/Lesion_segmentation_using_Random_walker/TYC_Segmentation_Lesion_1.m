clear
close all

img_r =im2double(imread('img_0010','JPEG'));
img= imcrop(img_r,[1030, 1390, 2180-1030, 2550-1390]);
img = imresize(img,[148 128]);
%img = imresize(img,[256 320]);

img = rgb2gray(img);
[X Y]=size(img);

%two seeds after resizing to 148 x 128
s1x=69; s1y=75; 
s2x=46; s2y=83; 



%two seeds after resizing to 256 x 320
%s1x=171; s1y=127; %Note that seed location is not central to object
%s2x=132; s2y=158; 

%Set two seeds
%s1x=472; s1y=572; %Note that seed location is not central to object
%s2x=350; s2y=671; 

% s1x=195; s1y=152; %Note that seed location is not central to object
% s2x=67; s2y=192; 

%Apply the random walker algorithms
[mask,probabilities] = random_walker(img,[sub2ind([X Y],s1y,s1x), ...
    sub2ind([X Y],s2y,s2x)],[1,2]);


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
%plot(s1x,s1y,'g.','MarkerSize',10)
%plot(s2x,s2y,'b.','MarkerSize',10)
title('Outlined mask');

%Get the coordinates of pixel segmented as lesion
[r, c] = find(segOutline==0);
coord_rc=cat(2,r,c);


% figure
% imagesc(probabilities(:,:,1))
% colormap('gray')
% axis equal
% axis tight
% axis off
% hold on
% plot(s1x,s1y,'g.','MarkerSize',10)
% plot(s2x,s2y,'b.','MarkerSize',10)
% title(strcat('Probability at each pixel that a random walker released ', ...
%     'from that pixel reaches the foreground seed'));