%% 05/20 by Tze-Yuan Cheng: 3D image point cloud generation using Kinect - the feasibility verification code
% Given a 2D image of face and its 3D surface profile data, generate the 3D image :
% This code integrate the data acquired by "ReconstructMe" with the calibration results from the literature 'Joint depth and color camera calibration with
% distortion correction'. We acquired the 3D reconstruction data of the face using kinect, and then acquire the image pair using kinect's two embedded camera 
% : depth camera and rgb camera. Since the transformation matrix from depth camera to rgb camera is already calibrated by the literature :'Joint depth and color camera calibration with
% distortion correction', we could just apply the external parameters provided by the code to realize the 3D imaging using the gray scale version of the 2D image. 


clear all; close all;
%% Toolbox codes loading 

% the toolbox for loading the .ply files of ReconstructMe in matlab 
% The toolbox can be downloaded from Matlab Central http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5355&objectType=FILE

%add dapa path: 'C:\Users\Heat Transfer Lab\Desktop\ReconstructMe\toolbox_graph\toolbox_graph';
pdata = strcat(pwd,'\toolbox_graph');
path(path, pdata);
%add dapa path: 'C:\Users\Heat Transfer Lab\Desktop\ReconstructMe\toolbox_graph\toolbox_graph\toolbox';
pdata = strcat(pwd,'\toolbox_graph\toolbox');
path(path, pdata);

% the toolbox for Kinect Calibration 
%add dapa path: 'C:\Users\Heat Transfer Lab\Documents\MATLAB\kinect_calibration_with_data_2_1\v2.1\Mydata';
pdata=strcat(pwd,'\kinect_calibration_with_data_2_1\v2.1\Mydata');
path(path, pdata);

%add dapa path: 'C:\Users\Heat Transfer Lab\Documents\MATLAB\kinect_calibration_with_data_2_1\v2.1\toolbox';
pdata=strcat(pwd,'\kinect_calibration_with_data_2_1\v2.1\toolbox');
path(path, pdata);


%% load the mesh
clear options;
filename='TYC_0521.ply';

[vertex,faces] = read_ply(filename);
options.name = filename; % useful for displaying


%% display the mesh

clf;
face_seg= faces(1:53000,:);
% figure(1)
% trisurf ( face_seg, vertex(:,1), vertex(:,2), vertex(:,3) );
% axis equal;
% shading interp;
% caxis([-750 -600])


%% You can zoom on the mesh and display its triangulated faces
% figure(2);
% clf;
% for i=1:4
%     subplot(2,2,i);
%     plot_mesh(vertex, face_seg);
%     shading faceted;
%     zoom(1.4^(i+1));
% end

%% compute the curvature & normal vector
options.curvature_smoothing = 10;
options.verb = 0;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,faces,options);



%%
figure(3);
clf;
% options.face_vertex_color = perform_saturation(Cgauss,1.2);
% plot_mesh(vertex,face_seg, options); shading interp; colormap jet(256);
% title('Gaussian curvature');
options.face_vertex_color = perform_saturation(abs(Cmin)+abs(Cmax),1.2);
plot_mesh(vertex,face_seg, options); shading interp; %colormap jet(256);
axis tight;
%title('Total curvature');
%ylim([-1000 100]);


%% Load the calibration results (from the path :'C:\Users\Heat Transfer Lab\Documents\MATLAB\kinect_calibration_with_data_2_1\')

do_load_calib(strcat(pwd,'\Camera 3D Calibration Parameters\small_set.mat'));
load(strcat(pwd,'\Camera 3D Calibration Parameters\Data_0814\Calibration 0822\final_results_0822_4_WD'));

%do_load_calib('C:\Users\Heat Transfer Lab\Desktop\ReconstructMe\20140103 All Figures Shifting according to TC1218 (MAIN)\Camera 3D Calibration Parameters\small_set.mat');
%load('C:\Users\Heat Transfer Lab\Desktop\ReconstructMe\20140103 All Figures Shifting according to TC1218 (MAIN)\Camera 3D Calibration Parameters\Data_0814\Calibration 0822\final_results_0822_4_WD');



%% read disparity image and C1 image

imd=read_disparity('TYC_0521.pgm');
imc1= imread('TYC_0521.jpg');

% im_size [2] : size of the resulting depthmap in c1 coordinate (using the 
% image size of c1: 1024 x 1280)

imsizeC1=[960,1280];

% splat_size [1] size of the splat kernel for each warped depth point
% (must be odd, default 1)
splat_size=3;

% conver the pixel in depthmap to the coordinate of rgb image (c1)
[depthmap,points_DtoC1]=compute_rgb_depthmap(imd,final_calib,imsizeC1,splat_size);

%% overlay the projected pixel from depth map to the rgb image coordinate, to see if the calbiration results is accurate

% 
% figure(4);
% imshow(points_DtoC1);
% hold on;
% h=imshow(imc1(:,95:1280)); % since the official calibration results was not made by our Kinect camera, the mapping is not totally accurate, 
% set(h, 'AlphaData', 0.5);  % so we made some manual correction for testing purpose
% 
% figure;
% imagesc(imd);

%% Extract the parameters from calibration for 3D to 2D transformation
% External parameters Tdc : camera coordinate transfomration from {D} to {C1}
Rdc=final_calib.dR;
tdc=final_calib.dt;
Tdc= [[Rdc;[0,0,0]],[tdc;1]];


%% internal parameters of C1 

% scalling factor and image center
fxy_uv0=final_calib.rK{1,1};

% distortion parameters
kc1=final_calib.rkc{1,1};
k1=kc1(1);
k2=kc1(2);
k3=kc1(3);
k4=kc1(4);
k5=kc1(5);
%% transformation from 3D camera coordinate in 3D --> 2D image coorinate in C1

% the segmented face profile in D camera coordinate

Xd_raw=vertex(1:28500,:)';

Xd=Xd_raw;
Xd(3,:)=-Xd(3,:);
Xd(4,:)=1;


%% the transformed 3D coordinate Xc in C1
Xc=Tdc*Xd;


% figure(5);
% subplot(1,2,1);
% scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),'.');
% axis equal;
% zlim([-1300 -1100]);
% title('Depth Camera Coordinate');
% 
% 
% subplot(1,2,2);
% scatter3(Xc(1,:),Xc(2,:),Xc(3,:),'.');
% axis equal;
% view(90,0);
% zlim([1100 1300]);
% title('RGB Camera Coordinate')

%%
for i=1:length (Xc);
xn(1:2,i)=Xc(1:2,i)/Xc(3,i);
xn(3,i)=1;
end

% image distortion computation
% (check page 3 of the paper : 'Joint depth and color camera calibration with distortion correction'
r_square= xn(1,:).^2+xn(2,:).^2;
xg=[2*k3*xn(1,:).*xn(2,:)+k4*(r_square+2*xn(1,:).^2);...
    2*k4*xn(1,:).*xn(2,:)+k3*(r_square+2*xn(2,:).^2);];

for i=1:length(xn)
xk(1:2,i)=(1+k1*r_square(1,i)+k2*r_square(1,i)^2+k5*r_square(1,i)^3)*xn(1:2,i) +xg(1:2,i);
xk(3,i)=1;
end

%% compute the pixel correspondence in C1 image
uv_C1= fxy_uv0*xk; 

for i=1: length(uv_C1)
    if (uv_C1(1,i)>1280 ||uv_C1(2,i)>1024);
    uv_C1(1:2,i)=[0;0];
    end
end

%% manual final correction to get the best image overlay between the "projected points" and "rgb image" 
uv_C1_m(1,:)=uv_C1(1,:)-170;
uv_C1_m(2,:)=(960.-uv_C1(2,:))-110;  % "960 - y" is converting the pixel coordinate to MATLAB fashision ( upper-left corner is the origin)  

% show the coordinate differences between "projected points" and "rgb image " 
% figure(6);
% subplot(1,2,1);
% scatter(uv_C1(1,:),uv_C1(2,:),'b.');
% subplot(1,2,2);
% imagesc(imc1); 
% 
% 
% % show the resulting image overlay after manual correction
% figure(7);
% h=imshow(imc1); 
% set(h, 'AlphaData', 0.8);
% hold on;
% scatter(uv_C1_m(1,:),uv_C1_m(2,:),'b.');
%% 2D interpolation in RGB 2D image
% Im_pro=[];
% imc1_gray=rgb2gray(imc1);
% for i=1:length(uv_C1_m);
% Im_pro(1,i)= interp2(imc1_gray,uv_C1_m(1,i),uv_C1_m(2,i));
% end

load('TYC_0521_Im_pro_C1.mat');

%% color the 3D point cloud using the interpolated pixel value in the depth camera coordinate
% figure(8);
% scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),200,Im_pro,'.');
% axis equal;
% zlim([-1300 -1100]);
% title('Depth Camera Coordinate');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------C1 to IR Computation-----------------------------------



%% Convert the projected point in {C1} to {IR} using the calibrated parameter given in
%"C:\Users\Heat Transfer Lab\Documents\MATLAB\kinect_calibration_with_data_2_1\Our Experiment\Data_0814\Calibration 0822\final_results_0822_4_WD"
%load('C:\Users\Heat Transfer Lab\Documents\MATLAB\kinect_calibration_with_data_2_1\Our Experiment\Data_0814\Calibration 0822\final_results_0822_4_WD');

% load the external transformation matrix "Tir_c"

Rir_c=final_calib.rR{1,2};
tir_c=final_calib.rt{1,2};
Tir_c= [[Rir_c;[0,0,0]],[tir_c;1]];

%transformation computation
Xir=Tir_c*Xc;


% Rir_c= [0.95805,-0.08309,0.27431;...
%         0.04680,0.98956,0.13628;...
%         -0.28277,-0.11772,0.95194;...
%         0,0,0];
% 
% tir_c=[0.01006,0.04877,0.06105,1]';    

figure(9);
subplot(1,2,1);
scatter3(Xc(1,:),Xc(2,:),Xc(3,:),'.');
axis equal;
view(90,0);
zlim([1100 1300]);
title('RGB Camera Coordinate')

subplot(1,2,2);

scatter3(Xir(1,:),Xir(2,:),Xir(3,:),'.');
axis equal;
view(90,0);
zlim([1050 1300]);
title('IR Camera Coordinate');



% Normalize with the thir compoent in Xc(3,i)
for i=1:length (Xir);
xn_ir(1:2,i)=Xir(1:2,i)/Xir(3,i);
xn_ir(3,i)=1;
end

% internal parameters 
% scalling factor and image center
fxy_uv0_ir=final_calib.rK{1,2};

% distortion parameters
kc2=final_calib.rkc{1,2};
k_ir1=kc2(1);
k_ir2=kc2(2);
k_ir3=kc2(3);
k_ir4=kc2(4);
k_ir5=kc2(5);


% compute the pixel correspondence in C1 image
uv_IR= fxy_uv0_ir*xn_ir; 

% for i=1: length(uv_IR)
%     if (uv_IR(1,i)>632 ||uv_IR(2,i)>476);
%     uv_IR(1:2,i)=[0;0];
%     end
% end


%% read IR image
stamp='TYC_0521';
[ImR]= Read_IR_raw_data_NEW(stamp,0);
[ImT]= TempConvert_NEW_Vig_Correct(ImR);

Temp_limit=[28.3 33.8];

%% shift all the tempearture value 3.3 degC higher to fit the measured temeprature by thermal couple
ImT=ImT+3.3;


% modify the projected point in IR according to image coordinate
uv_IR_m(1,:)=uv_IR(1,:)-225;
uv_IR_m(2,:)=(240.-uv_IR(2,:));  % "240 - y" is converting the pixel coordinate to MATLAB fashision ( upper-left corner is the origin)


% show the coordinate differences between "projected points" and "rgb image " 
figure(10);
subplot(1,2,1);
scatter(uv_IR(1,:),uv_IR(2,:),'g.');
subplot(1,2,2);
imagesc(ImT); 


% show the resulting image overlay after manual correction
figure(11);
h=imagesc(ImT); 
set(h, 'AlphaData', 0.8);
hold on;
scatter(uv_IR_m(1,:),uv_IR_m(2,:),'g.');
axis off;

%% 2D interpolation in IR image

Im_pro_IR=[];

% for i=1:length(uv_IR_m(1,:));
% Im_pro_IR(1,i)= interp2(ImT,uv_IR_m(1,i),uv_IR_m(2,i)) ;
% end
% save('TYC_0823_Im_pro_IR_SHIFT_TC.mat','Im_pro_IR');

% load('TYC_0823_Im_pro_IR.mat');
load('TYC_0823_Im_pro_IR_SHIFT_TC.mat');

%% color the 3D point cloud using the interpolated IR image pixel value 
figure(12);
scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),190,Im_pro_IR,'.');
axis equal;
zlim([-1300 -1100]);
title('IR pixel mapping in 3D space');
colormap('jet');
caxis(Temp_limit);
view(150,110);
colorbar;

%% store the "normal vector" for future viewing angle computation
options.normal = Normal;
% figure(13);
% plot_mesh(vertex,face_seg, options); shading interp; colormap jet(256);
% axis tight;

%% compute the viewing angle based on the normal vector of each vertex
Xd_normal=Normal(:,1:length(Xd_raw));

Im_pro_IR_corr=[];
for i= 1: length(Xd_normal)
Xd_angle(:,i) = acosd(dot(Xd_normal(:,i),[0,0,1])/norm(Xd_normal(:,i),2));
  
  % set the upperbound : 90 degree as the viewing angle cap
  if (Xd_angle(:,i)>90)
  Xd_angle(:,i)=90;
  end


% apply the viewing angle correction on the measured IR temperature
% and get a "corrected" interpolated pixel mapping : Im_pro_IR_corr

[Tc]= deltaT_by_skin_viewing_angle_ND(Xd_angle(:,i),Im_pro_IR(1,i),22);

Im_pro_IR_corr(:,i)= Tc; 

end 


%% show the viewing angle distribution over the skin surface
figure(13);
scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),190,Xd_angle,'.');
axis equal;
zlim([-1300 -1100]);
title('viewing angle mapping in 3D space (Max = 90 deg)  (unit:degree)');
colormap('jet');
caxis([0 90]);
view(150,110);
colorbar;


%% 3D to 2D back-interpolation 
[x_IR,y_IR] = meshgrid(1:632,1:476);


%% back project the viewing angle map to 2D image, for smoothing purpose
Angle_IR = griddata(uv_IR_m(1,:),uv_IR_m(2,:),Xd_angle,x_IR,y_IR);

% apply 2D "median filter" on 2D backprojected viewing angle map 
% to reduce the noise of viewing angle computed in 3D space
% return: 476x632 2D map: "Angle_IR_smooth"

Angle_IR_smooth = medfilt2(Angle_IR,[4,4]);


%% interpolate w.r.t "Angle_IR_smooth" processed in 2D image
% and project it back to 3D space
% return: 3xn vector : "Xd_angle_smooth"
% Xd_angle_smooth=[];
% 
% for i=1:length(uv_IR_m(1,:));    
% Xd_angle_smooth(1,i)= interp2(Angle_IR_smooth,uv_IR_m(1,i),uv_IR_m(2,i)) ;
% end
% save('TYC_1104_Xd_angle_smooth.mat','Xd_angle_smooth');

load('TYC_1104_Xd_angle_smooth.mat');

%Xd_angle_smooth=Xd_angle;

%% re-plot the filtered viewing angle distribution in 3D
figure(14);
scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),190,Xd_angle_smooth,'.');
axis equal;
zlim([-1300 -1100]);
title('Smoothed viewing angle mapping in 3D space (unit:degree)');
colormap('jet');
caxis([0 90]);
view(150,110);
colorbar;

%% apply the viewing angle correction AGAIN based on the smoothed viewing
% angle

% Select the Ambinet Temperature for correction formula 
Tamb=23; %unit:(degC)

for i= 1: length(Xd_angle_smooth)
[Tc_ND]= deltaT_by_skin_viewing_angle_ND(Xd_angle_smooth(:,i),Im_pro_IR(1,i),Tamb);
[Tc_D]= deltaT_by_skin_viewing_angle_D(Xd_angle_smooth(:,i),Im_pro_IR(1,i),Tamb);
[Tc_Cosine]= deltaT_by_skin_viewing_angle_Cosine(Xd_angle_smooth(:,i),Im_pro_IR(1,i),Tamb);

Im_pro_IR_corr_smooth_D(:,i)= Tc_D;
Im_pro_IR_corr_smooth_ND(:,i)= Tc_ND;
Im_pro_IR_corr_smooth_Cosine(:,i)= Tc_Cosine;

end 


%% IR thermography back projection
Temp_IR_uncorr = griddata(uv_IR_m(1,:),uv_IR_m(2,:),Im_pro_IR,x_IR,y_IR);
Temp_IR_corr_smooth_D = griddata(uv_IR_m(1,:),uv_IR_m(2,:),Im_pro_IR_corr_smooth_D,x_IR,y_IR);
Temp_IR_corr_smooth_ND = griddata(uv_IR_m(1,:),uv_IR_m(2,:),Im_pro_IR_corr_smooth_ND,x_IR,y_IR);
Temp_IR_corr_smooth_Cosine = griddata(uv_IR_m(1,:),uv_IR_m(2,:),Im_pro_IR_corr_smooth_Cosine,x_IR,y_IR);

% take a sampled point on the face with small viwing angle on both
% "Temp_IR_uncorr" & "Temp_IR_corr", and take their ratio as the scalling
% factor, in this case, the image at [y=219,x=312] is used
scaling_factor_D = Temp_IR_uncorr(219,312)/Temp_IR_corr_smooth_D(219,312);
scaling_factor_ND = Temp_IR_uncorr(219,312)/Temp_IR_corr_smooth_ND(219,312);
scaling_factor_Cosine = Temp_IR_uncorr(219,312)/Temp_IR_corr_smooth_Cosine(219,312);




%% color the 3D point cloud using the "2-steps Viewing-angle Corrected" interpolated IR image pixel value 

Temp_limit=[31.5 33.8];

% before correction
figure(15);
scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),190,Im_pro_IR,'.');
axis equal;
zlim([-1300 -1100]);
title('IR pixel mapping in 3D space (before correction)');
colormap('jet');
%caxis(Temp_limit);
caxis([28 34]);

view(150,110);
colorbar;

%% after smoothing ( 90 degree upper bound +  2D median filter)
figure(16);
scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),190,scaling_factor_D*Im_pro_IR_corr_smooth_D,'.');
axis equal;
zlim([-1300 -1100]);
title('IR pixel mapping in 3D space (Dielectric model correction)');
colormap('jet');
%caxis(Temp_limit);
caxis([28 34]);
view(150,110);
colorbar;
%%

figure(17);
scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),190,scaling_factor_ND*Im_pro_IR_corr_smooth_ND,'.');
axis equal;
zlim([-1300 -1100]);
title('IR pixel mapping in 3D space (Non-dielectric model correction)');
colormap('jet');
caxis(Temp_limit);
view(150,110);
colorbar;
%%
figure(18);
scatter3(Xd_raw(1,:),Xd_raw(2,:),Xd_raw(3,:),190,scaling_factor_Cosine*Im_pro_IR_corr_smooth_Cosine,'.');
axis equal;
zlim([-1300 -1100]);
title('IR pixel mapping in 3D space (Cosine model correction)');
colormap('jet');
caxis(Temp_limit);
view(150,110);
colorbar;

%% plot the back-projected 3D point to 2D image frame
%%s_Temp_IR_corr_smooth= scaling_factor * Temp_IR_corr_smooth;

figure(19);
imagesc(Angle_IR_smooth);
caxis([0 90]);
colorbar;
%% select the skin range to visualize the temperature profile using different emissivity model
% the range of skin in IR image
skin_range_R=80:320;
skin_range_C=222:415;
% the selected row to plot 1D tmperature profile
choosen_row1=50;
choosen_row2=110;
choosen_row3=135;
% the selected range of column 
choosen_C_range=30:150;
Start=30;
Middle=89;
End=149;

Temp_limit_back2D=[31.3 33.8];

crop_Temp_uncorr=Temp_IR_uncorr(skin_range_R,skin_range_C);
crop_Temp_corr_D=Temp_IR_corr_smooth_D(skin_range_R,skin_range_C);
crop_Temp_corr_ND=Temp_IR_corr_smooth_ND(skin_range_R,skin_range_C);
crop_Temp_corr_Cosine=Temp_IR_corr_smooth_Cosine(skin_range_R,skin_range_C);

%% Display the 3D to 2D back-projection thermography before/after correction

% noise filetering: remove the overshot singular value from the two image :  
% "crop_Temp_corr_D" and "crop_Temp_corr_ND"

Over_shot = 60; % choose the over-shot tempeature threshold for processing (deg C)
I_D = mat2gray(crop_Temp_corr_D, [0 Over_shot]); % convert matdata to grey level based on threshild
I_ND = mat2gray(crop_Temp_corr_ND, [0 Over_shot]);

% make the binary mask using the grey_level
BW_mask_D = im2bw(I_D,0.9);  %
BW_mask_ND = im2bw(I_ND,0.9);

% apply the inverted binary mask "(1-BW)" to filter the original corrected image
% and get the filtered image:"crop_Temp_corr_D_f" and "crop_Temp_corr_ND_f"
crop_Temp_corr_D_f=(1-BW_mask_D).*crop_Temp_corr_D;
crop_Temp_corr_ND_f=(1-BW_mask_D).*crop_Temp_corr_ND;

% 
% the uncorrected case
figure(20);
imagesc(crop_Temp_uncorr);
hold on;
plot(choosen_C_range,choosen_row1,'b.','MarkerSize',5);
plot(choosen_C_range,choosen_row2,'g.','MarkerSize',5);
scatter(Start,choosen_row1,80,'fill','k');
scatter(Middle,choosen_row1,80,'fill','k');
scatter(End,choosen_row1,80,'fill','k');
scatter(Start,choosen_row2,80,'fill','k');
scatter(Middle,choosen_row2,80,'fill','k');
scatter(End,choosen_row2,80,'fill','k');
%plot(choosen_C_range,choosen_row3,'k.','MarkerSize',10);
caxis(Temp_limit_back2D);
title('uncorrected IR image');
axis off;
colorbar;
%


% the corrected case using dielectric model (filtered out the over shot noise)
figure(21);
imagesc(crop_Temp_corr_D_f);
hold on;
plot(choosen_C_range,choosen_row1,'b.','MarkerSize',5);
plot(choosen_C_range,choosen_row2,'g.','MarkerSize',5);
scatter(Start,choosen_row1,80,'fill','k');
scatter(Middle,choosen_row1,80,'fill','k');
scatter(End,choosen_row1,80,'fill','k');
scatter(Start,choosen_row2,80,'fill','k');
scatter(Middle,choosen_row2,80,'fill','k');
scatter(End,choosen_row2,80,'fill','k');
%plot(choosen_C_range,choosen_row3,'k.','MarkerSize',4);
caxis(Temp_limit_back2D);
title('correction by dielectric model');
axis off;
colorbar;

% the corrected case using non-dielectric model (filtered out the over shot noise)
figure(22);
imagesc(crop_Temp_corr_ND_f);
hold on;
plot(choosen_C_range,choosen_row1,'b.','MarkerSize',5);
plot(choosen_C_range,choosen_row2,'g.','MarkerSize',5);
%plot(choosen_C_range,choosen_row3,'k.','MarkerSize',4);
scatter(Start,choosen_row1,80,'fill','k');
scatter(Middle,choosen_row1,80,'fill','k');
scatter(End,choosen_row1,80,'fill','k');
scatter(Start,choosen_row2,80,'fill','k');
scatter(Middle,choosen_row2,80,'fill','k');
scatter(End,choosen_row2,80,'fill','k');
caxis(Temp_limit_back2D);
title('correction by non-dielectric model');
axis off;
colorbar;


% figure(23);
% imagesc(crop_Temp_corr_Cosine);
% hold on;
% plot(choosen_C_range,choosen_row1,'b.','MarkerSize',4);
% plot(choosen_C_range,choosen_row2,'g.','MarkerSize',4);
% plot(choosen_C_range,choosen_row3,'k.','MarkerSize',4);
% caxis([25 30.5]);
% title('IR pixel mapping in 2D space (Cosine model)');
% colorbar;

%% plot the 1D profile alone the chosen row

temp_limit=[31.5 35];


figure(24);
%subplot(2,2,1);
plot(crop_Temp_uncorr(choosen_row1,choosen_C_range),'b');
hold on;
plot(crop_Temp_uncorr(choosen_row2,choosen_C_range),'g--');
hold on;
plot(Start-Start+1,crop_Temp_uncorr(choosen_row1,Start),'ko','Markersize',10);
plot(Middle-Start+1,crop_Temp_uncorr(choosen_row1,Middle),'ko','Markersize',10);
plot(End-Start+1,crop_Temp_uncorr(choosen_row1,End),'ko','Markersize',10);
plot(Start-Start+1,crop_Temp_uncorr(choosen_row2,Start),'ko','Markersize',10);
plot(Middle-Start+1,crop_Temp_uncorr(choosen_row2,Middle),'ko','Markersize',10);
plot(End-Start+1,crop_Temp_uncorr(choosen_row2,End),'ko','Markersize',10);
% plot(crop_Temp_uncorr(choosen_row3,choosen_C_range),'k');
grid on;
ylim(temp_limit);
xlim([0 120]);
ylabel('Temperature ( ^oC)');
title('uncorrected 3D thermography');

figure(25);
%subplot(2,2,2);
plot(crop_Temp_corr_D(choosen_row1,choosen_C_range),'b');
hold on;
plot(crop_Temp_corr_D(choosen_row2,choosen_C_range),'g--');
hold on;
plot(Start-Start+1,crop_Temp_corr_D(choosen_row1,Start),'ko','Markersize',10);
plot(Middle-Start+1,crop_Temp_corr_D(choosen_row1,Middle),'ko','Markersize',10);
plot(End-Start+1,crop_Temp_corr_D(choosen_row1,End),'ko','Markersize',10);
plot(Start-Start+1,crop_Temp_corr_D(choosen_row2,Start),'ko','Markersize',10);
plot(Middle-Start+1,crop_Temp_corr_D(choosen_row2,Middle),'ko','Markersize',10);
plot(End-Start+1,crop_Temp_corr_D(choosen_row2,End),'ko','Markersize',10);
% plot(crop_Temp_corr_D(choosen_row3,choosen_C_range),'k');
grid on;
ylim(temp_limit);
xlim([0 120]);
ylabel('Temperature ( ^oC)');
title('3D thermography corrected by dielectric model');

figure(26);
%subplot(2,2,3);
plot(crop_Temp_corr_ND(choosen_row1,choosen_C_range),'b');
hold on;
plot(crop_Temp_corr_ND(choosen_row2,choosen_C_range),'g--');
hold on;
plot(Start-Start+1,crop_Temp_corr_ND(choosen_row1,Start),'ko','Markersize',10);
plot(Middle-Start+1,crop_Temp_corr_ND(choosen_row1,Middle),'ko','Markersize',10);
plot(End-Start+1,crop_Temp_corr_ND(choosen_row1,End),'ko','Markersize',10);
plot(Start-Start+1,crop_Temp_corr_ND(choosen_row2,Start),'ko','Markersize',10);
plot(Middle-Start+1,crop_Temp_corr_ND(choosen_row2,Middle),'ko','Markersize',10);
plot(End-Start+1,crop_Temp_corr_ND(choosen_row2,End),'ko','Markersize',10);
%plot(crop_Temp_corr_ND(choosen_row3,choosen_C_range),'k');
grid on;
ylim(temp_limit);
xlim([0 120]);
ylabel('Temperature ( ^oC)');
title('3D thermography corrected by non-dielectric model');

figure(27)
%subplot(2,2,4);
plot(crop_Temp_corr_Cosine(choosen_row1,choosen_C_range),'b');
hold on;
plot(crop_Temp_corr_Cosine(choosen_row2,choosen_C_range),'g');
hold on;
% plot(crop_Temp_corr_Cosine(choosen_row3,choosen_C_range),'k');
grid on;
ylim([30 39]);
ylabel('Temperature ( ^oC)');
title('3D thermography corrected by cosine model');


%% 
% plot the 1D temperature proflies along blue line together (uncorr, D, ND)
figure(28);
plot(crop_Temp_uncorr(choosen_row1,choosen_C_range),'b--');
hold on;
plot(crop_Temp_corr_D(choosen_row1,choosen_C_range),'b-');
hold on;
plot(crop_Temp_corr_ND(choosen_row1,choosen_C_range),'b-*');
grid on;
plot(Start-Start+1,crop_Temp_uncorr(choosen_row1,Start),'ko','Markersize',10);
plot(Middle-Start+1,crop_Temp_uncorr(choosen_row1,Middle),'ko','Markersize',10);
plot(End-Start+1,crop_Temp_uncorr(choosen_row1,End),'ko','Markersize',10);
ylim(temp_limit);
xlim([0 120]);
ylabel('Temperature ( ^oC)');
xlabel('No. of Pixels');

legend('uncorrected','dielectric model','non-dielectric model' );


% plot the 1D temperature proflies along green line together (uncorr, D, ND)
figure(29);
plot(crop_Temp_uncorr(choosen_row2,choosen_C_range),'g--');
hold on;
plot(crop_Temp_corr_D(choosen_row2,choosen_C_range),'g-');
hold on;
plot(crop_Temp_corr_ND(choosen_row2,choosen_C_range),'g-*');
grid on;
plot(Start-Start+1,crop_Temp_uncorr(choosen_row2,Start),'ko','Markersize',20);
plot(Middle-Start+1,crop_Temp_uncorr(choosen_row2,Middle),'ko','Markersize',20);
plot(End-Start+1,crop_Temp_uncorr(choosen_row2,End),'ko','Markersize',20);
ylim(temp_limit);
xlim([0 120]);
xlabel('No. of Pixels');
ylabel('Temperature ( ^oC)');
legend('uncorrected','dielectric model','non-dielectric model' );

% plot the correction differences for the three cases (uncorr, D, ND)

figure(30);
plot(crop_Temp_corr_D(choosen_row1,choosen_C_range)-crop_Temp_uncorr(choosen_row1,choosen_C_range),'b-');
hold on;
plot(crop_Temp_corr_ND(choosen_row1,choosen_C_range)-crop_Temp_uncorr(choosen_row1,choosen_C_range),'b-*');
hold on;
plot(crop_Temp_corr_D(choosen_row2,choosen_C_range)-crop_Temp_uncorr(choosen_row2,choosen_C_range),'g-');
hold on;
plot(crop_Temp_corr_ND(choosen_row2,choosen_C_range)-crop_Temp_uncorr(choosen_row2,choosen_C_range),'g-*');
hold on;
legend('dielectric model','non-dielectric model','dielectric model','non-dielectric model' )
grid on;
xlabel('No. of Pixels');
xlim([0 120]);
ylabel('Temperature ( ^oC)');


%% plot the 3D profile of the back-projected IR image before/after correction
% pick the lower bound and upper bound of skin range in the Image frame
x_l=skin_range_R(1,1);
x_u=skin_range_R(1,length(skin_range_R));
y_l=skin_range_C(1,1);
y_u=skin_range_C(1,length(skin_range_C));


figure(31);
subplot(2,1,1);
[X,Y]=meshgrid(1:x_u-x_l+1,1:y_u-y_l+1);
contour3(X,Y,Temp_IR_uncorr(skin_range_R,skin_range_C)',10);
surface(X,Y,Temp_IR_uncorr(skin_range_R,skin_range_C)','FaceColor','interp','Edgecolor','none');
caxis(Temp_limit);
view(60,60);
colormap jet;
%set(gca,'Clim',cblim)
colorbar;


%% 2D mapping for before/after correction (with different model)
figure(32);
subplot(2,2,1);
scatter(uv_IR_m(1,:),uv_IR_m(2,:),20,Im_pro_IR,'filled');
axis equal;
title('uncorrected');
colormap('jet');
xlim([0 600]);
ylim([50 350]);
caxis(Temp_limit);
view(0,270);
colorbar;

%figure;
subplot(2,2,2);
scatter(uv_IR_m(1,:),uv_IR_m(2,:),20,scaling_factor_ND*Im_pro_IR_corr_smooth_ND,'filled');
axis equal;
title('correction by non-dielectric model');
colormap('jet');
xlim([0 600]);
ylim([50 350]);
caxis(Temp_limit);
view(0,270);
colorbar;

%figure;
subplot(2,2,3);
scatter(uv_IR_m(1,:),uv_IR_m(2,:),20,scaling_factor_D*Im_pro_IR_corr_smooth_D,'filled');
axis equal;
title('correction by dielectric model');
colormap('jet');
xlim([0 600]);
ylim([50 350]);
caxis(Temp_limit);
view(0,270);
colorbar;

%figure;
subplot(2,2,4);
scatter(uv_IR_m(1,:),uv_IR_m(2,:),20,scaling_factor_Cosine*Im_pro_IR_corr_smooth_Cosine,'filled');
axis equal;
title('IR pixel mapping in 2D space (Cosine model)');
colormap('jet');
xlim([0 600]);
ylim([50 350]);
caxis(Temp_limit);
view(0,270);
colorbar;
