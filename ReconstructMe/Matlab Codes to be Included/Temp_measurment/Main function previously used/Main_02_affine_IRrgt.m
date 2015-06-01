% Start with clean Matlab Workspace
clear all; close all; clc

% Compile the fast Lucas Kanade c-code tracking
mex LucasKanadeInverseAffine.c -v

%% load IR steady-state image & IR image file
% stamp='115445';  % practice trial TC- arm 2012/01/31 -trail 1 (steady-state in motion)
% W_file_name='hand_photo_marked';
% initial= 2;
% final= 30; 
% SSstamp='115445';  % steady state image for comparison
% SS_N= 1;

% stamp='120446';  % practice trial TC- arm 2012/01/31 -trail 2 
% W_file_name='hand_photo_marked';
% initial= 54;
% final= 114;
% SSstamp='120446';  % steady state image for comparison
% SS_N= 10;

% stamp='121540';  % practice trial TC- arm 2012/01/31 -trail 3
% W_file_name='hand_photo_marked';
% initial= 29;
% final= 83; 
% SSstamp='121253';  % steady state image for comparison
% SS_N= 2;

% stamp='122156';  % practice trial TC- arm 2012/01/31 -trail 4
% W_file_name='hand_photo_marked';
% initial= 46;
% final= 136; 
% SSstamp='122156';  % steady state image for comparison and registration
% SS_N= 2;

% stamp='123719';  % practice trial TC- arm 2012/01/31 -trail 5 (pin-point)
% W_file_name='hand_photo_marked';
% initial= 5;
% final= 30; 
% SSstamp='123719';  % steady state image for comparison
% SS_N= 4;

%---------------------------------------------
% stamp='125454';  % practice trial TC- arm 2012/02/02 -hand motion 1:translationg test
% W_file_name='IMG_0373';
% initial= 2;
% final= 45; 
% SSstamp='125454';  % steady state image for comparison
% SS_N= 1;

% stamp='125906';  % practice trial TC- arm 2012/02/02 -hand motion 2:roation test
% W_file_name='IMG_0373';
% initial= 3;
% final= 45; 
% SSstamp='125906';  % steady state image for comparison
% SS_N= 2;
% time_step=2;

%-------------------------------------
% stamp='115201';  % practice trial AKB- arm 2012/02/09 -hand motion_roation_SS_3 test- higher frame rate 2 frames/sec.
% W_file_name='IMG_0373';
% initial= 2;
% final= 180; 
% SSstamp='115201';  % steady state image for comparison
% SS_N= 1;
% time_step=0.5;

% stamp='114000';  % practice trial AKB- arm 2012/02/09 -hand motion_translation_SS_3 test- higher frame rate 2 frames/sec.
% W_file_name='IMG_0373';
% initial= 2;
% final= 180; 
% SSstamp='114000';  % steady state image for comparison
% SS_N= 1;
% time_step=0.5;

stamp='121001';  % practice trial 2012/02/09 -hand motion_roation_3 gel cooling- higher frame rate 2 frames/sec.
%W_file_name='IMG_0373';
initial= 105;
final= 225; 
SSstamp='121001';  % steady state image
SS_N= 1;
time_step=0.5;

%%
% Convert the pixel intensity to temperature scale
[ImR1]= Read_IR_raw_data(stamp,initial);
[ImT1]= TempConvert(ImR1);

% Get the first movie frame
I= mat2gray(ImT1);

% steady state image for comparison 
[ImRS]= Read_IR_raw_data(SSstamp,SS_N);
[ImTS]= TempConvert(ImRS);

% Get the steady-state frame
I_S= mat2gray(ImTS);


%% Registration from "Steady-State IR" to "IR frame series in motion"
% [CIR_S,CIR]=Corners_manual_selection_IRs2IR(ImRS,ImR1);
% save('CIR_S.mat','CIR_S');
% save('CIR.mat', 'CIR');

%% Selecting Template Image by 2 corners manually

%Show the movie frame for template selection
% figure, imshow(I,[]);
% colormap('jet');
% title('please select the boundary points of template');
% Temp= ginput(2);
% TempPos1=[Temp(:,2),Temp(:,1)]';
% save('TempPos1.mat','TempPos1');

%% Segment lesion outline in white light image using "Random Walker" algorithm
% [Lesion_IRs,imRs_gray]=lesion_outlines_IRs(ImRS);
% save('Lesion_IRs.mat','Lesion_IRs');
% save('imRs_gray.mat','imRs_gray');


%% load pre-segmented lesion outline in whitelight image if necessary 
load ('CIR.mat');
load ('CIR_S.mat');
load('TempPos1.mat');
load('Lesion_IRs.mat');
load('imRs_gray.mat');

%% Select the lesion point and healthy in the first IR image by mapping the lesion outline in IR steady-state image to IR images in motion 
% using homographic transformation; Alternatively, activate if using the "segmented lesion corrdinates" for observation, and only the healthy tissue should be selecte
[l_p,h_p,Lesion_B_IR_mapped,CIR_mapped,point_num]=lesion_points_mapping_auto(Lesion_IRs,ImT1,CIR_S,CIR);

% load('l_p.mat');
% load('h_p.mat');

%% Compute the mean value of the temperature in the steady-state image
mean_temp= median(median(ImTS));
Var_ROI= ImTS- mean_temp;

% display the temperature variation map in the steady state
figure;
imagesc(Var_ROI, [-1 max(max(Var_ROI))]);
hold on;
colorbar;
plot(Lesion_IRs(1,:), Lesion_IRs(2,:),'r.','MarkerSize',10);


%% Make a struct to store templates
TemplateData=struct;

point_num_l= size(l_p,1);
point_num_h= size(h_p,1);

% Pad the select templates with extra boundary pixels. These boundary
% pixels are not used for the actual template tracking. But to get
% more reliable image derivatives.
b=1;
padding=[-b,b;-b,b];
TempPos1=TempPos1+padding;


% -Set initial parameters of the templates
% -Set padded template image.
% 
% Backwards affine Transformation Matrix is used in Lucas Kanade Tracking
% with 6 parameters
% M    = [ 1+p(1) p(3)   p(5); 
%          p(2)   1+p(4) p(6); 
%          0      0      1];
%

center=[TempPos1(1,1)+TempPos1(1,2)-1 TempPos1(2,1)+TempPos1(2,2)-1]/2;
TemplateData(1).p=[0 0 0 0 center(1) center(2)];

TemplateData(1).weight=mat2gray(I(TempPos1(1,1):TempPos1(1,2),TempPos1(2,1):TempPos1(2,2)));
TemplateData(1).image=I(TempPos1(1,1):TempPos1(1,2),TempPos1(2,1):TempPos1(2,2));

c_to_lesion=[l_p(:,2)-center(1), l_p(:,1)-center(2)];
c_to_healthy=[h_p(:,2)-center(1), h_p(:,1)-center(2)];

% This weight function is used in the LK-Hessian and multiplied with the 
% error between in image and template. And is used to exclude unreliable pixels form
% the template tracking.

% LK Tracking Options (other options default)
Options.TranslationIterations=30;
Options.AffineIterations=0;
Options.RoughSigma=3;
Options.FineSigma=1.5;

% Make a colormap
cmap=hot(256);

% Matrix to store squared pixel error between template and ROI in
% movieframe after template tracking.
T_error=zeros(final-initial+1, length(TemplateData));


% -Loop through the movie frames

% load('healthy_p.mat');

for i= initial:1:final
    % Get the a movie frame        
    [ImR_i]= Read_IR_raw_data(stamp,i);
    [ImT_i]= TempConvert(ImR_i);
     
    I_i= mat2gray(ImT_i);

    % Do the tracking for all templates, using Lucas Kanade Inverse Affine
    for t=1:length(TemplateData)
        [TemplateData(t).p,ROIimage,T_error(i,t)]=LucasKanadeInverseAffine(I_i,TemplateData(t).p,TemplateData(t).image,TemplateData(t).weight,Options);
    end
    n_p=TemplateData(1).p ;
    
    % update the Affine Transformation Matrix
     M = [ 1+n_p(1),n_p(3),n_p(5); 
           n_p(2),   1+n_p(4),n_p(6); 
             0      0      1];
    % warp the position of lesion and healthy tissue to the coordinate of input frame
    for n=1:point_num;
      lesion_coor= M*[c_to_lesion(n,:)';1];
      healthy_coor= M*[c_to_healthy(n,:)';1];
      lesion_p{i-initial+1}(n,:)= lesion_coor(1:2,1)';
      healthy_p{i-initial+1}(n,:)= healthy_coor(1:2,1)';
    end
    
  
    % Show the movie frame
    if(i==initial), figure, %handle_imshow=imagesc(ImT_i,[25 35]); 
        handle_imshow=imshow(ImT_i,[25 35]); colormap('jet');
        hold on;
    else
        for t=1:length(TemplateData), delete(h(t)),delete(h1(t)); end; 
        %imagesc(ImT_i,[25 35]);
        imshow(ImT_i,[25 35]);
        colormap('jet');
    end
    colorbar;
    
    % Show the location of the templates in the movie frame
    for t=1:length(TemplateData)
          h(t)=plot(lesion_p{i-initial+1}(1:point_num,2),lesion_p{i-initial+1}(1:point_num,1),'r.','MarkerSize',8,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          h1(t)=plot(healthy_p{i-initial+1}(1:point_num,2),healthy_p{i-initial+1}(1:point_num,1),'bx','MarkerSize',8,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          drawnow;
    end
title(strcat('t=',num2str((i-initial+1)*time_step),'s'));

LT_sum=0;
HT_sum=0;

for p=1:point_num
    
LT(i-initial+1,p)= ImT_i(round(lesion_p{i-initial+1}(p,1)),round(lesion_p{i-initial+1}(p,2)));
HT(i-initial+1,p)= ImT_i(round(healthy_p{i-initial+1}(p,1)),round(healthy_p{i-initial+1}(p,2)));

LT_sum= LT(i-initial+1,p)+LT_sum; 
HT_sum= HT(i-initial+1,p)+HT_sum; 

end
Temp_rec_l(i-initial+1,1)=LT_sum/size(l_p,1);
Temp_rec_h(i-initial+1,1)=HT_sum/size(l_p,1);

F(i-initial+1)=getframe(gcf);         
end

%% Show the "averaged" recovery curve for both lesion and healthy tissue region
figure,
plot(0:time_step:(final-initial)*time_step,Temp_rec_l(:,1),'b*'); 
hold on;
grid on;
plot(0:time_step:(final-initial)*time_step,Temp_rec_h(:,1),'m.'); 
legend('lesion ','surrounding tissue ',2);
title('Recovery Temperature (with motion correction)');
xlabel('recovery time (seconds)');
ylabel('temperature (deg C)');

% Show the recovery curve for each lesion & healthy tissue point,
% respective, to observe the variation range
figure;
plot(0:time_step:(final-initial)*time_step,LT(:,1:point_num),'b*'); 
hold on;
grid on;
plot(0:time_step:(final-initial)*time_step,HT(:,1:point_num),'m.'); 
title('Recovery Temperature (9 points)');
xlabel('recovery time (seconds)');
ylabel('temperature (deg C)');

%% replay the movie
ht=figure;
movie(ht,F,30,2);





