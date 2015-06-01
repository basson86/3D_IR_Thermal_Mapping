% Start with clean Matlab Workspace
clear all; close all; clc

% Compile the fast Lucas Kanade c-code tracking
mex LucasKanadeInverseAffine.c -v

%% load IR steady-state image & IR image file
%020912-------------------------------------
% stamp='115201';  % practice trial AKB- arm 2012/02/09 -hand motion_roation_SS_3 test- higher frame rate 2 frames/sec.
% initial= 2;
% final= 180; 
% Update_Frame=70;
% SSstamp='115201';  % steady state image for comparison
% SS_N= 1;
% time_step=0.5;
% Frame_int= 3;

% stamp='121001';  % practice trial 2012/02/09 -hand motion_roation_3 gel cooling- higher frame rate 2 frames/sec.
% initial= 105;
% final= 350;%255; 
% Update_Frame=181;
% SSstamp='121001';  % steady state image
% SS_N= 1;
% time_step=0.5;
% Frame_int= 2;

stamp='120339';  % practice trial 2012/02/09 -Hand_tranlation_gelcooling_3- higher frame rate 2 frames/sec.
initial= 103;
final= 250; 
Update_Frame=181;
SSstamp='120339';  % steady state image
SS_N= 1;
time_step=0.5;
Frame_int= 2;


% 021612--------------------------------------------------------------
% stamp='132816';  % practice trial 2012/02/16 -Hand_2_markers_round_1
% initial= 101;
% final= 200; 
% Update_Frame=120;
% SSstamp='132816';  % steady state image
% SS_N= 1;
% time_step=1;
% Frame_int= 1;

% stamp='131830';  % practice trial 2012/02/16 -Hand_2_markers_square
% W_file_name='IMG_0449';
% initial= 70;
% final= 190; 
% Update_Frame=100;
% SSstamp='131830';  % steady state image 
% SS_N= 1;
% time_step=1;
% Frame_int= 1;
%%
% stamp='100057';  % Patient 04
% initial= 49;
% final= 139; 
% Update_Frame=100;
% SSstamp='100023';  % steady state image
% SS_N= 1;
% time_step=2;
% Frame_int= 1;

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
% [CIR_S,CIR]=Corners_manual_selection_IRs2IR(ImRS,ImT1);
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
TemplateData(1).p=[0 0 0 0 center(1) center(2)]; % center(1) = row, center (2) = column 

TemplateData(1).weight=mat2gray(I(TempPos1(1,1):TempPos1(1,2),TempPos1(2,1):TempPos1(2,2)));
TemplateData(1).image=I(TempPos1(1,1):TempPos1(1,2),TempPos1(2,1):TempPos1(2,2));

c_to_lesion=[l_p(:,2)-center(1), l_p(:,1)-center(2)];
c_to_healthy=[h_p(:,2)-center(1), h_p(:,1)-center(2)];

% the relative position of Template Corners
c_to_tempPos= [TempPos1(1,:)-[center(1),center(1)]; TempPos1(2,:)-[center(2),center(2)]];


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


for i= initial:Frame_int:final
    % Get the a movie frame        
    [ImR_i]= Read_IR_raw_data(stamp,i);
    [ImT_i]= TempConvert(ImR_i);
     
    I_i= mat2gray(ImT_i);
    
    
    % Do the tracking for all templates, using Lucas Kanade Inverse Affine
    for t=1:length(TemplateData)
        [TemplateData(t).p,ROIimage,T_error(i,t)]=LucasKanadeInverseAffine(I_i,TemplateData(t).p,TemplateData(t).image,TemplateData(t).weight,Options);
    end
    n_p=TemplateData(1).p ;
    
    % Update the Affine Transformation Matrix
     M = [ 1+n_p(1),n_p(3),n_p(5); 
           n_p(2),   1+n_p(4),n_p(6); 
             0      0      1];
    % Warp the position of lesion and healthy tissue to the coordinate of input frame
    for n=1:point_num;
      lesion_coor= M*[c_to_lesion(n,:)';1];
      healthy_coor= M*[c_to_healthy(n,:)';1];
      lesion_p{(i-initial)/Frame_int+1}(n,:)= lesion_coor(1:2,1)';
      healthy_p{(i-initial)/Frame_int+1}(n,:)= healthy_coor(1:2,1)';
    end
    
    
    % Update the new template position at "Update_Frame" toward steady
    % state
    if (i==Update_Frame)
    TempPos1_coor1= M*[c_to_tempPos(1,:)';1];
    TempPos1_coor2= M*[c_to_tempPos(2,:)';1];
    
    TempPos1_coor= M*[c_to_tempPos;1,1];
    
    TempPos1(1,:)=TempPos1_coor(1,1:2);
    TempPos1(2,:)=TempPos1_coor(2,1:2);
    % Update the new template image based on new template corner location
    TemplateData(1).weight=mat2gray(I_i(TempPos1(1,1):TempPos1(1,2),TempPos1(2,1):TempPos1(2,2)));
    TemplateData(1).image=I_i(TempPos1(1,1):TempPos1(1,2),TempPos1(2,1):TempPos1(2,2));
    end 
      
    % Show the movie frame
        subplot(1,2,1);
        imagesc(TemplateData(1).image);
        subplot(1,2,2);
 %       imagesc(ImT_i,[10 35]);
         imshow(ImT_i,[25 35]);
        colormap('jet');
        hold on;
    
    % Show the location of the templates in the movie frame
    for t=1:length(TemplateData)
          plot(lesion_p{(i-initial)/Frame_int+1}(1:point_num,2),lesion_p{(i-initial)/Frame_int+1}(1:point_num,1),'k.','MarkerSize',5,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          plot(healthy_p{(i-initial)/Frame_int+1}(1:point_num,2),healthy_p{(i-initial)/Frame_int+1}(1:point_num,1),'mo','MarkerSize',2,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          drawnow;
    end
    title(strcat('frame No.=',num2str(i),'  Recovery time=',num2str((i-initial+1)*time_step),'sec'));

LT_sum=0;
HT_sum=0;

for p=1:point_num
    
LT((i-initial)/Frame_int+1,p)= ImT_i(round(lesion_p{(i-initial)/Frame_int+1}(p,1)),round(lesion_p{(i-initial)/Frame_int+1}(p,2)));
HT((i-initial)/Frame_int+1,p)= ImT_i(round(healthy_p{(i-initial)/Frame_int+1}(p,1)),round(healthy_p{(i-initial)/Frame_int+1}(p,2)));
LT_sum= LT((i-initial)/Frame_int+1,p)+LT_sum; 
HT_sum= HT((i-initial)/Frame_int+1,p)+HT_sum; 

end
Temp_rec_l((i-initial)/Frame_int+1,1)=LT_sum/size(l_p,1);
Temp_rec_h((i-initial)/Frame_int+1,1)=HT_sum/size(l_p,1);

F((i-initial)/Frame_int+1)=getframe(gcf);         
end

%% Show the "averaged" recovery curve for both lesion and healthy tissue region
figure,
plot(0:time_step*Frame_int:(final-initial)*time_step,Temp_rec_l(:,1),'b*'); 
hold on;
grid on;
plot(0:time_step*Frame_int:(final-initial)*time_step,Temp_rec_h(:,1),'m.'); 
legend('lesion ','surrounding tissue ',2);
title('Recovery Temperature (with motion correction)');
xlabel('recovery time (seconds)');
ylabel('temperature (deg C)');

% Show the recovery curve for each lesion & healthy tissue point,
% respectively, to observe the variation range
figure;
plot(0:time_step*Frame_int:(final-initial)*time_step,LT(:,1:point_num),'b*'); 
hold on;
grid on;
plot(0:time_step*Frame_int:(final-initial)*time_step,HT(:,1:point_num),'m.'); 
title('Recovery Temperature (9 points)');
xlabel('recovery time (seconds)');
ylabel('temperature (deg C)');

%% replay the movie
ht=figure;
movie(ht,F,30,2);





