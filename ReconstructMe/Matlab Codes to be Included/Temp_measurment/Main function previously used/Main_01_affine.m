% Start with clean Matlab Workspace
clear all; close all; clc

% Compile the fast Lucas Kanade c-code tracking
mex LucasKanadeInverseAffine.c -v

%% load White light & IR image file

% Hemagioma- S4
stamp='164735';  % patient S4
W_file_name='IMG_0573';
initial= 40;
final= 80; 
SSstamp='164524';  % steady state image for comparison
SS_N= 1;
Frame_int=1;
time_step=1;

lowT=30;
highT=37;


% Convert the pixel intensity to temperature scale
[ImR1]= Read_IR_raw_data(stamp,initial);
[ImT1]= TempConvert(ImR1);

% steady state image for comparison 
[ImRS]= Read_IR_raw_data(SSstamp,SS_N);
[ImTS]= TempConvert(ImRS);

% Get the steady-state frame
I_S= mat2gray(ImTS);

% Get the first movie frame
I= mat2gray(ImT1);

%% Compute the mean value of the temperature in the steady-state image
med_temp= median(median(ImTS));
Var_ROI= ImTS- med_temp;

% display the temperature variation map in the steady state
figure;
%imshow(Var_ROI, [-1 max(max(Var_ROI))]);
imagesc(ImTS,[lowT highT]);
colormap('jet');
hold on;
colorbar;


%% Select Corners Correspondence betweeen White Light and IR image Manually
%[CW,CIR]=Corners_manual_selection_Wh2IR(W_file_name,ImT1);
% save('CW.mat','CW');
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
% [Lesion_B_W,imW_resize]=lesion_outlines_BW(W_file_name);
% 
% save('Lesion_B_W.mat','Lesion_B_W');
% save('imW_resize.mat','imW_resize');
%% load pre-segmented lesion outline in whitelight image if necessary 
load ('CW.mat');
load ('CIR40.mat');
% % 
load('TempPos40_68.mat');
load('Lesion_B_W.mat');
load('imW_resize.mat');
% 
% load('Lesion_B_IR_mapped.mat');
% load('CIR_mapped.mat');

%% Alternatively, activate if using the "segmented lesion corrdinates" for observation, and only the healthy tissue should be selecte
[l_p,h_p,Lesion_B_IR_mapped,CIR_mapped,point_num]=lesion_points_mapping_manual(Lesion_B_W,ImT1,CW,CIR,lowT,highT);


% load('l_p.mat');
% load('h_p.mat');

%%
% Display the segmented lesion boundary in white-light image
% figure;
% imshow(imW_resize);
% hold on;
% plot(Lesion_B_W(1,:),Lesion_B_W(2,:),'r.');
% title('the segmented lesion boundary in the white-light image');

%%
% Make a struct to store templates
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

figure;
for i= initial:final
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
        %imagesc(TemplateData(1).image,[lowT highT]);
        imagesc(ImT_i,[lowT highT]);
        colormap('jet');
        hold on;
        colorbar;
           
    % Show the location of the templates in the movie frame
    for t=1:length(TemplateData)
          h(t)=plot(lesion_p{i-initial+1}(1:point_num,2),lesion_p{i-initial+1}(1:point_num,1),'go','MarkerSize',5,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          h1(t)=plot(healthy_p{i-initial+1}(1:point_num,2),healthy_p{i-initial+1}(1:point_num,1),'ko','MarkerSize',5,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          %drawnow;
    end
title(strcat('t=',num2str((i-initial)*time_step),'s'));

LT_sum=0;
HT_sum=0;

for p=1:point_num
    
LT(i-initial+1,p)= ImT_i(round(lesion_p{i-initial+1}(p,1)),round(lesion_p{i-initial+1}(p,2)));
HT(i-initial+1,p)= ImT_i(round(healthy_p{i-initial+1}(p,1)),round(healthy_p{i-initial+1}(p,2)));

LT_sum= LT(i-initial+1,p)+LT_sum; 
HT_sum= HT(i-initial+1,p)+HT_sum; 

% LT_sum= ImT_i(round(lesion_p{i-initial+1}(t,1)),round(lesion_p{i-initial+1}(t,2)))+LT_sum; 
% HT_sum= ImT_i(round(healthy_p{i-initial+1}(t,1)),round(healthy_p{i-initial+1}(t,2)))+HT_sum; 
% 
end
Temp_rec_l(i-initial+1,1)=LT_sum/size(l_p,1);
Temp_rec_h(i-initial+1,1)=HT_sum/size(l_p,1);

F(i-initial+1)=getframe(gcf);         
end

%% Show the "averaged" recovery curve for both lesion and healthy tissue region
figure,
plot(0:time_step:((final-initial)/Frame_int)*time_step,Temp_rec_l(:,1),'g*'); 
hold on;
grid on;
plot(0:time_step:((final-initial)/Frame_int)*time_step,Temp_rec_h(:,1),'k.'); 
legend('lesion ','surrounding tissue ',2);
title('Recovery Temperature (with motion tracking)');
xlabel('recovery time (seconds)');
ylim([33,36]);
ylabel('temperature (deg C)');

% Show the recovery curve for each lesion & healthy tissue point,
% respective, to observe the variation range
figure;
plot(0:time_step:((final-initial)/Frame_int)*time_step,LT(:,1:point_num),'g*'); 
hold on;
grid on;
plot(0:time_step:((final-initial)/Frame_int)*time_step,HT(:,1:point_num),'k.'); 
title('Recovery Temperature (9 points)');
xlabel('recovery time (seconds)');
ylabel('temperature (deg C)');


%% Show the thermal contrast in recovery phase
figure;
DeltaT= Temp_rec_l(:,1)-Temp_rec_h(:,1);
plot(0:2:(final-initial+1)*2-1,DeltaT,'b*','MarkerSize',5); 
hold on;
grid on;

%% replay the movie
ht=figure;
movie(ht,F,10,2);

%% Save the patient-dependent computed variables  
% 
% save('TempPos1.mat','TempPos1');
% save('Lesion_B_W.mat','Lesion_B_W');
% save('Lesion_B_IR_mapped.mat','Lesion_B_IR_mapped');
% save('CIR_mapped.mat','CIR_mapped');
% save('imW_resize.mat','imW_resize');

%% Save the lesion & healthy tissue coordinate and their recovery curve 
% save('l_p.mat','l_p');
% save('h_p.mat','h_p');
% save('lesion_temp_E_14C60s.mat','Temp_rec_l');
% save('surrounding_temp_E_14C60s.mat','Temp_rec_h');
% save('9pt_lesion_temp_E_14C60s.mat','LT');
% save('9pt_surrounding_temp_E_14C60s.mat','HT');

