% Start with clean Matlab Workspace
clear all; close all; clc

% Compile the fast Lucas Kanade c-code tracking
mex LucasKanadeInverseAffine.c -v

%% load White light & IR image file
stamp='100057';  % patient 04
W_file_name='IMG_0010';
initial= 49;
final= 139; 
SSstamp='100023';  % steady state image for comparison
SS_N= 1;


lowT=10;
highT=37;

% Convert the pixel intensity to temperature scale
[ImR1]= Read_IR_raw_data(stamp,initial);
[ImT1]= TempConvert(ImR1);

% Get the first movie frame
I= mat2gray(ImT1);

%% Select Corners Correspondence betweeen White Light and IR image Manually
[CW,CIR]=Corners_manual_selection_Wh2IR(W_file_name,ImR1);
% save('CW.mat','CW');
% save('CIR.mat', 'CIR');

figure;
imshow(ImT1,[lowT highT]),colormap('hot');
hold on;
scatter(CIR(1,:),CIR(2,:), 80, 'g', 'filled');




%% Selecting Template Image by 2 corners manually

%Show the movie frame for template selection
figure, imshow(I,[]);
colormap('jet');
title('please select the boundary points of template');
Temp= ginput(2);
TempPos1=[Temp(:,2),Temp(:,1)]';

%% Segment lesion outline in white light image using "Random Walker" algorithm
[Lesion_B_W,imW_resize]=lesion_outlines_BW(W_file_name);
hold on;
scatter(CW(1,:),CW(2,:), 40, 'r', 'filled');

%% load pre-segmented lesion outline in whitelight image if necessary 
load ('CW.mat');
load ('CIR.mat');

load('TempPos1.mat');
load('Lesion_B_W.mat');
load('imW_resize.mat');
% 
load('Lesion_B_IR_mapped.mat');
load('CIR_mapped.mat');


%% Select the lesion point and healthy in the first IR image by mapping the lesion outline in white light image to IR image 
% using homographic transformation, 

%"lesion-healthy pair" number to observe thermal contrast
% point_num=20;
%[l_p,h_p,Lesion_B_IR_mapped,CIR_mapped]=lesion_points_mapping(Lesion_B_W,ImT1,point_num);

%% Alternatively, activate if using the "segmented lesion corrdinates" for observation, and only the healthy tissue should be selecte
[l_p,h_p,Lesion_B_IR_mapped,CIR_mapped,point_num]=lesion_points_mapping_auto(Lesion_B_W,ImT1,CW,CIR,lowT,highT);
% load('l_p.mat');
% load('h_p.mat');

%%
% Display the segmented lesion boundary in white-light image
% figure;
% imshow(imW_resize);
% hold on;
% plot(Lesion_B_W(1,:),Lesion_B_W(2,:),'r.');
% title('the segmented lesion boundary in the white-light image');


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


for i= initial:final
    % Get the a movie frame        
    [ImR_i]= Read_IR_raw_data(stamp,i);
    [ImT_i]= TempConvert(ImR_i);
     
    I_i= mat2gray(ImT_i);

    % Do the tracking for all templates, using Lucas Kanade Inverse Affine
    for t=1:length(TemplateData)
        [TemplateData(t).p,ROIimage,T_error(i,t)]=LucasKanadeInverseAffine(I_i,TemplateData(t).p,TemplateData(t).image,TemplateData(t).weight,Options);
    end
    
    %
    for n=1:point_num;
    lesion_p{i-initial+1}(n,:)=TemplateData(1).p(5:6) +  c_to_lesion(n,:);
    healthy_p{i-initial+1}(n,:)=TemplateData(1).p(5:6) +  c_to_healthy(n,:);
    end
    %
    
    % Show the movie frame
    if(i==initial), figure, handle_imshow=imagesc(ImT_i); clims = [12 35]; hold on
    else
        for t=1:length(TemplateData), delete(h(t)),delete(h1(t)); end; 
        clims = [12 35];
        imagesc(ImT_i);
        %set(handle_imshow,'Cdata',ImT_i);        
    end
    colorbar;
    
    % Show the location of the templates in the movie frame
    for t=1:length(TemplateData)
          h(t)=plot(lesion_p{i-initial+1}(1:point_num,2),lesion_p{i-initial+1}(1:point_num,1),'bo','MarkerSize',5,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          h1(t)=plot(healthy_p{i-initial+1}(1:point_num,2),healthy_p{i-initial+1}(1:point_num,1),'mo','MarkerSize',5,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          drawnow;
    end
title(strcat('t=',num2str((i-initial+1)*2),'s'));

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
plot(0:2:(final-initial+1)*2-1,Temp_rec_l(:,1),'b*'); 
hold on;
grid on;
plot(0:2:(final-initial+1)*2-1,Temp_rec_h(:,1),'m.'); 
legend('lesion ','surrounding tissue ',2);
title('Recovery Temperature (with motion correction)');
xlabel('recovery time (seconds)');
ylabel('temperature (deg C)');

% Show the recovery curve for each lesion & healthy tissue point,
% respective, to observe the variation range
figure;
plot(0:2:(final-initial+1)*2-1,LT(:,1:point_num),'b*'); 
hold on;
grid on;
plot(0:2:(final-initial+1)*2-1,HT(:,1:point_num),'m.'); 
title('Recovery Temperature (9 points)');
xlabel('recovery time (seconds)');
ylabel('temperature (deg C)');

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

