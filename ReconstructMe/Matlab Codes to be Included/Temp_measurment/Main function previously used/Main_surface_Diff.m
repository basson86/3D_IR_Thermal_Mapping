% Start with clean Matlab Workspace
clear all; close all; clc

% Compile the fast Lucas Kanade c-code tracking
mex LucasKanadeInverseAffine.c -v

%% load White light & IR image file
% stamp='100057';  % patient 04
% W_file_name='IMG_0010';
% initial= 49;
% final= 99;
% r_ext= 0;
% c_ext= 0;
% TempD_scale= 20;
% SSstamp='100023';  % steady state image for comparison
% SS_N= 1;



% stamp='084737';  % patient 12
% W_file_name='DSC00769';
% initial= 47;
% final= 137;
% r_ext= 0;
% c_ext= 0;
% TempD_scale= 15;


% stamp='122128'; % patient 14
% W_file_name='DSC00778';
% initial= 27;
% final= 117;
% r_ext= -4;
% c_ext= 0;
% TempD_scale= 15;

% stamp='084521'; % patient 30
% W_file_name='DSC00834';
% initial= 31;
% final= 121;
% r_ext= 0;
% c_ext= 0;
% TempD_scale= 10;
% 

% stamp='090640'; % patient 20
% W_file_name='DSC00796';
% initial= 28;
% final= 118;
% r_ext= 0;
% c_ext= 0;
% TempD_scale= 10;
% SSstamp='090557';  % steady state image for comparison
% SS_N= 1;


stamp='101112'; % patient 37
W_file_name='DSC00848';
initial= 52;
final= 140;
r_ext= 80;
c_ext= -3;
TempD_scale= 15;
SSstamp='101049';  % steady state image for comparison
SS_N= 1;


% Conver the pixel intensity to temperature scale
[ImR1]= Read_IR_raw_data(stamp,initial);
[ImT1]= TempConvert(ImR1);

% Read steady-state image for reference
[SS_R]= Read_IR_raw_data(SSstamp,SS_N);
[SS_T]= TempConvert(SS_R);

%% If the pre-computed template position and the outline are available, load it for use is available, load it for use
% patient-dependent variable
load('CW.mat');
load('CIR.mat');
load('TempPos1.mat');
load('TempPos_S.mat');
load('Lesion_B_W.mat');
load('imW_resize.mat');
load('Lesion_B_IR_mapped.mat');
load('CIR_mapped.mat');

% lesion & healthy tissue location ( case-dependent)
load('l_p.mat');
load('h_p.mat');

% Get the number of selected points
point_num= size(l_p,1);

% Get the first movie frame
I= mat2gray(ImT1);


%% Show the white light image for reference
figure;
imshow(imW_resize);
hold on;
scatter(CW(1,:), CW(2,:), 30, 'k','filled');
scatter(Lesion_B_W(1,:), Lesion_B_W(2,:), 8, 'r','filled');
colormap('gray');

% Display the steady-state image for comparison 
figure;
imagesc(SS_T(TempPos_S(1,1):TempPos_S(1,2),TempPos_S(2,1):TempPos_S(2,2)));
colormap('jet');
colorbar;


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

c_to_LB=[Lesion_B_IR_mapped(2,:)'-center(1), Lesion_B_IR_mapped(1,:)'-center(2)];

%Center_1= TemplateData(1).p(5:6);
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
%T_error=zeros(final-initial+1, length(TemplateData));


%% -Loop through the movie frames

% Take the first Temperature frame as the reference image for subtraction
[ImR_1]= Read_IR_raw_data(stamp,initial);
[ImT_1]= TempConvert(ImR_1);

ROI_temp_1=ImT_1(TempPos1(1,1):TempPos1(1,2)+r_ext,TempPos1(2,1):TempPos1(2,2)+c_ext);     


for i= initial:final
    % Get the a movie frame        
    [ImR_i]= Read_IR_raw_data(stamp,i);
    [ImT_i]= TempConvert(ImR_i);
     
    % Take the temperature mapping for the ROI
    I_i= mat2gray(ImT_i);
    
    % Do the tracking for all templates, using Lucas Kanade Inverse Affine
    for t=1:length(TemplateData)
        [TemplateData(t).p,ROIimage,T_error(i,t)]=LucasKanadeInverseAffine(I_i,TemplateData(t).p,TemplateData(t).image,TemplateData(t).weight,Options);
    end
       
    % Update the new template position after tracking
    TempPos1(1,:)= c_to_tempPos(1,:) + TemplateData(1).p(5);
    TempPos1(2,:)= c_to_tempPos(2,:) + TemplateData(1).p(6);
    
    
    C1= TempPos1(:,1)';
    C2= TempPos1(:,2)';
    
    %Locate the lesion and healthy tissue point by relative vector in
    %between
    for n=1:point_num;
    lesion_p{i-initial+1}(n,:)=TemplateData(1).p(5:6) +  c_to_lesion(n,:);
    healthy_p{i-initial+1}(n,:)=TemplateData(1).p(5:6) +  c_to_healthy(n,:);
    lesion_p_d{i-initial+1}(n,:)=TemplateData(1).p(5:6) +  c_to_lesion(n,:)-C1;
    healthy_p_d{i-initial+1}(n,:)=TemplateData(1).p(5:6) +  c_to_healthy(n,:)-C1;
    
    end
    
    for n= 1: size(Lesion_B_IR_mapped,2);
      LB_s(n,:)=TemplateData(1).p(5:6) + c_to_LB(n,:)-C1;    
    end
    
    % Image Subtraction using the first frame to observe the temperature elevation    
    ROI_temp_i=ImT_i(TempPos1(1,1):TempPos1(1,2)+r_ext,TempPos1(2,1):TempPos1(2,2)+c_ext);     
    ROI_temp_i_d=ROI_temp_i-ROI_temp_1;
    
    % Compute the mean value of the temperature over the domain of interest
    mean_temp= median(median(ROI_temp_i));
    Var_ROI= ROI_temp_i- mean_temp;
    
        
% Show the movie frame
    if(i==initial), figure, 
       handle_imshow=imagesc(ROI_temp_i_d, [0 TempD_scale]);
       %handle_imshow=mesh(ROI_temp_i_d, [0 TempD_scale]);
       %view(10,60);
       %handle_imshow=imagesc(Var_ROI, [-10 TempD_scale]);
       colormap('jet');  hold on;
    else
        for t=1:length(Var_ROI); end; 
        clims = [0 TempD_scale];
        imagesc(ROI_temp_i_d); 
        %mesh(ROI_temp_i_d); 
        %view(10,60)
        %imagesc(Var_ROI);
    end
    colorbar;
    plot(LB_s(:,2), LB_s(:,1),'r.','MarkerSize',3);

        

 % Show the location of the templates in the movie frame
    for t=1:length(TemplateData)
          h(t)=plot(lesion_p_d{i-initial+1}(1:point_num,2),lesion_p_d{i-initial+1}(1:point_num,1),'b.','MarkerSize',8,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          h1(t)=plot(healthy_p_d{i-initial+1}(1:point_num,2),healthy_p_d{i-initial+1}(1:point_num,1),'m.','MarkerSize',8,'MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
          drawnow;
    end
    
    
title(strcat('t=',num2str((i-initial+1)*2),'s'));

% plot the mean value of recovery curves from both lesion and healthy
% points group
LT_sum=0;
HT_sum=0;
for p=1:point_num

LT(i-initial+1,p)= ROI_temp_i_d(round(lesion_p_d{i-initial+1}(p,1)),round(lesion_p_d{i-initial+1}(p,2)));
HT(i-initial+1,p)= ROI_temp_i_d(round(healthy_p_d{i-initial+1}(p,1)),round(healthy_p_d{i-initial+1}(p,2)));

LT_sum= LT(i-initial+1,p)+LT_sum; 
HT_sum= HT(i-initial+1,p)+HT_sum;     
     
end
Temp_rec_l(i-initial+1,1)=LT_sum/size(l_p,1);
Temp_rec_h(i-initial+1,1)=HT_sum/size(l_p,1);

F(i-initial+1)=getframe(gcf); 


end

%% plot the "mean" recovery curves
figure,
plot(0:2:(final-initial+1)*2-1,Temp_rec_l(:,1),'b*'); 
hold on;
grid on;
plot(0:2:(final-initial+1)*2-1,Temp_rec_h(:,1),'m.'); 
legend('lesion ','surrounding tissue ',2);
title('Recovery Temperature (mean value)');
xlabel('recovery time (seconds)');
ylabel('temperature elevation(deg C)');

% Show the recovery curve for each lesion & healthy tissue point,
% respective, to observe the variation range
figure;
plot(0:2:(final-initial+1)*2-1,LT(:,1:point_num),'b*'); 
hold on;
grid on;
plot(0:2:(final-initial+1)*2-1,HT(:,1:point_num),'m.'); 
title('Recovery Temperature (9 points)');
xlabel('recovery time (seconds)');
ylabel('temperature elevation (deg C)');

%% Dimensionless Temperature Recovery Curve
l_max=max(Temp_rec_l);
h_max=max(Temp_rec_h);

figure;
plot(0:2:(final-initial+1)*2-1,Temp_rec_l(:,1)./h_max,'b*'); 
hold on;
grid on;
plot(0:2:(final-initial+1)*2-1,Temp_rec_h(:,1)./h_max,'m.');
title('Dimensionless Recovery Temperature (mean value)');
xlabel('recovery time (seconds)');
ylabel('dimensionless temperature (deg C)');



%% replay the movie
ht=figure;
movie(ht,F,5,2);


