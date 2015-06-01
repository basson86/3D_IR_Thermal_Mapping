% Start with clean Matlab Workspace
clear all; close all; clc

% Compile the fast Lucas Kanade c-code tracking
mex LucasKanadeInverseAffine.c -v

%% load White light & IR image file
stamp='100057';  % patient 04
W_file_name='IMG_0010';
initial= 49;
final= 139;

% stamp='122128'; % patient 14
% W_file_name='DSC00778';
% initial= 27;
% final= 117;

% stamp='084521'; % patient 30
% W_file_name='DSC00834';
% initial= 31;
% final= 121;

% stamp='090640'; % patient 20
% W_file_name='DSC00796';
% initial= 28;
% final= 78;



% Conver the pixel intensity to temperature scale
[ImR1]= Read_IR_raw_data(stamp,initial);
[ImT1]= TempConvert(ImR1);


%% If the pre-computed template position and the outline are available, load it for use is available, load it for use
% patient-dependent variable
load('CW.mat');
load('CIR.mat');
load('TempPos1.mat');
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


%% Display the segmented lesion boundary in white-light image, mapped lesion & corners in IR image
%  and dislpay the selected Template Image, and the selected lesion & healthy location in IR image

figure;
subplot(2,2,1);
imshow(imW_resize);
hold on;
scatter(CW(1,:), CW(2,:), 30, 'k','filled');
scatter(Lesion_B_W(1,:), Lesion_B_W(2,:), 8, 'r','filled');
colormap('gray');

subplot(2,2,2);
imcontour(ImT1,25); 
hold on;
scatter(round(CIR_mapped(1,:)), round(CIR_mapped(2,:)), 40, 'k','filled');
scatter(Lesion_B_IR_mapped(1,:), Lesion_B_IR_mapped(2,:), 8, 'r','filled');
scatter(l_p(:,1), l_p(:,2), 8, 'r');

subplot(2,2,3);
imshow(I(TempPos1(1,1):TempPos1(1,2),TempPos1(2,1):TempPos1(2,2)));

subplot(2,2,4);
imshow(I);
hold on;
scatter(l_p(:,1), l_p(:,2), 20, 'r','MarkerFaceColor',[1,1,1]);
scatter(h_p(:,1), h_p(:,2), 8, 'm','MarkerFaceColor',[1,1,1]);
colormap('jet');



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


%% -Loop through the movie frames


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

% plot the mean value of recovery curves from both lesion and healthy
% points group
LT_sum=0;
HT_sum=0;
for p=1:point_num

LT(i-initial+1,p)= ImT_i(round(lesion_p{i-initial+1}(p,1)),round(lesion_p{i-initial+1}(p,2)));
HT(i-initial+1,p)= ImT_i(round(healthy_p{i-initial+1}(p,1)),round(healthy_p{i-initial+1}(p,2)));

LT_sum= LT(i-initial+1,p)+LT_sum; 
HT_sum= HT(i-initial+1,p)+HT_sum;     
    
%LT_sum= ImT_i(round(lesion_p{i-initial+1}(t,1)),round(lesion_p{i-initial+1}(t,2)))+LT_sum; 
%HT_sum= ImT_i(round(healthy_p{i-initial+1}(t,1)),round(healthy_p{i-initial+1}(t,2)))+HT_sum; 
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
movie(ht,F,5,1);

%% save the resulting recovery curve
% save('lesion_temp_E_14C60s.mat','Temp_rec_l');
% save('surrounding_temp_E_14C60s.mat','Temp_rec_h');
% save('9pt_lesion_temp_E_14C60s.mat','LT');
% save('9pt_surrounding_temp_E_14C60s.mat','HT');
