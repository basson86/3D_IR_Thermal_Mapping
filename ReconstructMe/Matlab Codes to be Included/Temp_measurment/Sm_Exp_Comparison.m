clear all;
close all;

stamp=100057;
stnum=49;
endnum=140;
step=1;

recovery_time=180;
%% load temperature from simulation
fname1='14c';
frame2='60s';
frame3=' Const T cooling ';
Fthick='C2f2_w0.045';
lesion_temp_S=load (strcat(fname1,'_Rec_Const',frame2,'_skin_surface_lesion_',Fthick,'.txt'));
surrounding_temp_S=load (strcat(fname1,'_Rec_Const',frame2,'_skin_surface_surrounding_',Fthick,'.txt'));

% temperature difference in recovery phase (simulation)
SDelT(:,2)= lesion_temp_S(:,2)-surrounding_temp_S(:,2);
SDelT(:,1)= lesion_temp_S(:,1);

%% load temperature curve from experimental data ( after motion-correction processing)
load('l_p.mat');
load('h_p.mat');

load('lesion_temp_E_14C60s.mat');
load('surrounding_temp_E_14C60s.mat');
lesion_temp_E=Temp_rec_l;
surrounding_temp_E=Temp_rec_h;


ExDelT_mean(:,1)= 0:step*2:(endnum-stnum)*2-1;
ExDelT_mean(:,2)=lesion_temp_E-surrounding_temp_E;


%%
figure;
plot(lesion_temp_S(:,1),lesion_temp_S(:,2),'r');
hold on;
plot(surrounding_temp_S(:,1), surrounding_temp_S(:,2),'b');
hold on;
plot(0:step*2:(endnum-stnum)*2-1,lesion_temp_E,'r.');
hold on;
plot(0:step*2:(endnum-stnum)*2-1,surrounding_temp_E,'b.');
grid on;
xlabel('time (sec)');
ylabel('recovery temperature( ^oC)');
xlim([0,recovery_time]);
ylim([10,35]);
h = legend('lesion (Simulation)','healthy tissue','lesion (Experiment)','healthy tissue',2);

figure;
plot(SDelT(:,1),SDelT(:,2),'b-');
hold on;
plot(ExDelT_mean(:,1),ExDelT_mean(:,2),'b.','MarkerSize',14);
hold on;
grid on;
xlabel('time (sec)');
ylabel('temperature difference(degC)');
xlim([0,recovery_time]);
ylim([0,3]);
h = legend('Simulation','Mean of Experimental Data',2);
ht=title(strcat('Temperature difference between lesion & surrounding skin during recovery','(',fname1,' ',frame2,' ',Fthick,',',frame3,')'));

%%
point_num=size(l_p,1);
[ImR_r]= Read_IR_raw_data(num2str(stamp),num2str(stnum));
[ImT_r]= TempConvert(ImR_r);

figure;
imagesc(ImT_r); clims = [12 35]; 
hold on; colorbar;
plot(l_p(1:point_num,1),l_p(1:point_num,2),'ro','MarkerSize',5,'MarkerFaceColor',[1,1,1]);
plot(h_p(1:point_num,1),h_p(1:point_num,2),'mo','MarkerSize',5,'MarkerFaceColor',[1,1,1]);
