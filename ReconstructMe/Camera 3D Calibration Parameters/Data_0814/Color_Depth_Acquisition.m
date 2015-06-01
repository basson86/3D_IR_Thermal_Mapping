%% 
clear all; close all;
info = imaqhwinfo('kinect');

%%
%get kinect information
info = imaqhwinfo('kinect');
info.DeviceInfo(1) 
info.DeviceInfo(2)

%% Create the objects for the color and depth sensors.Device 1 is the color sensor and Device 2 is the depth sensor.
%  %RGB camera:
vid = videoinput('kinect',1,'RGB_1280x960');
% 
%  %Depth camera:
vid2 = videoinput('kinect',2,'Depth_640x480');

%% set-up the parameters for Kinect camera


%%Get the source properties for the depth device.
srcDepth = getselectedsource(vid2);


%%Set the frames per trigger for both devices to 1.
vid.FramesPerTrigger = 1;
vid2.FramesPerTrigger = 1;


%%Set the trigger repeat for both devices to 200, inorder to acquire 201 frames from both the color sensor and the depthsensor.
vid.TriggerRepeat = 10000;
vid2.TriggerRepeat = 10000;


%%Configure the camera for manual triggering for bothsensors.
triggerconfig([vid vid2],'manual');


%%Start both video objects.
start([vid vid2]);


%%Trigger the devices, then get the acquired data.
% Trigger 200 times to get the frames.
% P.S 100 is around 10 sec


figure;
for i = 1:10   
    % Trigger both objects.
    trigger([vid vid2])
    % Get the acquired frames and metadata.
    [imgColor, ts_color, metaData_Color] = getdata(vid);
    [imgDepth, ts_depth, metaData_Depth] = getdata(vid2);
    

    % Normalize Depth map to 0..1 range
    D = imgDepth; D = double(D); 
    D= max(D(:))-D; 
    %D = D-min(D(:)); 
    D = D./(max(D(:)) + eps);
    
    % update the image pair acquired
    subplot(1,2,1), h1 = imshow(flipdim(imgColor,2));
    subplot(1,2,2), h2 = imshow(flipdim(D,2));
    
    pause(0.01);
    drawnow;    
    

end


%% Save the final updated image pair

% [C1_filename, pathname]=uiputfile({'*.jpg;*.tif;*.png;*.pgm','All Image Files';...
%           '*.*','All Files' },'Save Image')
%           
% imgColor_f=flipdim(imgColor,2);
% imwrite(imgColor_f,'0047.jpg','jpg');
% 
% imgDepth_f=flipdim(imgDepth,2);
% imwrite(imgDepth_f,'0047.pgm','pgm');


%%
delete(vid);
clear vid;
delete(vid2);
clear vid;
