%% Given the border points of the lesion, output the region pixels coordinate that fill up the lesion
function[Region_Points]=Border_to_Region_points(Border_Points)

%Border_Points=Lesion_B_IR_mapped(1:2,:);

%% sort the border points cloud by its Y (column) coordinate
Border_Points_sort_by_Y=sort(Border_Points,2);

% rounding the point coordinate to integers in image 
Border_P_floor = floor(Border_Points_sort_by_Y);

%initialize variables for region point collection
x_index=Border_P_floor(1,1); % the leftmost x coordinate of lesion points
y_index=  [];

%the x-range and corresponding y-range of lesion point
x_range = [];
y_range = [];
%%
for i= 1:length (Border_P_floor)
    
    % at the same x coordinate, find the range of y-coordinate
    if Border_P_floor(1,i) == x_index
        
        y_index = cat(1,y_index,Border_P_floor(2,i));
        
    else % if the x coordinate is different, store the range found from above
        x_range = cat(1,x_range,x_index);
        y_range = cat(1,y_range,[min(y_index),max(y_index)]);
        
        % increment th x-index to find the next y-range at that x
        % coordinate
        y_index = Border_P_floor(2,i);
        x_index = Border_P_floor(1,i);
    end 
    
    % at the last iteration, store the final x and y-range found
    if i==length (Border_P_floor)
    x_range = cat(1,x_range,x_index);
    y_range = cat(1,y_range,[min(y_index),max(y_index)]);
    end
    
end

Region_Points = [];
%% deploy the defined "lesion region points" by the pair of x-range & y-range parameters
for i= 1: length(x_range)
    
    for j= y_range(i,1):1:y_range(i,2)
    
    Region_Points = [Region_Points,[x_range(i);j]];  
    
    end
        
end
        
%% for de-bugging use
% figure;
% imagesc(ImT1,[lowT highT]),colormap('jet');
% colorbar;
% hold on;
% scatter(CIR_mapped(1,:), CIR_mapped(2,:), 40, 'g', 'filled');
% scatter(Lesion_B_IR_mapped(1,:), Lesion_B_IR_mapped(2,:), 2, 'k', 'filled');
% scatter(Region_Points(1,:), Region_Points(2,:), 2, 'r', 'filled');    
