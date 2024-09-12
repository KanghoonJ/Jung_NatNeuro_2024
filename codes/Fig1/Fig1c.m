%% Code for Fig1c
% Created by Kanghoon Jung
% Please contact kanghoon.jung AT alleninstitute.org with any bug reports, questions or feature requests for this program. 

clear all
close all

data_folder = 'C:\Users\jungk\Google Drive\Github\Jung_NatNeuro_2024\codes\Fig1\';
cd(data_folder)

%% __init__
db = [];
load("fig1c_data.mat")

Sorted_Event_Sequence_Table = Gen_Event_Sequence_Table(db);

Threat_on_time = Sorted_Event_Sequence_Table.Event_time(find(ismember(Sorted_Event_Sequence_Table.Event_type,'Threat_on')));
Threat_on_Frames = find(ismember(Sorted_Event_Sequence_Table.Event_type,'Threat_on'));

ROI_on_time = []; %% Later one among Movement_on and Threat_on
ROI_off_time = []; %% Sooner one among Shelter_in and Threat_off

ROI_Frames = find(ismember(Sorted_Event_Sequence_Table.Event_type,'Threat_on'));
for(i=1:size(Threat_on_Frames,1))    
    j = Threat_on_Frames(i,1);    
    ROI_on_time(i,1) = Sorted_Event_Sequence_Table.Event_time(j);
    while j<=size(Sorted_Event_Sequence_Table,1)
        if(ismember(Sorted_Event_Sequence_Table.Event_type(j), {'Movement_on','Flight_on'}))
            if(Sorted_Event_Sequence_Table.Event_time(j)-ROI_on_time(i,1)<10)
                ROI_on_time(i,1) = Sorted_Event_Sequence_Table.Event_time(j);     
            end
            break;
        end
        j = j+1;    
    end

    k = Threat_on_Frames(i,1);
    while k<=size(Sorted_Event_Sequence_Table,1)
        if(ismember(Sorted_Event_Sequence_Table.Event_type(k), 'Threat_off'))
            ROI_off_time(i,1) = Sorted_Event_Sequence_Table.Event_time(k);     
            break;
        elseif(ismember(Sorted_Event_Sequence_Table.Event_type(k), 'Shelter_in'))
            ROI_off_time(i,1) = Sorted_Event_Sequence_Table.Event_time(k);    
            break;
        end
        k = k+1;    
    end
end

ROI_time_range_set = [ROI_on_time, ROI_off_time];
Conc_ROI_Frames = [];
for(k=1:size(Threat_on_Frames,1))    
    ROI_time_range = [ROI_time_range_set(k,1) ROI_time_range_set(k,2)];
    ROI_Behavior_data_table = Extract_ROI_Track(db, ROI_time_range);
    %% Plots
    % 1) Traces on Arena with direction
    Plot_Traces(db, ROI_time_range)

    % 2) Normalized Traces on Arena with direction
    Plot_Norm_dist2shelter(db, ROI_time_range)

    ROI_Frames = Find_Frames_in_TS(ROI_time_range, db.Behavior_data.Behavior_data_table.Timestamp);
    Conc_ROI_Frames = [Conc_ROI_Frames;ROI_Frames];
end

    % 3-4) Plot WindRose Goal Angle
    % Flight under threat    
    Plot_WindRose_Goal_angle(db, intersect(Conc_ROI_Frames, find(db.Behavior_data.Behavior_data_table.Flight_Adjusted)))
    set(gcf,'position',[200 200 300 300])

           
    % Normal exploration
    Movement_Frames = find(db.Behavior_data.Behavior_data_table.Movement_Adjusted);
    Shelter_out_Frames = find(db.Behavior_data.Behavior_data_table.Shelter_inout_Adjusted==0);
    Threat_off_Frames = find(db.Behavior_data.Behavior_data_table.Threat_onoff==0);
    Movement_Outside_Not_during_Threat_Frames = intersect(intersect(Movement_Frames,Shelter_out_Frames),Threat_off_Frames);
    
    Plot_WindRose_Goal_angle(db, Movement_Outside_Not_during_Threat_Frames)
    set(gcf,'position',[200 200 300 300])


%% Extract_ROI_Track
function ROI_Behavior_data_table = Extract_ROI_Track(db, ROI_time_range)
    ROI_Frames = Find_Frames_in_TS(ROI_time_range, db.Behavior_data.Behavior_data_table.Timestamp);
    ROI_Behavior_data_table = db.Behavior_data.Behavior_data_table(ROI_Frames,:);
end

%% Plot Traces w/ Goal angle
function Plot_Traces(db, ROI_time_range)
    Maze_X0Y0 = db.DAQ.Maze_X0Y0;
    Maze_XY_Coords = db.DAQ.Maze_XY_Coords;
    
    Maze_XY_Coords_z = zeros(size(Maze_XY_Coords,1),1);
    Maze_XYZ = [Maze_XY_Coords, Maze_XY_Coords_z];
    
    Maze_color = [220 220 220]/255;

    Shelter_X0Y0 = db.DAQ.Shelter_X0Y0;
    Shelter_XY_Coords = db.DAQ.Shelter_XY_Coords;   
    Shelter_color = [130 208 255]/255;
    % Shelter
    r = 5;   

    ROI_Frames = Find_Frames_in_TS(ROI_time_range, db.Behavior_data.Behavior_data_table.Timestamp);
    ROI_Behavior_data_table = db.Behavior_data.Behavior_data_table(ROI_Frames,:);
    Timestamp = ROI_Behavior_data_table.Timestamp;
    Coord_Data = ROI_Behavior_data_table.Body_Coords;
    ROI_Frames = 1:size(Coord_Data,1);
    End_point = db.DAQ.Shelter_X0Y0;
    Col_data = ROI_Behavior_data_table.Speed;
    Heading_Abs_Angle = ROI_Behavior_data_table.Heading_Abs_angle;
    Heading_Abs_Angle_Vec = Convert_Angle2CatesianCoords(Heading_Abs_Angle, 0);
        
    % Trajectory
    %% Plot Trajectory 2D
    x = Coord_Data(:,1);
    y = Coord_Data(:,2);
    
    [X,Y,Z] = cylinder(r,50);
    Shelter_cylinder_X = X+Shelter_X0Y0(1);
    Shelter_cylinder_Y = Y+Shelter_X0Y0(2);
    
    Diff_ROI_Coords = [diff(Coord_Data); nan nan];
    scale_factor1 = 1;
    Arrow1 = Diff_ROI_Coords.*scale_factor1;   
    Arrow2 = Heading_Abs_Angle_Vec;
    Arrow2_c = [6, 152 234]/255;
    figure(1), 
    set(gcf,'color','w','Position', [200, 200, 600, 600]);
    plot(Maze_XY_Coords(:,1), Maze_XY_Coords(:,2),'k'); hold on;
    plot(Shelter_X0Y0(:,1),Shelter_X0Y0(:,2),'bd'); hold on; 
    plot(Shelter_XY_Coords(:,1),Shelter_XY_Coords(:,2), 'color',Shelter_color); hold on;
    quiver(Coord_Data(:,1),Coord_Data(:,2),Arrow1(:,1),Arrow1(:,2),'autoscale','on','color',[0.5 0.5 0.5]); hold on;
    quiver(Coord_Data(:,1),Coord_Data(:,2),Arrow2(:,1),Arrow2(:,2),'color',Arrow2_c); hold on;    
    set(gca,'xtick',[],'ytick',[],'box','off','visible','on','ylim',[-0.1, 1.1],'tickdir','out');
    axis equal    
end


%% Plot Norm_dist2shelter w/ Heading angle
function Plot_Norm_dist2shelter(db, ROI_time_range)
    ROI_Frames = Find_Frames_in_TS(ROI_time_range, db.Behavior_data.Behavior_data_table.Timestamp);
    ROI_Behavior_data_table = db.Behavior_data.Behavior_data_table(ROI_Frames,:);
    Coord_Data = ROI_Behavior_data_table.Body_Coords;
    ROI_Frames = 1:size(Coord_Data,1);
    End_point = db.DAQ.Shelter_X0Y0;
    Col_data = ROI_Behavior_data_table.Speed;
    Heading_Abs_Angle = ROI_Behavior_data_table.Heading_Abs_angle;
    Ref_Angle = ROI_Behavior_data_table.Abs_angle;
    Heading_Relative_Angle_Vec = Convert_Angle2CatesianCoords(Heading_Abs_Angle-Ref_Angle(1),90);    

    ROI_Coord_Data = Coord_Data(ROI_Frames,:);
    Coord2Shelter = ROI_Coord_Data-End_point;
    Start_point = Coord2Shelter(1,:);

    % Transform
    Norm_Dist = norm(Start_point);
    theta = -(atan2d(Start_point(1,2),Start_point(1,1))+90);
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    ROI_Norm_Vector = Coord2Shelter./Norm_Dist;
    Rot_ROI_Norm_Vector = (R*ROI_Norm_Vector')';
    Rot_ROI_Norm_Vector(:,2) = Rot_ROI_Norm_Vector(:,2) + 1;


    Norm_ROI_Coord_Data = Rot_ROI_Norm_Vector;          
   
    Diff_ROI_Coords = [diff(Norm_ROI_Coord_Data); nan nan];
    scale_factor1 = 1;
    Arrow1 = Diff_ROI_Coords.*scale_factor1;   
    Arrow2 = Heading_Relative_Angle_Vec;
    Arrow2_c = [6 152 234]/255;
    figure(2) 
    plot(0,0,'bs'); hold on;
    plot(0,1,'kd'); hold on; 
    set(gcf,'color','w','Position', [300, 300, 300, 600]);
    quiver(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),Arrow1(:,1),Arrow1(:,2),'autoscale','on','color',[0.5 0.5 0.5]); hold on;
    quiver(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),Arrow2(:,1),Arrow2(:,2),'color',Arrow2_c); hold on;
    axis equal
    set(gca,'xtick',[],'ytick',[],'box','off','visible','on','ylim',[-0.1, 1.1],'tickdir','out');

end



