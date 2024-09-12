%% Fig4c
% Created by Kanghoon Jung
% Contact: kanghoon.jung AT alleninstitute.org for bug reports, questions, or feature requests.

close all 
clear all

data_folder = 'C:\Users\jungk\Google Drive\Github\Jung_NatNeuro_2024\codes\Fig4\';
cd(data_folder)
% Load the dataset from the .mat file
% load("Fig4c_data.mat");

   
T_Subject_Data_folder = 'H:\Project #5\Experiments\DREADD ret_fhM4Di@NAc CaMKII-Cre@vHPC\Behavior\hM4Di_CNO_First3';
% T_Subject_Data_folder = 'H:\Project #5\Experiments\DREADD ret_fhM4Di@NAc CaMKII-Cre@vHPC\Behavior\hM4Di_PBS_First3';
% T_Subject_Data_folder = 'H:\Project #5\Experiments\Behavior\Experiment 1_Normal flights\Data\New Dataset_Visual';

cd(T_Subject_Data_folder)    
Base_folder = T_Subject_Data_folder; 

%% subdirectory 
Sub_Folders = dir('*');
isub = [Sub_Folders(:).isdir];
Sub_Folders = Sub_Folders(isub);
Sub_Folders(ismember({Sub_Folders.name}',{'.','..'})) = [];
Num_Sub_Folders = size(Sub_Folders,1);

figure(2) 
set(gcf,'color','w','position',[200 200 200 500])
for(nSubj=1:Num_Sub_Folders)
% for(nSubj=1:1)    
    Subject_Data_folder = [Sub_Folders(nSubj,:).folder '\' Sub_Folders(nSubj,:).name];    
    cd(Subject_Data_folder)
    %% Import individual behavior dataset
    db_file = dir('db_beh*.mat');
    load(db_file.name)

    Behavior_data_T = db.Behavior_data.Behavior_data_table;
    Threat_on_time = 0;
    Threat_on_Frames = find(Behavior_data_T.Timestamp==Threat_on_time,1);
    Movement_on_Frame = find(Behavior_data_T.Movement_Adjusted==1,1);
    Flight_on_Frame = find(Behavior_data_T.Flight_Adjusted==1,1);
    Shelter_in_Frame = find(Behavior_data_T.Shelter_inout_Adjusted==1,1);

    ROI_on_time = []; %% Later one among Movement_on and Threat_on
    ROI_off_time = []; %% Sooner one among Shelter_in and Threat_off

    ROI_Frames = Behavior_data_T.Frame;

    ROI_on_time = Threat_on_time;
    if(~isempty(Shelter_in_Frame) & Shelter_in_Frame < ROI_Frames(end))
        ROI_off_time = Behavior_data_T.Timestamp(Shelter_in_Frame);
    else
        ROI_off_time = 10;
    end
    ROI_time_range = [ROI_on_time, ROI_off_time];       


%         Plot_Norm_dist2shelter1(db, [0 10])
    Plot_Norm_dist2shelter2(db, ROI_time_range)
    xlim([-0.7 1.2])
    ylim([-0.3 1.2])
%         axis equal


    %% Plot Normalized trace
%         [Cdb_t] = Plot_Stim(db);
%         [Cdb_t2] = Plot_Flight(db);
% 
%         for(i=1:size(Cdb_t,2))
%             Cdb{i}.x{1,k} = Cdb_t{i}.x;
%             Cdb{i}.y{1,k} = Cdb_t{i}.y;
%             Cdb2{i}.x{1,k} = Cdb_t2{i}.x;
%             Cdb2{i}.y{1,k} = Cdb_t2{i}.y;
%         end
%         k = k+1;
end
plot(0,0,'bs'); hold on;
plot(0,1,'ks'); hold on; 
    
    % end
% cd(T_Base_folder)

% % During Stim 
% T = Convert_CellData2Table(Cdb{1}); % Dist
% var.x = T.x_avg;
% var.y = T.mean;
% var.y_sem = T.sem;
% Plot_Dist2Shelter(var)
% % 
% T2 = Convert_CellData2Table(Cdb{2}); % Norm Dist
% var.x = T2.x_avg;
% var.y = T2.mean;
% var.y_sem = T2.sem;
% Plot_Dist2Shelter(var)
% % 
% % During Flight
% T = Convert_CellData2Table(Cdb2{1}); % Dist
% var.x = T.x_avg;
% var.y = T.mean;
% var.y_sem = T.sem;
% Plot_Dist2Shelter(var)
% % 
% T2 = Convert_CellData2Table(Cdb2{2}); % Norm Dist
% var.x = T2.x_avg;
% var.y = T2.mean;
% var.y_sem = T2.sem;
% Plot_Dist2Shelter(var)


% T = Convert_CellData2Table(Cdb{6}); % Norm Trace
% Plot_NormTrace(var)
% Base_folder = 'H:\Project #5\Experiments\Behavior\Experiment 1_Normal flights\Data\New Dataset\kjwt28_092220\kjwt28_092220_T3_AHO';
% db_file = dir('db_beh*.mat');
% load(db_file.name)
% Plot_Flight(db)

%% Plot Norm_dist2shelter w/ Heading angle
function Plot_Norm_dist2shelter1(db, ROI_time_range)
    ROI_Frames = Find_Frames_in_TS(ROI_time_range, db.Behavior_data.Behavior_data_table.Timestamp);
    ROI_Behavior_data_table = db.Behavior_data.Behavior_data_table(ROI_Frames,:);
    Coord_Data = ROI_Behavior_data_table.Body_Coords;
    ROI_Frames = 1:size(Coord_Data,1);
    End_point = db.DAQ.Shelter_X0Y0;
    Col_data = ROI_Behavior_data_table.Speed;
    Heading_Abs_Angle = ROI_Behavior_data_table.Heading_Abs_angle;
    Ref_Angle = ROI_Behavior_data_table.Abs_angle;
    Heading_Relative_Angle_Vec = Convert_Angle2CatesianCoords(Heading_Abs_Angle-Ref_Angle(1),90);
    
%     Goal_Angle = ROI_Behavior_data_table.Goal_Angle;
%     Goal_Angle_Vec = Convert_Angle2CatesianCoords(Goal_Angle, 90);
    
    
    ROI_Coord_Data = Coord_Data(ROI_Frames,:);
    Coord2Shelter = ROI_Coord_Data-End_point;
    Start_point = Coord2Shelter(1,:);

    % Transform
    Norm_Dist = norm(Start_point);
    theta = -(atan2d(Start_point(1,2),Start_point(1,1))+90);
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    ROI_Norm_Vector = Coord2Shelter./Norm_Dist;
    % ROI_Norm_Shelter = Shelter_XY_Coords./Norm_Dist;
    Rot_ROI_Norm_Vector = (R*ROI_Norm_Vector')';
    Rot_ROI_Norm_Vector(:,2) = Rot_ROI_Norm_Vector(:,2) + 1;


    Norm_ROI_Coord_Data = Rot_ROI_Norm_Vector;          
     % %% Plot Normalized distance
    % figure, 
    % set(gcf,'color','w','Position', [300, 300, 761, 569]);
    % plot(Rot_ROI_Norm_Vector(:,1),Rot_ROI_Norm_Vector(:,2),'g'); hold on;
%     figure, 
%     col = mat2gray(Col_data(ROI_Frames));  % This is the color, vary with x in this case.
%     scatter(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),25,col,'filled'); hold on;
%     plot(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),'color',[0.5 0.5 0.5]); hold on;
%     % plot(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2)); hold on;
%     caxis([min(col) max(col)])
%     colormap(gca, 'redblue')
%     axis equal
%     set(gca,'xtick',[],'ytick',[],'box','off','visible','on','ylim',[-0.1, 1.1],'tickdir','out');
%     colorbar;    
    
    Diff_ROI_Coords = [diff(Norm_ROI_Coord_Data); nan nan];
    scale_factor1 = 1;
    Arrow1 = Diff_ROI_Coords.*scale_factor1;   
    Arrow2 = Heading_Relative_Angle_Vec;
    Arrow2_c = [6 152 234]/255;

%     set(gcf,'color','w','Position', [300, 300, 300, 600]);
%     h1 = quiver(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),Arrow1(:,1),Arrow1(:,2),'autoscale','on','color',[0.5 0.5 0.5]); hold on;
%     h2 = quiver(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),Arrow2(:,1),Arrow2(:,2),'color',Arrow2_c); hold on;
%     h1.Color(4) = 0.3;
%     h2.Color(4) = 0.3;

    %     alpha(h1,0.1)
%     alpha(h2,0.1)
%     figure,
%     p = plot(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),'color',[0.2 0.2 0.2],'linewidth',1); hold on;
    p = plot(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),'color',[0.0235 .8 .95],'linewidth',0.4); hold on;
    p = plot(Norm_ROI_Coord_Data(end,1),Norm_ROI_Coord_Data(end,2),'o','color','b','linewidth',0.4); hold on;

    p.Color(4) = 0.3;
    set(gca,'box','off','visible','on','tickdir','out');
    
%     colorbar;    

%     [Norm_ROI_Coord_Data] = Norm_movement_trace_func(Coord_Data, ROI_Frames, End_point, Col_data); 
end

function Plot_Norm_dist2shelter2(db, ROI_time_range)
    ROI_Frames = Find_Frames_in_TS(ROI_time_range, db.Behavior_data.Behavior_data_table.Timestamp);
    ROI_Behavior_data_table = db.Behavior_data.Behavior_data_table(ROI_Frames,:);
    Coord_Data = ROI_Behavior_data_table.Body_Coords;
    ROI_Frames = 1:size(Coord_Data,1);
    End_point = db.DAQ.Shelter_X0Y0;
    Col_data = ROI_Behavior_data_table.Speed;
    Heading_Abs_Angle = ROI_Behavior_data_table.Heading_Abs_angle;
    Ref_Angle = ROI_Behavior_data_table.Abs_angle;
    Heading_Relative_Angle_Vec = Convert_Angle2CatesianCoords(Heading_Abs_Angle-Ref_Angle(1),90);
    
%     Goal_Angle = ROI_Behavior_data_table.Goal_Angle;
%     Goal_Angle_Vec = Convert_Angle2CatesianCoords(Goal_Angle, 90);
    
    
    ROI_Coord_Data = Coord_Data(ROI_Frames,:);
    Coord2Shelter = ROI_Coord_Data-End_point;
    Start_point = Coord2Shelter(1,:);

    % Transform
    Norm_Dist = norm(Start_point);
    theta = -(atan2d(Start_point(1,2),Start_point(1,1))+90);
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    ROI_Norm_Vector = Coord2Shelter./Norm_Dist;
    % ROI_Norm_Shelter = Shelter_XY_Coords./Norm_Dist;
    Rot_ROI_Norm_Vector = (R*ROI_Norm_Vector')';
    Rot_ROI_Norm_Vector(:,2) = Rot_ROI_Norm_Vector(:,2) + 1;


    Norm_ROI_Coord_Data = Rot_ROI_Norm_Vector;          
     % %% Plot Normalized distance
    % figure, 
    % set(gcf,'color','w','Position', [300, 300, 761, 569]);
    % plot(Rot_ROI_Norm_Vector(:,1),Rot_ROI_Norm_Vector(:,2),'g'); hold on;
%     figure, 
%     col = mat2gray(Col_data(ROI_Frames));  % This is the color, vary with x in this case.
%     scatter(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),25,col,'filled'); hold on;
%     plot(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),'color',[0.5 0.5 0.5]); hold on;
%     % plot(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2)); hold on;
%     caxis([min(col) max(col)])
%     colormap(gca, 'redblue')
%     axis equal
%     set(gca,'xtick',[],'ytick',[],'box','off','visible','on','ylim',[-0.1, 1.1],'tickdir','out');
%     colorbar;    
    
    Diff_ROI_Coords = [diff(Norm_ROI_Coord_Data); nan nan];
    scale_factor1 = 1;
    Arrow1 = Diff_ROI_Coords.*scale_factor1;   
    Arrow2 = Heading_Relative_Angle_Vec;
    Arrow2_c = [6 152 234]/255;

%     set(gcf,'color','w','Position', [300, 300, 300, 600]);
%     h1 = quiver(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),Arrow1(:,1),Arrow1(:,2),'autoscale','on','color',[0.5 0.5 0.5]); hold on;
%     h2 = quiver(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),Arrow2(:,1),Arrow2(:,2),'color',Arrow2_c); hold on;
%     h1.Color(4) = 0.3;
%     h2.Color(4) = 0.3;

    %     alpha(h1,0.1)
%     alpha(h2,0.1)
%     figure,
%     p = plot(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),'color',[0.2 0.2 0.2],'linewidth',1); hold on;
    p = plot(Norm_ROI_Coord_Data(:,1),Norm_ROI_Coord_Data(:,2),'color',Arrow2_c,'linewidth',2); hold on;
    p = plot(Norm_ROI_Coord_Data(end,1),Norm_ROI_Coord_Data(end,2),'.','color','b','linewidth',2); hold on;

    p.Color(4) = 0.3;
    set(gca,'box','off','visible','on','tickdir','out');
    
%     colorbar;    

%     [Norm_ROI_Coord_Data] = Norm_movement_trace_func(Coord_Data, ROI_Frames, End_point, Col_data); 
end




function [Cdb] = Plot_Stim(db)
% Key behavior datasets
% Tracking
Maze_X0Y0 = db.DAQ.Maze_X0Y0;
Shelter_X0Y0 = db.DAQ.Shelter_X0Y0;
Maze_XY_Coords = db.DAQ.Maze_XY_Coords;
Track_T = db.Behavior_data.Behavior_data_table;
Timestamp = Track_T.Timestamp;
Video_Timestamp = Timestamp;
Body_Coords = Track_T.Body_Coords;
Speed = Track_T.Speed;
Acc = Track_T.Acc;
Body_Coords2Shelter = Track_T.Body_Coords2Shelter;
Dist2Shelter = Track_T.Dist2Shelter;
Goal_Angle = Track_T.Goal_Angle;
S = db.Behavior_data.Analyzed.S;

i = 1;
Shelter_in_time = S{1, 1}.ROI_Track_T.ROI_Timestamp_adjusted(find(S{1, 1}.ROI_Track_T.Shelter_inout_Adjusted==1,1),1);
Shelter_in_Dist2Shelter = S{1, 1}.ROI_Track_T.Dist2Shelter(find(S{1, 1}.ROI_Track_T.Shelter_inout_Adjusted==1,1),1);
%% Plot Flight trace with speed color coding
% Draw Arena and Shelter
figure(1),
set(gcf,'color','w','Position', [200, 200, 600, 600]);
imshow(db.DAQ.Video_Info.I,db.DAQ.Video_Info.RA.cm); set(gca,'YDir','normal'); hold on;
plot(Maze_XY_Coords(:,1), Maze_XY_Coords(:,2),'y'); hold on;
plot(Shelter_X0Y0(1), Shelter_X0Y0(2),'bo'); hold on;
% Color coding
x = Body_Coords(:,1)';
y = Body_Coords(:,2)';
z = zeros(size(x));
% col = mat2gray(Speed)'; % This is the color, vary with x in this case.
col = Speed'; % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',5); 
set(gca,'box','off','TickDir','out');
% Colorbar
% caxis([min(col) 1])
caxis([min(Speed), max(Speed)])
% caxis([min(col) 0.0227])
c = colorbar('TickLength',0);
c.Label.String = 'Speed (cm/s)';
c.Label.FontSize = 12;                    
% set(c,'Ytick',[min(Speed), max(Speed)]); % Colorbar tick
colormap(viridis)


%% Plot Normalized Traces    
figure(2)
set(gcf,'color','w','Position', [3117, 250, 200, 600]);
plot(0,0,'r.'); hold on;    
plot(0,1,'c.'); hold on; 
axis auto; 
x = S{i,1}.Norm_Trace.Rot_ROI_Norm_Vector(:,1)';
y = S{i,1}.Norm_Trace.Rot_ROI_Norm_Vector(:,2)';
z = zeros(size(x));

Cdb{6}.x = x;
Cdb{6}.y = y;
% col = mat2gray(Speed)'; % This is the color, vary with x in this case.
col = Speed'; % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',5); 
ylim([-0.2, 1.2])    
xlim([-1, 1])    
c_range = [0, max(Speed)];
caxis(c_range)
c = colorbar('location','EastOutside','Ticks',floor(c_range),'TickLength',0);
c.Label.String = 'Speed (cm/s)';
c.Label.FontSize = 10;  
cpos = c.Position; 
cpos(2) = 1.1*cpos(2);
cpos(4) = 0.5*cpos(4);
c.Position = cpos;
colormap(viridis)
set(gca,'box','off','TickDir','out');
     
%% Plot Dist2Shelter
Timestamp = Track_T.Timestamp;
Dist2Shelter = Track_T.Dist2Shelter;
x = S{i,1}.ROI_Track_T.ROI_Timestamp_adjusted;
y = S{i,1}.ROI_Track_T.Dist2Shelter;
x_line = Shelter_in_time;

figure(3)
set(gcf,'color','w','Position', [2110, 607, 450, 370]);
plot(x,y); hold on;
if(~isempty(Shelter_in_time))
    y_line = Shelter_in_Dist2Shelter;
    plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
    yl = ylim;
    plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
set(gca,'box','off','tickdir','out')
xlabel('Time (s)')
ylabel('Dist2Shelter')

Cdb{1}.x = x;
Cdb{1}.y = y;

%% Plot Norm Dist2Shelter
i=1;
% Norm Dist
Norm_Dist2Shelter = Dist2Shelter./Dist2Shelter(find(~isnan(Dist2Shelter),1));
x = Timestamp;
y = Norm_Dist2Shelter;

% y_line = Shelter_in_Dist2Shelter;
figure(4)
set(gcf,'color','w','Position', [2600, 607, 450, 370]);
plot(x,y); hold on;
% plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
set(gca,'ylim',[0, 1.2]);
if(~isempty(Shelter_in_time))
    x_line = Shelter_in_time;
    yl = ylim;
    plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
set(gca,'box','off','tickdir','out')
xlabel('Time (s)')
ylabel('Norm Dist2Shelter')

Cdb{2}.x = x;
Cdb{2}.y = y;

%% Plot Speed
figure(5)
set(gcf,'color','w','Position', [2110, 147, 450, 370]);
x = Timestamp;
y = Speed;
y_line = 10;
plot(x,y); hold on;
% plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
yl = ylim;
if(~isempty(Shelter_in_time))
    plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
set(gca,'box','off','tickdir','out')
xlabel('Time (s)')
ylabel('Speed (m/s)')

Cdb{3}.x = x;
Cdb{3}.y = y;

%% Plot Head-Shelter Angle with respect to time
figure(6)
set(gcf,'color','w','Position', [2110, -312, 450, 370]);
x = Timestamp;
y = cosd(Goal_Angle); 
% y_line = 10;
plot(x,y); hold on;
yl = ylim;
if(~isempty(Shelter_in_time))
    plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
% plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
set(gca,'box','off','tickdir','out')
xlabel('Time (s)')
ylabel('SI')


Cdb{4}.x = x;
Cdb{4}.y = y;


%% Plot Head-Shelter Angle with respect to Dist2Shelter
figure(7)
set(gcf,'color','w','Position', [2600, 147, 450, 370]);
x = -S{i,1}.ROI_Track_T.Dist2Shelter;
y = cosd(S{i,1}.ROI_Track_T.Goal_Angle); 
% y_line = 10;
plot(x,y); hold on;
yl = ylim;
xline = Shelter_in_Dist2Shelter;
if(~isempty(Shelter_in_time))
    plot(-[x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
% plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
set(gca,'box','off','tickdir','out')
xlabel('Dist2Shelter (cm)')
ylabel('SI')


Cdb{5}.x = x;
Cdb{5}.y = y;

%% Plot Dist2shelter and Flight displacement
figure(8)
set(gcf,'color','w','Position', [2600, -312, 450, 370]);
plot(linspace(0,70),linspace(0,70));  hold on;
s = scatter(S{i,1}.Dist2Shelter_Displacement(1), S{i,1}.Dist2Shelter_Displacement(2),20,'k','filled'); hold on;
% %     elseif(Events(i,2)==4)
% %         s = scatter(S{i,1}.Dist2Shelter_Displacement(1), S{i,1}.Dist2Shelter_Displacement(2),20,'c','filled'); hold on;
% %     end
s.Marker = 'o';
set(gca,'box','off','tickdir','out')
xlabel('Dist2Shelter');
ylabel('Flight Displacement');
end


function [Cdb] = Plot_Flight(db)
% Key behavior datasets
% Tracking
Maze_X0Y0 = db.DAQ.Maze_X0Y0;
Shelter_X0Y0 = db.DAQ.Shelter_X0Y0;
Maze_XY_Coords = db.DAQ.Maze_XY_Coords;
S = db.Behavior_data.Analyzed.S;

i = 1;
ROI_Track_T = S{i, 1}.ROI_Track_T;
Flight_onset_frame = find(ROI_Track_T.Flight_Adjusted==1,1);
Flight_offset_frame = find(ROI_Track_T.Flight_Adjusted==1,1,'last');

if(~isempty(Flight_onset_frame))
    Onset_frame = Flight_onset_frame;
    Flight_Track_T = ROI_Track_T(Onset_frame:end,:);
    Flight_Track_T.ROI_Frames =  Flight_Track_T.Frame - Flight_Track_T.Frame(1)+1;
    Flight_Track_T.ROI_Timestamp_adjusted  =  Flight_Track_T.Timestamp - Flight_Track_T.Timestamp(1);

    Shelter_in_frame = Flight_Track_T.ROI_Frames(find(Flight_Track_T.Shelter_inout_Adjusted==1,1),1);
    % Shelter_in_Dist2Shelter = S{1, 1}.ROI_Track_T.Dist2Shelter(find(S{1, 1}.ROI_Track_T.Shelter_inout_Adjusted==1,1),1);
    if(~isempty(Shelter_in_frame))
        Offset_frame = Shelter_in_frame;
        Flight_Track_T = Flight_Track_T(1:Offset_frame,:);
    end
Shelter_in_time = Flight_Track_T.Timestamp(Shelter_in_frame);
Shelter_in_Dist2Shelter = Flight_Track_T.Dist2Shelter(Shelter_in_frame);
% Track_T = db.Behavior_data.Behavior_data_table;
Timestamp = Flight_Track_T.Timestamp;
Video_Timestamp = Timestamp;
Body_Coords = Flight_Track_T.Body_Coords;
Speed = Flight_Track_T.Speed;
Acc = Flight_Track_T.Acc;
Body_Coords2Shelter = Flight_Track_T.Body_Coords2Shelter;
Dist2Shelter = Flight_Track_T.Dist2Shelter;
Goal_Angle = Flight_Track_T.Goal_Angle;



%% Plot Flight trace with speed color coding
% Draw Arena and Shelter
figure(1),
set(gcf,'color','w','Position', [200, 200, 600, 600]);
imshow(db.DAQ.Video_Info.I,db.DAQ.Video_Info.RA.cm); set(gca,'YDir','normal'); hold on;
plot(Maze_XY_Coords(:,1), Maze_XY_Coords(:,2),'y'); hold on;
plot(Shelter_X0Y0(1), Shelter_X0Y0(2),'bo'); hold on;
% Color coding
x = Body_Coords(:,1)';
y = Body_Coords(:,2)';
z = zeros(size(x));
% col = mat2gray(Speed)'; % This is the color, vary with x in this case.
col = Speed'; % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',5); 
set(gca,'box','off','TickDir','out');
% Colorbar
% caxis([min(col) 1])
caxis([min(Speed), max(Speed)])
% caxis([min(col) 0.0227])
c = colorbar('TickLength',0);
c.Label.String = 'Speed (cm/s)';
c.Label.FontSize = 12;                    
% set(c,'Ytick',[min(Speed), max(Speed)]); % Colorbar tick
colormap(viridis)


%% Plot Normalized Traces    
figure(2)
set(gcf,'color','w','Position', [3117, 250, 200, 600]);
plot(0,0,'r.'); hold on;    
plot(0,1,'c.'); hold on; 
axis auto; 
x = S{i,1}.Norm_Trace.Rot_ROI_Norm_Vector(:,1)';
y = S{i,1}.Norm_Trace.Rot_ROI_Norm_Vector(:,2)';
z = zeros(size(x));

Cdb{6}.x = x;
Cdb{6}.y = y;
% col = mat2gray(Speed)'; % This is the color, vary with x in this case.
col = Speed'; % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',5); 
ylim([-0.2, 1.2])    
xlim([-1, 1])    
c_range = [0, max(Speed)];
caxis(c_range)
c = colorbar('location','EastOutside','Ticks',floor(c_range),'TickLength',0);
c.Label.String = 'Speed (cm/s)';
c.Label.FontSize = 10;  
cpos = c.Position; 
cpos(2) = 1.1*cpos(2);
cpos(4) = 0.5*cpos(4);
c.Position = cpos;
colormap(viridis)
set(gca,'box','off','TickDir','out');
     
%% Plot Dist2Shelter
% Timestamp = Track_T.Timestamp;
% Dist2Shelter = Track_T.Dist2Shelter;
x = Timestamp;
y = Dist2Shelter;
x_line = Shelter_in_time;

figure(3)
set(gcf,'color','w','Position', [2110, 607, 450, 370]);
plot(x,y); hold on;
if(~isempty(Shelter_in_time))
    y_line = Shelter_in_Dist2Shelter;
    plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
    yl = ylim;
    plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
set(gca,'box','off','tickdir','out')
xlabel('Time (s)')
ylabel('Dist2Shelter')

Cdb{1}.x = x;
Cdb{1}.y = y;

%% Plot Norm Dist2Shelter
i=1;
% Norm Dist
Norm_Dist2Shelter = Dist2Shelter./Dist2Shelter(find(~isnan(Dist2Shelter),1));
x = Timestamp;
y = Norm_Dist2Shelter;

% y_line = Shelter_in_Dist2Shelter;
figure(4)
set(gcf,'color','w','Position', [2600, 607, 450, 370]);
plot(x,y); hold on;
% plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
set(gca,'ylim',[0, 1.2]);
if(~isempty(Shelter_in_time))
    x_line = Shelter_in_time;
    yl = ylim;
    plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
set(gca,'box','off','tickdir','out')
xlabel('Time (s)')
ylabel('Norm Dist2Shelter')

Cdb{2}.x = x;
Cdb{2}.y = y;

%% Plot Speed
figure(5)
set(gcf,'color','w','Position', [2110, 147, 450, 370]);
x = Timestamp;
y = Speed;
y_line = 10;
plot(x,y); hold on;
% plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
yl = ylim;
if(~isempty(Shelter_in_time))
    plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
set(gca,'box','off','tickdir','out')
xlabel('Time (s)')
ylabel('Speed (m/s)')

Cdb{3}.x = x;
Cdb{3}.y = y;

%% Plot Head-Shelter Angle with respect to time
figure(6)
set(gcf,'color','w','Position', [2110, -312, 450, 370]);
x = Timestamp;
y = cosd(Goal_Angle); 
% y_line = 10;
plot(x,y); hold on;
yl = ylim;
if(~isempty(Shelter_in_time))
    plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
% plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
set(gca,'box','off','tickdir','out')
xlabel('Time (s)')
ylabel('SI')


Cdb{4}.x = x;
Cdb{4}.y = y;


%% Plot Head-Shelter Angle with respect to Dist2Shelter
figure(7)
set(gcf,'color','w','Position', [2600, 147, 450, 370]);
x = -S{i,1}.ROI_Track_T.Dist2Shelter;
y = cosd(S{i,1}.ROI_Track_T.Goal_Angle); 
% y_line = 10;
plot(x,y); hold on;
yl = ylim;
xline = Shelter_in_Dist2Shelter;
if(~isempty(Shelter_in_time))
    plot(-[x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
end
% plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
set(gca,'box','off','tickdir','out')
xlabel('Dist2Shelter (cm)')
ylabel('SI')


Cdb{5}.x = x;
Cdb{5}.y = y;

%% Plot Dist2shelter and Flight displacement
figure(8)
set(gcf,'color','w','Position', [2600, -312, 450, 370]);
plot(linspace(0,70),linspace(0,70));  hold on;
s = scatter(S{i,1}.Dist2Shelter_Displacement(1), S{i,1}.Dist2Shelter_Displacement(2),20,'k','filled'); hold on;
% %     elseif(Events(i,2)==4)
% %         s = scatter(S{i,1}.Dist2Shelter_Displacement(1), S{i,1}.Dist2Shelter_Displacement(2),20,'c','filled'); hold on;
% %     end
s.Marker = 'o';
set(gca,'box','off','tickdir','out')
xlabel('Dist2Shelter');
ylabel('Flight Displacement');
else
    for i=1:6
        Cdb{i}.x = [];
        Cdb{i}.y = [];
    end
end
end




