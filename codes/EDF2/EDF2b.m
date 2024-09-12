%% Fiber photometry homing analysis code
%% Code for EDF2b
% Created by Kanghoon Jung
% Contact: kanghoon.jung AT alleninstitute.org for bug reports, questions, or feature requests.

clear all 
close all

data_folder = 'C:\Users\jungk\Google Drive\Github\Jung_NatNeuro_2024\codes\EDF2\';
cd(data_folder)

% Load the dataset from the .mat file
load("EDF2b_data.mat");

%% Analysis
Maze_X0Y0 = db.DAQ.Maze_X0Y0;
Shelter_X0Y0 = db.DAQ.Shelter_X0Y0;
Maze_XY_Coords = db.DAQ.Maze_XY_Coords;
T = db.Matched_data.Data_table;

Shelter_inout_Adjusted = T.Shelter_inout_Adjusted;    
Timestamp = T.Synced_Photometry_Timestamp;
DFoF_G = T.DFoF_ChG;
T.zscore_DFoF_ChG = zscore(T.DFoF_ChG);
Dist2Shelter = T.Dist2Shelter;
Body_Coords = T.Body_Coords;
Speed = T.Speed;
Acc = T.Acc;
Body_Coords2Shelter = T.Body_Coords2Shelter;
Goal_Angle = T.Goal_Angle;  

minimum_frame_num = 15;
[Shelter_inout_Frames, Adjusted_States] = Convert_States2OnOffFrames(T.Shelter_inout,minimum_frame_num);
[Movement_onoff_Frames, Adjusted_States]  = Convert_States2OnOffFrames(T.Movement_Adjusted,minimum_frame_num);

Shelter_inout_times = Timestamp(Shelter_inout_Frames);    

Movement_states = T.Movement;
Shelter_states = T.Shelter_inout_Adjusted;

[Movement_Bout_Table] = Convert_States2Bouts(T.Movement);
[Shelter_Bout_Table] = Convert_States2Bouts(T.Shelter_inout_Adjusted);

% Status1
T.Status(1:size(T,1)) = {'NA'};
T.Status(find(Shelter_inout_Adjusted)) = {'Shelter'};
T.Status(find(~Shelter_inout_Adjusted)) = {'Exploration'};

% Status2      
for(i=1:size(Movement_Bout_Table,1))
    Bout_Frames = cell2mat(Movement_Bout_Table.Frames(i,1));
    Shelter_values = unique(Shelter_states(Bout_Frames));
    if(Movement_Bout_Table.Value(i,1)==0)
        if(all(ismember(unique(Shelter_values),1)))
            T.Status2(Bout_Frames) = {'Still-Shelter'};
        elseif(all(ismember(unique(Shelter_values),0)))
            T.Status2(Bout_Frames) = {'Still-Outside'};
        else
            T.Status2(Bout_Frames) = {'Undefined'};    
        end
    elseif(Movement_Bout_Table.Value(i,1)==1)
        if(all(ismember(unique(Shelter_values),1)))
            T.Status2(Bout_Frames) = {'Move-Shelter'};
        elseif(all(ismember(unique(Shelter_values),0)))
            T.Status2(Bout_Frames) = {'Move-Outside'};
        else            
            T.Status2(Bout_Frames) = {'Transition'};
        end
    end
end
    
% Status3
Transition_states = zeros(size(T,1),1);
Transition_states(find(contains(T.Status2,'Transition')),1) = 1;
[Transition_Bout_Table] = Convert_States2Bouts(Transition_states);

T.Status3(1:size(T,1)) = {'NA'};
for(i=1:size(Transition_Bout_Table,1))
    if(Transition_Bout_Table.Value(i,1)==1)
        Bout_Frames = cell2mat(Transition_Bout_Table.Frames(i,1));
        Init_values = T.Status{Bout_Frames(1)};
        End_values = T.Status{Bout_Frames(end)};
        if(contains(Init_values, 'Shelter') & contains(End_values, 'Shelter'))
            [M, I] = max(T.Dist2Shelter(Bout_Frames));
            Max_Dist2Shelter_Frame = Bout_Frames(I);
            T.Status3(Bout_Frames(1):Max_Dist2Shelter_Frame) = {'ToMaxDist'};
            T.Status3(Max_Dist2Shelter_Frame+1:Bout_Frames(end)) = {'FromMaxDist'};       
        elseif(contains(Init_values, 'Shelter') & contains(End_values, 'Exploration'))
            [M, I] = max(T.Dist2Shelter(Bout_Frames));
            Max_Dist2Shelter_Frame = Bout_Frames(I);
            T.Status3(Bout_Frames(1):Max_Dist2Shelter_Frame) = {'LeavingToMaxDist'};
            T.Status3(Max_Dist2Shelter_Frame+1:Bout_Frames(end)) = {'LeavingFromMaxDist'};       
        elseif(contains(Init_values, 'Exploration') & contains(End_values, 'Shelter'))
            [M, I] = max(T.Dist2Shelter(Bout_Frames));
            Max_Dist2Shelter_Frame = Bout_Frames(I);
            T.Status3(Bout_Frames(1):Max_Dist2Shelter_Frame) = {'HomingToMaxDist'};
            T.Status3(Max_Dist2Shelter_Frame+1:Bout_Frames(end)) = {'HomingFromMaxDist'};                   
        end 
    end
end

% 1) Leaving Shleter
Leaving_states = zeros(size(T,1),1); 
Leaving_states(find(contains(T.Status3,'ToMaxDist')|contains(T.Status3,'LeavingToMaxDist')),1) = 1;
[Leaving_Bout_Table] = Convert_States2Bouts(Leaving_states);
Leaving_Bout_Table = Leaving_Bout_Table(find(Leaving_Bout_Table.Value==1),:);

% Plot Leaving Movement traces with color codes
ROI_Table = Leaving_Bout_Table;
Num_movements = size(ROI_Table,1);
ROI_Leaving_trials = [];
figure(1), 
set(gcf,'color','w','position',[600, 100, 800, 800])
imshow(db.DAQ.Video_Info.I,db.DAQ.Video_Info.RA.cm); set(gca,'YDir','normal'); hold on;    
for(i=1:Num_movements)        
    ROI_Frames = ROI_Table.Frames{i,1};
    plot(Shelter_X0Y0(:,1), Shelter_X0Y0(:,2),'ko'); hold on;
            
   if(Dist2Shelter(ROI_Frames(1),1)<5 & Dist2Shelter(ROI_Frames(end),1)>15)
        ROI_Leaving_trials = [ROI_Leaving_trials; i];
        % Color coding
        x = Body_Coords(ROI_Frames,1)';
        y = Body_Coords(ROI_Frames,2)';
        z = zeros(size(x));
        col = T.zscore_DFoF_ChG(ROI_Frames)'; % This is the color, vary with x in this case.

        surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',3); hold on;
        c = colorbar('TickLength',0);
        c.Label.String = 'DFoF';
        c.Label.FontSize = 12;                    
        colormap("parula")
        drawnow    
        caxis([-1 1])
   end
    title('Leaving shelter DFoF')    
end
    
    
% 2) Homing
Homing_states = zeros(size(T,1),1); 
Homing_states(find(contains(T.Status3,'FromMaxDist')|contains(T.Status3,'HomingFromMaxDist')),1) = 1;
[Homing_Bout_Table] = Convert_States2Bouts(Homing_states);
Homing_Bout_Table = Homing_Bout_Table(find(Homing_Bout_Table.Value==1),:);

% Plot Homing Movement traces with color codes
ROI_Table = Homing_Bout_Table;
Num_movements = size(ROI_Table,1);
ROI_Homing_trials = [];
figure(2), 
set(gcf,'color','w','position',[600, 100, 800, 800])
imshow(db.DAQ.Video_Info.I,db.DAQ.Video_Info.RA.cm); set(gca,'YDir','normal'); hold on;    
plot(Shelter_X0Y0(:,1), Shelter_X0Y0(:,2),'ko'); hold on;
for(i=1:Num_movements)        
    ROI_Frames = ROI_Table.Frames{i,1};
   if(Dist2Shelter(ROI_Frames(1),1)>10 & Dist2Shelter(ROI_Frames(end),1)<5)
       if(~isempty(find(T.Threat_onoff(ROI_Frames)>0)))
       ROI_Homing_trials = [ROI_Homing_trials; i];
        % Color coding
        x = Body_Coords(ROI_Frames,1)';
        y = Body_Coords(ROI_Frames,2)';
        z = zeros(size(x));
        col = zscore(T.zscore_DFoF_ChG(ROI_Frames))'; % This is the color, vary with x in this case.

        surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',3); hold on;
        c = colorbar('TickLength',0);
        c.Label.String = 'DFoF';
        c.Label.FontSize = 12;                   

        colormap("parula")
        drawnow    
        caxis([-1 1]);
       end
   end
    title('Homing DFoF')    
end
    

















