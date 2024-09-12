%% Shelter preference
% Created by Kanghoon Jung
% Contact: kanghoon.jung AT alleninstitute.org for bug reports, questions, or feature requests.

function [db] = GP_shelter_preference_func(db)
Shelter_In_vs_Out_stats.In_vs_Out = [];
Shelter_In_vs_Out_stats.In_vs_Rand = [];
Shelter_In_vs_Out_stats.Out_vs_Rand = [];

Field_name = ['GP_shelter_preference']

Sampling_rate = db.Matched_data.Matched_Neural_data.Sampling_rate;
Raster_DFoF = db.Matched_data.Matched_Neural_data.Raster_DFoF;    
ROI_cell_F = Raster_DFoF*Sampling_rate;
Full_ROI_cell_F = ROI_cell_F;


if(~isempty(db.Matched_data.Matched_Behavior_data.Shelter_in_out_time))
    Shelter_in_out_time = db.Matched_data.Matched_Behavior_data.Shelter_in_out_time;
    Raster_DFoF = db.Matched_data.Matched_Neural_data.Raster_DFoF;    

    NFrames = size(Full_ROI_cell_F,1);
    ROI_NFrames = size(ROI_cell_F,1);
    ncells = size(ROI_cell_F,2);
    %% Shelter_State
    Shelter_States = zeros(NFrames,1);
    TS = db.Matched_data.Matched_Neural_data.Timestamp;
    for i=1:size(Shelter_in_out_time,1)    
        t_range = [Shelter_in_out_time(i,1) Shelter_in_out_time(i,2)];
        ROI_Frames = Find_Frames_in_TS(t_range, TS);
        Shelter_States(ROI_Frames,1) = 1;
    end

else
    Shelter_States = db.Matched_data.Matched_Track_T.Shelter_inout_Adjusted;
end
    
Shelter_In_Frames = find(Shelter_States(:,1)==1);
Shelter_Out_Frames = find(Shelter_States(:,1)==0);   
Shelter_In_Frames(find(Shelter_In_Frames>ROI_NFrames)) = [];
Shelter_Out_Frames(find(Shelter_Out_Frames>ROI_NFrames)) = [];

% Time in Shelter 
Time_Shelter = numel(Shelter_In_Frames)/Sampling_rate;

% Shuffled frames;
Shuffled_Frames = cell(1,1);
Shuffled_Shelter_In_Frames = cell(1,1);
Shuffled_Shelter_Out_Frames = cell(1,1);
Avg_ER_Shuffled_Shelter_In_Frames = [];
Avg_ER_Shuffled_Shelter_Out_Frames = [];
nSim = 100;
for nS = 1:nSim
    Shuffled_Frames{nS,1} = randperm(ROI_NFrames);
    Shuffled_Shelter_In_Frames{nS,1} = Shuffled_Frames{nS,1}(Shelter_In_Frames);
    Shuffled_Shelter_Out_Frames{nS,1} = Shuffled_Frames{nS,1}(Shelter_Out_Frames);
    Avg_ER_Shuffled_Shelter_In_Frames(nS,:) = mean(ROI_cell_F(Shuffled_Shelter_In_Frames{nS,1},:));    
    Avg_ER_Shuffled_Shelter_Out_Frames(nS,:) = mean(ROI_cell_F(Shuffled_Shelter_Out_Frames{nS,1},:));    
end

% Shelter preferred firing test (In vs out)
Shelter_In_vs_Out_stats.In_vs_Out =[];
Avg_ER_Shelter_In_Frames = [];
Avg_ER_Shelter_Out_Frames = [];
for nC=1:ncells
    Shelter_In_vs_Out_stats.In_vs_Out(1:2,nC) = [mean(ROI_cell_F(Shelter_In_Frames,nC)); mean(ROI_cell_F(Shelter_Out_Frames,nC))];
    Avg_ER_Shelter_In_Frames(nC,1) = mean(ROI_cell_F(Shelter_In_Frames,nC));
    Avg_ER_Shelter_Out_Frames(nC,1) = mean(ROI_cell_F(Shelter_Out_Frames,nC));
    % In vs Out
    [h1 pvalue1] = ttest2(ROI_cell_F(Shelter_In_Frames,nC),ROI_cell_F(Shelter_Out_Frames,nC));
    Shelter_In_vs_Out_stats.In_vs_Out(3:4,nC) = [h1; pvalue1];
end    
Shelter_preference_ratio_index = log(Avg_ER_Shelter_In_Frames./Avg_ER_Shelter_Out_Frames);

Shelter_Neurons = [find(((Shelter_In_vs_Out_stats.In_vs_Out(1,:)-Shelter_In_vs_Out_stats.In_vs_Out(2,:))>0) & (Shelter_In_vs_Out_stats.In_vs_Out(3,:)==1))]';    
Outside_Neurons = [find(((Shelter_In_vs_Out_stats.In_vs_Out(1,:)-Shelter_In_vs_Out_stats.In_vs_Out(2,:))<0) & (Shelter_In_vs_Out_stats.In_vs_Out(3,:)==1))]';    
Nonspecific_Neurons = [find(((Shelter_In_vs_Out_stats.In_vs_Out(1,:)-Shelter_In_vs_Out_stats.In_vs_Out(2,:))==0) | (Shelter_In_vs_Out_stats.In_vs_Out(3,:)==0))]';

% Significance test using shuffled test
Stats_In_Frames_Emp_vs_Rand = [];
Stats_Out_Frames_Emp_vs_Rand = [];
for nC=1:ncells
    [Stats_In_Frames_Emp_vs_Rand(nC,1) Stats_In_Frames_Emp_vs_Rand(nC,2)] = ttest2(Avg_ER_Shelter_In_Frames(nC,1),Avg_ER_Shuffled_Shelter_In_Frames(:,nC));
    [Stats_Out_Frames_Emp_vs_Rand(nC,1) Stats_Out_Frames_Emp_vs_Rand(nC,2)] = ttest2(Avg_ER_Shelter_Out_Frames(nC,1),Avg_ER_Shuffled_Shelter_Out_Frames(:,nC));
end
True_Shelter_Neurons = intersect(Shelter_Neurons,find(Stats_In_Frames_Emp_vs_Rand==1));
True_Outside_Neurons = intersect(Outside_Neurons,find(Stats_Out_Frames_Emp_vs_Rand==1));

Shelter_Neurons_Population = [numel(Shelter_Neurons) numel(True_Shelter_Neurons) numel(Outside_Neurons) numel(True_Outside_Neurons) numel(Nonspecific_Neurons) ncells]

%% Data to save
db.Analysis.(Field_name).ROI_cell_F = ROI_cell_F;
db.Analysis.(Field_name).Shelter_States = Shelter_States;
db.Analysis.(Field_name).Shelter_In_Frames = Shelter_In_Frames;
db.Analysis.(Field_name).Shelter_Out_Frames = Shelter_Out_Frames;   
db.Analysis.(Field_name).Avg_ER_Shelter_In_Frames = Avg_ER_Shelter_In_Frames;
db.Analysis.(Field_name).Avg_ER_Shelter_Out_Frames = Avg_ER_Shelter_Out_Frames;
db.Analysis.(Field_name).Avg_ER_Shuffled_Shelter_In_Frames = Avg_ER_Shuffled_Shelter_In_Frames;
db.Analysis.(Field_name).Avg_ER_Shuffled_Shelter_Out_Frames = Avg_ER_Shuffled_Shelter_Out_Frames;
db.Analysis.(Field_name).Shelter_preference_ratio_index = Shelter_preference_ratio_index;
db.Analysis.(Field_name).Shelter_In_vs_Out_stats.In_vs_Out = Shelter_In_vs_Out_stats.In_vs_Out;
db.Analysis.(Field_name).Shelter_Neurons = Shelter_Neurons;
db.Analysis.(Field_name).Outside_Neurons = Outside_Neurons;
db.Analysis.(Field_name).Nonspecific_Neurons = Nonspecific_Neurons;
db.Analysis.(Field_name).True_Shelter_Neurons = True_Shelter_Neurons;
db.Analysis.(Field_name).True_Outside_Neurons = True_Outside_Neurons;
db.Analysis.(Field_name).Shelter_Neurons_Population = Shelter_Neurons_Population;       
        
        
% %% Draw Avg. ER rate In vs Out
% figure, 
% set(gcf,'color','w')
% set(gcf,'Position', [100, 50, 500, 500]); % x, y, width, height
% scatter(Avg_ER_Shelter_In_Frames, Avg_ER_Shelter_Out_Frames,40,'MarkerEdgeColor','none', 'MarkerFaceColor',[0.5 .5 .5]); hold on;
% scatter(Avg_ER_Shelter_In_Frames(True_Shelter_Neurons),Avg_ER_Shelter_Out_Frames(True_Shelter_Neurons),40,'MarkerEdgeColor','none', 'MarkerFaceColor','g'); hold on;
% scatter(Avg_ER_Shelter_In_Frames(True_Outside_Neurons),Avg_ER_Shelter_Out_Frames(True_Outside_Neurons),40,'MarkerEdgeColor','none', 'MarkerFaceColor','r'); hold on;
% % axis equal    
% set(gca,'tickdir','out','xtick',[0:0.1:0.3],'ytick',[0:0.1:0.3],'xlim',[0, 0.35],'ylim',[0, 0.35]);
% 

