%% Code for Fig2l and EDF5j
% Created by Kanghoon Jung
% Contact: kanghoon.jung AT alleninstitute.org for bug reports, questions, or feature requests.

clear all
close all
data_folder = 'C:\Users\jungk\Google Drive\Github\Jung_NatNeuro_2024\codes\Fig2\';
cd(data_folder)

% Load the dataset from the .mat file
load("Fig2l_EDF5j_data.mat");

% Create Table from the struct data
db.Matched_data.Matched_T.B_Timestamp = db.Matched_data.Matched_Behavior_data.Timestamp;
db.Matched_data.Matched_T.F_Timestamp = db.Matched_data.Matched_Neural_data.Timestamp;
db.Matched_data.Matched_T.DFoF = db.Matched_data.Matched_Neural_data.DFoF;
db.Matched_data.Matched_T.Raster_DFoF = db.Matched_data.Matched_Neural_data.Raster_DFoF;
T = [struct2table(db.Matched_data.Matched_Behavior_data.Data), struct2table(db.Matched_data.Matched_T)];

%% Shelter preference Analysis
Field_name = ['GP_shelter_preference'];
[db] = GP_shelter_preference_func(db);

Avg_ER_Shelter_In_Frames = db.Analysis.(Field_name).Avg_ER_Shelter_In_Frames;
Avg_ER_Shelter_Out_Frames = db.Analysis.(Field_name).Avg_ER_Shelter_Out_Frames;
Shelter_preference_ratio_index = log(Avg_ER_Shelter_In_Frames./Avg_ER_Shelter_Out_Frames);
Shelter_Neurons = db.Analysis.(Field_name).Shelter_Neurons;
Outside_Neurons = db.Analysis.(Field_name).Outside_Neurons;
Nonspecific_Neurons = db.Analysis.(Field_name).Nonspecific_Neurons;
True_Shelter_Neurons = db.Analysis.(Field_name).True_Shelter_Neurons;
True_Outside_Neurons = db.Analysis.(Field_name).True_Outside_Neurons;
Shelter_Neurons_Population = db.Analysis.(Field_name).Shelter_Neurons_Population;    
disp(['Ncells: ' num2str(Shelter_Neurons_Population(1,6))])

%Find Right Shelter Neurons
Shelter_Neurons = db.Analysis.GP_shelter_preference.Shelter_Neurons;
Outside_Neurons = db.Analysis.GP_shelter_preference.Outside_Neurons;
Nonspecific_Neurons = db.Analysis.GP_shelter_preference.Nonspecific_Neurons;
Shelter_Shelter_preference_ratio = db.Analysis.GP_shelter_preference.Shelter_preference_ratio_index(Shelter_Neurons);         
Outside_Shelter_preference_ratio = db.Analysis.GP_shelter_preference.Shelter_preference_ratio_index(Outside_Neurons);         
Nonspecific_Shelter_preference_ratio = db.Analysis.GP_shelter_preference.Shelter_preference_ratio_index(Nonspecific_Neurons);         

[Sort_Shelter_Shelter_preference_ratio, Shelter_I] = sort(Shelter_Shelter_preference_ratio,'descend');
[Sort_Outside_Shelter_preference_ratio, Outside_I] = sort(Outside_Shelter_preference_ratio,'descend');
[Sort_Nonspecific_Shelter_preference_ratio, Nonspecific_I] = sort(Nonspecific_Shelter_preference_ratio,'descend');

DFoF = T.DFoF;
zscore_DFoF = zscore(DFoF,0,1);
Dist2Shelter = T.Distance_to_point;
Timestamp = db.Matched_data.Matched_Neural_data.Timestamp;

Patch_x = db.Matched_data.Matched_Neural_data.Timestamp;  
Patch_y = Dist2Shelter;
Shelter_states = zeros(numel(Timestamp),1);
Shelter_states(find(Dist2Shelter<=8)) = 1;
minimum_frame_num = 15;
[OnOff_Frames, Adjusted_States] = Convert_States2OnOffFrames(Shelter_states ,minimum_frame_num);

Event_on_off = OnOff_Frames;
roi_x = [Patch_x(1), Patch_x(end)];
figure
set(gcf,'color','w','Position',[100, 100, 800, 400])
% Plot Shelter Distance
ax1 = subplot(9,1,1:2); % Patch Shelter
plot(Patch_x, Patch_y,'k'); hold on;
yl = ylim;
% Draw Stim time
for(i=1:numel(db.Matched_data.Matched_Behavior_data.All_stim_time))
    line([db.Matched_data.Matched_Behavior_data.All_stim_time(i,1) db.Matched_data.Matched_Behavior_data.All_stim_time(i,1)],ylim,'color','r'); hold on;
end        
set(gca,'tickdir','out','box','on','YTick',[0 50],'xtick',[])
hold off;
drawnow
xlim(roi_x)

% Plot Neural Activity
ax2 = subplot(9,1,3:9);
ncells = size(DFoF,2);
imagesc(Patch_x,1:ncells,zscore_DFoF(:,[Shelter_Neurons(Shelter_I); Nonspecific_Neurons(Nonspecific_I); Outside_Neurons(Outside_I)]')'); hold on;
caxis([-2 2]);
colormap(viridis);
set(gca,'YDir','normal','YTick',[1 numel(Shelter_Neurons) numel(Shelter_Neurons)+numel(Nonspecific_Neurons) ncells],'box','off', 'fontsize',12,'FontWeight','bold','tickdir','out');
xlabel('Time (s)');    
ylabel('Neuron');
xlim(roi_x)        
linkaxes([ax1, ax2],'x')


%% Shelter In Demo Figure
Shelter_states = zeros(size(T,1),1);
Shelter_states(find(Dist2Shelter<5))=1;
minimum_frame_num = 5;
[OnOff_Frames, Adjusted_States] = Convert_States2OnOffFrames(Shelter_states,minimum_frame_num);

Shelter_in_out_time = Timestamp(OnOff_Frames);
ROI_row = [21,23,28,32,42];  % Define ROIs for specific neurons
for(i=1:numel(ROI_row))    
    Shelter_in_time = Shelter_in_out_time(ROI_row(i),1);
    xrange = Shelter_in_time + [-2 10];
    
    figure
    set(gcf,'color','w','Position',[300, 300, 150, 600])
    % Plot Shelter Distance
    ax1 = subplot(10,1,1:3) % Patch Shelter
    plot(Patch_x, Patch_y,'k'); hold on;
    xl = xlim, 
    line([xl(1) xl(2)], [5 5]); hold on;
    ylim([0 50]);    
    set(gca,'tickdir','out','box','off','YTick',[0 50],'xtick',[])
    hold off;
    drawnow
    xlim(xrange)
    
    % Plot Neural Activity
    ax2 = subplot(10,1,4:7);
    ncells = size(DFoF,2);
    imagesc(Patch_x,1:ncells,zscore_DFoF(:,[Shelter_Neurons(Shelter_I); Nonspecific_Neurons(Nonspecific_I); Outside_Neurons(Outside_I)]')'); hold on;
    caxis([-2 2]);
    colormap(viridis);
    set(gca,'YDir','normal','xTick',[], 'YTick',[1 numel(Shelter_Neurons) numel(Shelter_Neurons)+numel(Nonspecific_Neurons) ncells],'box','off', 'fontsize',12,'FontWeight','bold','tickdir','out');
    xlim(xrange)        
    
    % Population activity
    ROI_Neurons_g1 = Outside_Neurons;
    ROI_Neurons_g2 = Nonspecific_Neurons;
    ROI_Neurons_g3 = Shelter_Neurons;
    x = Timestamp;  
    y1 = nanmean(zscore_DFoF(:,ROI_Neurons_g1),2);
    y2 = nanmean(zscore_DFoF(:,ROI_Neurons_g2),2);
    y3 = nanmean(zscore_DFoF(:,ROI_Neurons_g3),2);
    
    ax3 = subplot(10,1,8:10);
    g1_color = [210 52 131]/255;
    g2_color = [128 128 128]/255;
    g3_color = [59 197 124]/255;
    plot(x, y1,'color',g1_color); hold on;
    plot(x, y2,'color',g2_color); hold on;
    plot(x, y3,'color',g3_color); hold on;
    
    set(gca,'tickdir','out','box','off','xtick',[ceil(xrange(1)):5:ceil(xrange(2))])
    xlim(xrange)        
    linkaxes([ax1, ax2, ax3],'x')
end