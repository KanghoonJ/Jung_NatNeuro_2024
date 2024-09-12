%% Fiber photometry data analysis code
%% Code for EDF2e
% Created by Kanghoon Jung
% Contact: kanghoon.jung AT alleninstitute.org for bug reports, questions, or feature requests.

clear all 
close all

data_folder = 'C:\Users\jungk\Google Drive\Github\Jung_NatNeuro_2024\codes\EDF2\';
cd(data_folder)

% Load the dataset from the .mat file
load("EDF2e_data.mat");
Num_Frames = [];
Num_Threats = [];
for(nSubj=1:size(T_Indv_ROI_DFoF,1))    
    Num_Frames(nSubj,1) = size(T_Indv_ROI_DFoF{nSubj,1},1);
    Num_Threats(nSubj,1) = size(T_Indv_ROI_DFoF{nSubj,1},2);
end
l_col = parula(max(Num_Threats));
nSubj=1;
Norm_DFoF = normalize(T_Indv_ROI_DFoF{nSubj,1},2,'range');

% Plot
figure, 
set(gcf,'color','w','position',[200 200 300 500])
subplot(311)
y = 1:size(Norm_DFoF,2);
z = T_Indv_ROI_Dist2Shelter{nSubj,1}';
imagesc(x,y,z);
xlim(x_window_range)        
xlabel('Time from Event onset')
title('Dist2Shelter')
caxis([0 50]);
sgtitle(Case)
subplot(312)
z = T_Indv_ROI_Speed{nSubj,1}';
imagesc(x,y,z);
caxis([0 30]);
xlim(x_window_range)        
title('Speed')
subplot(313)
z = T_Indv_ROI_DFoF{nSubj,1}';
imagesc(x,y,z);
xlim(x_window_range)   
title('Normalized DFoF');
saveas(gcf,'Threat_DFoF_by_trials.pdf')
