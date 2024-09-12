%% RNAscope data classification analysis
% Created by Kanghoon Jung
% Contact: kanghoon.jung AT alleninstitute.org for bug reports, questions, or feature requests.

close all 
clear all

data_folder = 'C:\Users\jungk\Google Drive\Github\Jung_NatNeuro_2024\codes\Fig7\';
cd(data_folder)
% Load the dataset from the .mat file
load("Fig7i_data.mat");


%% Prediction using Classification_CART_conc_data_table function
Conc_Data_table = Classification_CART_conc_data_table(Conc_Data_table);

%% Cross-validation analysis to calculate misclassification rate
Manual_set_rows = (Conc_Data_table.Role=="Manual");
Predict_set_rows = (Conc_Data_table.Role=="Predict");
    
Pred_CV_Data_table = Conc_Data_table(Predict_set_rows,:);
input_Names = {'Chat_F_rescaled'; 'DAPI_F_rescaled'; 'Tdt_F_rescaled'; 'D1_F_rescaled'; 'D2_F_rescaled'};
target_Names = [];
for(i=1:5)
    target_Names{i,1} = [input_Names{i} '_Prediction_DA'];
end

% Perform cross-validation and calculate the misclassification rate
Pred_Misclassification_rate = CrossValidation(Pred_CV_Data_table, input_Names, target_Names);
Pred_Misclassification_rate*100

%% Classification Analysis
Classification_results = [];
Classification_array = [];
input_Names = {'Chat_F_rescaled'; 'DAPI_F_rescaled'; 'Tdt_F_rescaled'; 'D1_F_rescaled'; 'D2_F_rescaled'};
% Loop through each input feature to perform classification
for(i=1:5)
    target_Names{i,1} = [input_Names{i} '_Prediction_DA'];
    Classification_results.(input_Names{i}).positive = find(Conc_Data_table.(target_Names{i,1})==1);
    Classification_results.(input_Names{i}).negative = find(Conc_Data_table.(target_Names{i,1})==0);    
    Classification_results.(input_Names{i}).positive_num = numel(find(Conc_Data_table.(target_Names{i,1})==1));
    Classification_results.(input_Names{i}).negative_num = numel(find(Conc_Data_table.(target_Names{i,1})==0));    
    Classification_results.(input_Names{i}).percent = Classification_results.(input_Names{i}).positive_num/size(Conc_Data_table,1)*100;    
    Classification_array(:,i) = [Classification_results.(input_Names{i}).positive_num; Classification_results.(input_Names{i}).negative_num; Classification_results.(input_Names{i}).percent];
end

% Extract specific cell types based on classification results
Chat_cells = Classification_results.Chat_F_rescaled.positive;
Tdt_cells = Classification_results.Tdt_F_rescaled.positive;
NAc_cells = Classification_results.DAPI_F_rescaled.positive;
D1_cells = Classification_results.D1_F_rescaled.positive;
D2_cells = Classification_results.D2_F_rescaled.positive;


%% Combine and identify intersections between cell types
Combined = unique([Chat_cells; D1_cells;D2_cells]);
Other_cells = setdiff(NAc_cells,Combined);

% Find intersections between different cell types
Intersect_D1D2_cells = intersect(D2_cells,D1_cells)
Intersect_D1Chat_cells = intersect(D1_cells,Chat_cells)
Intersect_D2Chat_cells = intersect(D2_cells,Chat_cells)
Intersect_D1D2Chat_cells = intersect(Intersect_D1D2_cells,Chat_cells)

% Calculate percentages of intersections

Intersect_D1D2_cells_perc = numel(Intersect_D1D2_cells)/numel(union(D1_cells,D2_cells))*100
Intersect_D1Chat_cells_perc = numel(Intersect_D1Chat_cells)/numel(union(D1_cells,Chat_cells))*100
Intersect_D2Chat_cells_perc = numel(Intersect_D2Chat_cells)/numel(union(D2_cells,Chat_cells))*100
Intersect_D1D2Chat_cells_perc = numel(Intersect_D1D2Chat_cells)/numel(Combined)*100

%%
x = 1; 
y = [3688 4735 180 20 6 12 0]./sum([3688 4735 180 20 6 12 0])*100; 
b = bar(x,y,"stacked");



% Percentage of active neurons in NAc neurons
Active_NAc_neurons.num =  numel(find(intersect(Tdt_cells, NAc_cells)));
Active_NAc_neurons.perc =  numel(find(intersect(Tdt_cells, NAc_cells)))/numel(Tdt_cells)*100;

% Percentage of active neurons in Chat neurons
Active_Chat_neurons.num = numel(find(intersect(Tdt_cells, Chat_cells)));
Active_Chat_neurons.perc =  numel(find(intersect(Tdt_cells, Chat_cells)))/numel(Tdt_cells)*100;
Chat_Active_neurons.num =  numel(find(intersect(Tdt_cells, Chat_cells)));
Chat_Active_neurons.perc =  numel(find(intersect(Tdt_cells, Chat_cells)))/numel(Chat_cells)*100;

% Percentage of active neurons in D1 neurons
Active_D1_neurons.num =  numel(find(intersect(Tdt_cells, D1_cells)));
Active_D1_neurons.perc =  numel(find(intersect(Tdt_cells, D1_cells)))/numel(Tdt_cells)*100;
D1_Active_neurons.num =  numel(find(intersect(Tdt_cells, D1_cells)));
D1_Active_neurons.perc =  numel(find(intersect(Tdt_cells, D1_cells)))/numel(D1_cells)*100;

% Percentage of active neurons in D2 neurons
Active_D2_neurons.num =  numel(find(intersect(Tdt_cells, D2_cells)));
Active_D2_neurons.perc =  numel(find(intersect(Tdt_cells, D2_cells)))/numel(Tdt_cells)*100;
D2_Active_neurons.num =  numel(find(intersect(Tdt_cells, D2_cells)));
D2_Active_neurons.perc =  numel(find(intersect(Tdt_cells, D2_cells)))/numel(D2_cells)*100;

% Percentage of active neurons in other neurons
Active_Other_neurons.num =  numel(find(intersect(Tdt_cells, Other_cells)));
Active_Other_neurons.perc =  numel(find(intersect(Tdt_cells, Other_cells)))/numel(Tdt_cells)*100;
Other_Active_neurons.num =  numel(find(intersect(Tdt_cells, Other_cells)));
Other_Active_neurons.perc =  numel(find(intersect(Tdt_cells, Other_cells)))/numel(Other_cells)*100;

% Summary of active neurons
Active_neurons.num = [Active_Chat_neurons.num Active_D1_neurons.num Active_D2_neurons.num Active_NAc_neurons.num-sum([Active_Chat_neurons.num Active_D1_neurons.num Active_D2_neurons.num])]';
Active_neuron_in_each_celltype = [Chat_Active_neurons.num, numel(Chat_cells); D1_Active_neurons.num, numel(D1_cells); D2_Active_neurons.num, numel(D2_cells); Other_Active_neurons.num, numel(Other_cells)];

%% Scatter plot for data visualization
T = [Conc_Data_table.D1_F_rescaled, Conc_Data_table.D2_F_rescaled, Conc_Data_table.Chat_F_rescaled, Conc_Data_table.Tdt_F_rescaled]; 
T_D1 = T(D1_cells,:);
T_D2 = T(D2_cells,:);
T_Chat = T(Chat_cells,:);
T_Other = T(Other_cells,:);

% Assign index to different cell types for scatter plot
T_Index = [];
T_Index(D1_cells,1) = 1;
T_Index(D2_cells,1) = 2;
T_Index(Chat_cells,1) = 3;
T_Index(Other_cells,1) = 4;

T2 = [T_Index, Conc_Data_table.Tdt_F_rescaled];
T2(find(T2(:,1)==0),:) = [];
% x2 = Conc_Data_table.D2_F_rescaled; 
% x3 = Conc_Data_table.Chat_F_rescaled; 
dot_size = 20;
figure, 
set(gcf,'color','w')
h = scatter3(T_D1(:,1),T_D1(:,2),T_D1(:,3),dot_size, [230 220 50]/255,'filled'); hold on;
h = scatter3(T_D2(:,1),T_D2(:,2),T_D2(:,3),dot_size, 'm','filled'); hold on;
h = scatter3(T_Chat(:,1),T_Chat(:,2),T_Chat(:,3),dot_size, [229.5 19.89 46.92]/255,'filled'); hold on;
h = scatter3(T_Other(:,1),T_Other(:,2),T_Other(:,3),dot_size, [0.5 0.5 0.5],'filled'); hold on;

alpha = 0.8;
set(h, 'MarkerEdgeColor','none', 'MarkerFaceAlpha', alpha)
xlabel('Drd1');
ylabel('Drd2');
zlabel('Chat');
view(-60,60)
 

%% Statistical analysis
% One-way ANOVA
[p, tbl, stats] = anova1(T2(:,2), T2(:,1))
[c,m,h,nms] = multcompare(stats,'alpha',.05,'ctype','bonferroni');

%% Kruskal-Wallis H test (non-parametric ANOVA)
T2 = [T_Index, Conc_Data_table.Tdt_F_rescaled];
g1_DATA = T2(T2(:,1)==1,2);
g2_DATA = T2(T2(:,1)==2,2);
g3_DATA = T2(T2(:,1)==3,2);
g4_DATA = T2(T2(:,1)==4,2);
% Combine data for all groups
g_DATA = [g1_DATA; g2_DATA; g3_DATA; g4_DATA];
group = [ones(numel(g1_DATA),1);2*ones(numel(g2_DATA),1); 3*ones(numel(g3_DATA),1); 4*ones(numel(g4_DATA),1)];

% Perform Kruskal-Wallis test and perform multiple comparisons using Bonferroni correction
[p, tbl, stats] = kruskalwallis(g_DATA, group, 'off');
[c, m, h, nms] = multcompare(stats, 'ctype', 'bonferroni');
