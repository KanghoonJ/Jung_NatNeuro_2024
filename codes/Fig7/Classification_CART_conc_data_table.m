%% CART model
% % Classification Tree
% load ionosphere
% tc = fitctree(X,Y)
% clear all
% load carsmall
% X = [Weight grp2idx(Origin)];
% a = ClassificationTree.fit(X,Cylinders,'cat',2);
% view(a,'mode','graph')
% 
% % Regression Tree
% tree = fitrtree([Weight, Cylinders],MPG,'CategoricalPredictors',2,'MinParentSize',20,'PredictorNames',{'W','C'})
% MPG4Kpred = predict(tree,[4000 4; 4000 6; 4000 8])
function Data_table = Classification_CART_conc_data_table(Data_table)

%     Data_table = Conc_Data_table;
    %% RNAscope data
    % load('ML_Data_table.mat')
    VariableNames = Data_table.Properties.VariableNames;
    Cell_types = {'Chat', 'DAPI', 'Tdt', 'D1', 'D2'};
    % Input_var_names = {'Chat_F', 'DAPI_F', 'Tdt_F', 'D1_F', 'D2_F'};
    Input_var_names = {'Chat_F_rescaled', 'DAPI_F_rescaled', 'Tdt_F_rescaled', 'D1_F_rescaled', 'D2_F_rescaled'};
    for(i=1:5)
        Input_var_columns(1,i) = find(strcmp(VariableNames, Input_var_names{i}));
    end
    Manual_rows = (Data_table.Role=="Manual");
    Manual_sample_num = numel(find(Manual_rows==true));
    Manual_set_Input = [Data_table(Manual_rows,Input_var_columns)];

    Test_rows = (Data_table.Role=="Predict");
    Test_set_Input = [Data_table(Test_rows,Input_var_columns)];
    for(i=1:5)
        Input_Name = Input_var_names{i};
        % Manual classification results
        Manual_set_Output = Data_table(Manual_rows,:).(Cell_types{i});
        Target = Manual_set_Output;

        % 2) Predict labels using discriminant analysis classification model
        Mdl = fitcdiscr(Manual_set_Input, Manual_set_Output);
        Test_set_Output = predict(Mdl,Test_set_Input);
        Label = Test_set_Output;
        Test_set_Output_Name = [Input_Name '_Prediction_DA'];
        Data_table.(Test_set_Output_Name)(Test_rows,1) = Test_set_Output;
    %     Performance(1,i) = sum(Label(Manual_rows,:)==Target)/numel(Label(Manual_rows,:)); % performance in the range of 0 to 1;
end



% Performance = [];
% for(i=1:5)
%     Input_Name = VariableNames{find(strcmp(VariableNames,'Chat'))-1+i};
%     Manual_set_Output = Data_table(Manual_rows,:).(Input_Name);
%     Target = Manual_set_Output;
% % Manual_set_a = ClassificationTree.fit(Manual_set_Input,Manual_set_Output,'cat',2);
% % view(Manual_set_a,'mode','graph')
% 
%     % 1) Regression Tree
%     % tree = fitrtree(Manual_set_Input, Manual_set_Output,'CategoricalPredictors',2,'MinParentSize',20,'PredictorNames',{'D1_F','D2_F'});
%     tree = fitrtree(Manual_set_Input, Manual_set_Output,'MinParentSize',floor(Manual_sample_num/3),'PredictorNames',VariableNames(2:6));
%     Test_set_Input = [Data_table(:,2:6)];
%     Test_set_Output = predict(tree, Test_set_Input);
%     Label = Test_set_Output;
%     Test_set_Output_Name = [Input_Name '_Prediction_RT'];
%     Data_table.(Test_set_Output_Name) = Test_set_Output;
%     Performance(1,i) = sum(Label(Manual_rows,:)==Target)/numel(Label(Manual_rows,:)); % performance in the range of 0 to 1;
%     
% end

% 
% 
% % Example
% %% Manual Side
% % random three class data with target matrix -- [9X3] 9 observation with 3 features
% data = [10 0 0;
%         10 0 1;
%         10 1 0;
%         2 10 0;
%         2 10 1;
%         2 11 0;
%         5 0 10;
%         5 0 11;
%         5 1 10];
% target = [1;1;1;2;2;2;3;3;3];
% % create a decision tree
% dtree = fitctree(data,target,'MinParentSize',2,'CrossVal','on'); % dtree is the trained model. save it at end for doing testing
% view(dtree.Trained{1},'Mode','graph')
% 
% % save('dtree.mat','dtree');
% % train performance
% label = predict(dtree,data);
% perf=sum(label==target)/size(label,1) % performance in the range of 0 to 1
% 
% 
% 
% 
% %% Testing Side
% % for testing load the trained model
% load('dtree.mat');
% testdata = [5 3 10; 10 3 1]; % take 1 new unknown observation and give to trained model
% Group = predict(dtree,testdata);
% 
% %% Empirical data
% data = Manual_set(:,3:6);
% target = Manual_set(:,2);
% dtree = fitctree(data,target,'MinParentSize',2); % dtree is the trained model. save it at end for doing testing
% save('dtree.mat','dtree');
% % train performance
% label = predict(dtree,data);
% perf=sum(label==target)/size(label,1) % performance in the range of 0 to 1
% testdata = F(:,1:100)';
% Group = predict(dtree,testdata);

