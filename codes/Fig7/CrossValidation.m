function Misclassification_rate = CrossValidation(CV_Data_table, input_Names, target_Names)
% Cross validation
% CV_Data_table = Data_table(ROI_Cells,:);

Misclassification_rate = [];
X = table2array(CV_Data_table(:,2:6));
for(i=1:5)    
    y = num2str(CV_Data_table.(target_Names{i}));
    % Create a random partition for a stratified 10-fold cross validation
    c = cvpartition(y,'KFold',10);
    fun = @(xTrain,yTrain,xTest,yTest)(sum(~strcmp(yTest,classify(xTest,xTrain,yTrain))));    
    Misclassification_rate(1,i) = sum(crossval(fun,X,y,'partition',c))/sum(c.TestSize);
end
% 
% %% Test
% CV_Data_table = Data_table(ROI_Cells,:);
% input_Names = {'Chat'; 'DAPI'; 'Tdt'; 'D1'; 'D2'};
% target_Names = [];
% for(i=1:5)
%     target_Names{i,1} = [input_Names{i} '_Prediction'];
% end
% 
% Misclassification_rate = [];
% X = table2array(CV_Data_table(:,2:6));
% for(i=1:5)    
%     y = num2str(CV_Data_table.(target_Names{i}));
%     % Create a random partition for a stratified 10-fold cross validation
%     c = cvpartition(y,'KFold',10);
%     fun = @(xTrain,yTrain,xTest,yTest)(sum(~strcmp(yTest,classify(xTest,xTrain,yTrain))));    
%     Misclassification_rate(1,i) = sum(crossval(fun,X,y,'partition',c))/sum(c.TestSize);       
% end
% Misclassification_rate

% %% Plot test
% S_x1 = table2array(CV_Data_table(:,4));
% S_x2 = table2array(CV_Data_table(:,5));
% group = table2array(CV_Data_table(:,11));
% h1 = gscatter(S_x1,S_x2,group,'rb','v^',[],'off');
% set(h1,'LineWidth',2)
% legend('D1','D2','Location','NW')
% 
% classify(xTest,xTrain,yTrain)

