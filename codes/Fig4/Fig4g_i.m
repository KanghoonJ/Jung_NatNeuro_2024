clear all
close all

col_outside = [210 52 131]/255; 
col_shelter = [59 197 124]/255; 
col_nonspecific = [128 128 128]/255;

% 'VariableNames',["Subject","Shelter_Neurons","Outside_Neuron","Nonspecific_Neurons","Shelter_Neurons_Perc","Outside_Neurons_Perc","Nonspecific_Neurons_Perc"]);
for(nG=3:3)
    Pop_table = [];
    if(nG==1)
        Group = 'NAc';
        Base_folder = 'G:\Project 5 - GP_Raw Data\Project #5-GP Rev\Experiments\Miniscope\NAc';
    elseif(nG==2)
        Group = '6-OHDA';
        Base_folder = 'G:\Project 5 - GP_Raw Data\Project #5-GP Rev\Experiments\Miniscope\6-OHDA@NAc';
    elseif(nG==3)
        Group = 'vHPC-CNO';
        Base_folder = 'G:\Project 5 - GP_Raw Data\Project #5-GP Rev\Experiments\Miniscope\hG6s_rghM4Di@NAc CaMKII-Cre@vHPC\CNO';
    elseif(nG==4)
        Group = 'vHPC-DMSO';
        Base_folder = 'G:\Project 5 - GP_Raw Data\Project #5-GP Rev\Experiments\Miniscope\hG6s_rghM4Di@NAc CaMKII-Cre@vHPC\DMSO';
    end

    cd(Base_folder)
    A = dir;
    A = A(~ismember({A.name},{'.','..'}));
    Folders = A([A.isdir]);
    Num_Subj = size(Folders,1);
    
    % Conc_Shelter_preference_Table = [];
    % Conc_Shelter_preference_Table.Time = [300:300:3600];
    % C_Shelter_preference_ratio = [];
    % k = 1;
    Conc_T = [];
    % for nSubj = 5:5
    for nSubj = 1:Num_Subj        
        T = table();
        cd(Base_folder)
        cd(Folders(nSubj).name)        
        disp('Loading matched db')
        db_file = dir('*matched.mat');
        
        Date = {(db_file(:).date)}';
        [Sorted_Date, I] = Sort_date(Date);
        % Most recent one
        db_file(I(end)).name;        
        % Load most recent analysis datafile
        load(db_file(I(end)).name);
        display(db.DAQ.Subj_id)
    
        Shelter_Neurons = db.Analysis.Miniscope.True_Shelter_Neurons;
        Outside_Neurons = db.Analysis.Miniscope.True_Outside_Neurons;
        Nonspecific_Neurons = db.Analysis.Miniscope.Nonspecific_Neurons;
        NCells = size(db.Analysis.Miniscope.ROI_cell_F,2);
        
        Pop_table.Subject{nSubj,1} = db.DAQ.Subj_id;
        Pop_table.Shelter_Neurons{nSubj,1} = numel(Shelter_Neurons);
        Pop_table.Outside_Neurons{nSubj,1} = numel(Outside_Neurons);
        Pop_table.Nonspecific_Neurons{nSubj,1} = numel(Nonspecific_Neurons);
        Pop_table.Shelter_Neurons_Perc{nSubj,1} = numel(Shelter_Neurons)/(numel(Shelter_Neurons)+numel(Outside_Neurons)+numel(Nonspecific_Neurons))*100;
        Pop_table.Outside_Neurons_Perc{nSubj,1} = numel(Outside_Neurons)/(numel(Shelter_Neurons)+numel(Outside_Neurons)+numel(Nonspecific_Neurons))*100;
        Pop_table.Nonspecific_Neurons_Perc{nSubj,1} = numel(Nonspecific_Neurons)/(numel(Shelter_Neurons)+numel(Outside_Neurons)+numel(Nonspecific_Neurons))*100;



        T.Subject = nSubj*ones(NCells,1);
        T.Cell = [1:NCells]';
        Response_type = nan(NCells,1);
        T.Response_type(Shelter_Neurons,1) = 1;
        T.Response_type(Outside_Neurons,1) = -1;
        T.Response_type(Nonspecific_Neurons,1) = 0;
        
        Shelter_preference_Table = db.Analysis.Miniscope.Shelter_preference_Table;       
        Shelter_preference_ratio = nan(NCells,13);
        Frames = 1:size(Shelter_preference_Table,1);
        Shelter_preference_ratio(:,Frames) = Shelter_preference_Table.Shelter_preference_ratio';
        Shelter_preference_ratio(~isfinite(Shelter_preference_ratio)) = nan;
        % Shelter_preference_ratio = fillmissing(Shelter_preference_ratio,'nearest',2);

        T.Shelter_preference_ratio = Shelter_preference_ratio;
        Conc_T = [Conc_T; T];
    %     Conc_Shelter_preference_Table.Time(Frames) = Shelter_preference_Table.Time(Frames);
    % %     Conc_Shelter_preference_Table.Shelter_preference_ratio(Frames,:);
    %     for(nC=1:size(Shelter_preference_Table.Shelter_preference_ratio,2))
    %         
    %         C_Shelter_preference_ratio(Frames,k) = Shelter_preference_Table.Shelter_preference_ratio(Frames,nC);
    %         k = k+1;
    %     end
    end
    Conc_T(find(all(isnan(Conc_T.Shelter_preference_ratio),2)),:) = []; 
    init_y = [];
    for(i=1:size(Conc_T,1))
        % init_y(i,1) = Conc_T.Shelter_preference_ratio(i,find(~isnan(Conc_T.Shelter_preference_ratio(i,:)),1));
        init_y(i,1) = Conc_T.Shelter_preference_ratio(i,1);        
    end    
    init_y(find(isnan(init_y))) = 0;
    Conc_T.Shelter_preference_ratio_change = Conc_T.Shelter_preference_ratio - init_y;

    % Shelter_dy = dy(find(Conc_T.Response_type==1),10)
    % Shelter_dy(~isfinite(Shelter_dy)) = [];
    
    cd(Base_folder)       
    save('Shelter_preference_ratio_over_time.mat', "Conc_T")
    % Conc_Shelter_preference_Table.C_Shelter_preference_ratio = C_Shelter_preference_ratio;    


    fig_1 = figure;
    set(gcf,'color','w','Position', [200, 200, 400, 400]);
    x = [300:300:3900]/60;
    y = Conc_T.Shelter_preference_ratio;
    y(~isfinite(y)) = nan;
    y = y - y(:,1);
    for(t=-1:1:1)
        ROI_y = y(find(Conc_T.Response_type==t),:);
        if(t==-1)
            col = [210 52 131]/255; %col_outside
        elseif(t==0)
            col = [128 128 128]/255; %col_nonspecific
        elseif(t==1)
            col = [59 197 124]/255; %col_shelter
        end    
    %     subplot(131)
    %     for(i=1:size(ROI_y,1))
    %         plot(x, ROI_y(i,:),'color',col); hold on;
    %     end
    %     xlim([0,55])
    %     set(gca,'box','off','xtick',[0:10:50],'tickdir','out')    
    %     subplot(132)    
        shadedErrorBar(x(1:10), nanmean(ROI_y(:,1:10)),nansem(ROI_y(:,1:10)),'lineProps',{'color',col,'markerfacecolor',col}); hold on;
    %     shadedErrorBar(x,  'color',col); hold on;
        ylim([-0.7 1.2])
        xlim([0,55])
        set(gca,'box','off','xtick',[0:10:50],'tickdir','out')
    end
    Fig_data.Shelter_preference_ratio_over_time.x = x(1:10);
    Fig_data.Shelter_preference_ratio_over_time.y = nanmean(ROI_y(:,1:10));
    Fig_data.Shelter_preference_ratio_over_time.y_sem = nansem(ROI_y(:,1:10));
    xlabel('Time (min)')
    ylabel('Shelter preference index')
    cd(Base_folder)
    saveas(fig_1,'Shelter_preference_ratio_over_time.pdf')
    
    %% Stacked chart
    col_outside = [210 52 131]/255; 
    col_shelter = [59 197 124]/255; 
    col_nonspecific = [128 128 128]/255;
    
    Shelter_Neurons_Population_T = [numel(find(Conc_T.Response_type==-1)) numel(find(Conc_T.Response_type==0)) numel(find(Conc_T.Response_type==1))]
    % figure
    % set(gcf,'color','w')
    % subplot(121)
    % b1 = bar(1,table2array(Shelter_Neurons_Population_T(1,:)),'stacked','Facecolor','flat');
    % b1(1).BarWidth = 0.25;
    % b1(1).CData = col_outside;
    % b1(2).CData = col_nonspecific;
    % b1(3).CData = col_shelter;
    % set(gca,'ydir','reverse','box','off')
    % subplot(122)
    
    fig_2 = figure
    b2 = bar(2,Shelter_Neurons_Population_T/sum(Shelter_Neurons_Population_T)*100,'stacked','Facecolor','flat');
    b2(1).BarWidth = 0.25;
    b2(1).CData = col_outside;
    b2(2).CData = col_nonspecific;
    b2(3).CData = col_shelter;
    set(gca,'ydir','reverse','box','off')
    sgtitle(Group)
    
    cd(Base_folder)
    saveas(fig_2,'Shelter_preference_population_bar.pdf')
    
    fig_3 = figure
    set(gcf,'color','w','Position', [500, 200, 450, 450]);
    ax = gca();
    h = pie(ax,[Shelter_Neurons_Population_T]);
    ax.Colormap = [col_outside;col_nonspecific;col_shelter]; 
    % Create legend
    labels = {'Outside','Non-specific','Shelter'};
    lgd = legend(labels,'Location','southoutside')
    sgtitle(Group)
    cd(Base_folder)
    saveas(fig_3,'Shelter_preference_population_chart.pdf')
    
    
    fig_4 = figure; 
    set(gcf,'color','w','Position', [200, 200, 450, 450]);
    x = [300:300:3900];
    y(~isfinite(y)) = nan;
    
    Data_cdf = [];
    % col_map = zeros(size(y,2));
    col_map(1:10,:) = viridis(10);
    % for(t=1:size(y,2))
    for(t=1:10)
        if(t<=10)
            ROI_y = y(:,t);
            Data_cdf{t} = gen_cdf(ROI_y);
            plot(Data_cdf{t}.x, Data_cdf{t}.y,'color',col_map(t,:)); hold on;
        end
    end
    colormap(col_map)
    colorbar('location','southoutside');
    xlim([-2.2 3])
    sgtitle(Group)
    set(gca,'box','off','tickdir','out','ytick',[0:0.2:1])
    
    cd(Base_folder)
    saveas(fig_4,'Shelter_preference_ratio_over_time_cdf.pdf')
    
    fig_5 = figure;
    set(gcf,'color','w','Position', [200, 200, 400, 400]);
    x = [300:300:3900]/60;
    y = Conc_T.Shelter_preference_ratio;
    y(~isfinite(y)) = nan;
    % y = y - y(:,1);
    y = y;
    Conc_T.Shelter_preference_ratio_change = y;
    for(t=-1:1:1)
        ROI_y = y(find(Conc_T.Response_type==t),:);
        if(t==-1)
            col = [210 52 131]/255; %col_outside
        elseif(t==0)
            col = [128 128 128]/255; %col_nonspecific
        elseif(t==1)
            col = [59 197 124]/255; %col_shelter
        end    
        subplot(121)
        for(i=1:size(ROI_y,1))
            plot(x, ROI_y(i,:),'color',col); hold on;
        end
        xlim([0,55])
        set(gca,'box','off','xtick',[0:10:50],'tickdir','out')    
        subplot(122)    
        shadedErrorBar(x(1:10), nanmean(ROI_y(:,1:10)),nansem(ROI_y(:,1:10)),'lineProps',{'color',col,'markerfacecolor',col}); hold on;    
        ylim([-2, 2])
        xlim([0,55])
        set(gca,'box','off','xtick',[0:10:50],'tickdir','out')
    end
    sgtitle(Group)
    Fig_data.Shelter_preference_ratio_over_time.x = x(1:10);
    Fig_data.Shelter_preference_ratio_over_time.y = nanmean(ROI_y(:,1:10));
    Fig_data.Shelter_preference_ratio_over_time.y_sem = nansem(ROI_y(:,1:10));
    xlabel('Time (min)')
    ylabel('Shelter preference index')
    cd(Base_folder)
    saveas(fig_5,'Shelter_preference_ratio_over_time_indiv.pdf')


    Fig_data.Shelter_preference_ratio_cdf.Data_cdf = Data_cdf;
    Fig_data.Shelter_preference_ratio_cdf.Data_cdf = Data_cdf;
    
    save('Fig_data.mat','Fig_data')
end
Pop_table = struct2table(Pop_table)
% 
% 
% 
%% Average
% close all
for(nG=4:4)
    if(nG==1)
        Group = 'NAc';
        Base_folder = 'G:\Project 5 - GP_Raw Data\Project #5-GP Rev\Experiments\Miniscope\NAc';
    elseif(nG==2)
        Group = '6-OHDA';
        Base_folder = 'G:\Project 5 - GP_Raw Data\Project #5-GP Rev\Experiments\Miniscope\6-OHDA@NAc';
    elseif(nG==3)
        Group = 'vHPC-CNO';
        Base_folder = 'G:\Project 5 - GP_Raw Data\Project #5-GP Rev\Experiments\Miniscope\hG6s_rghM4Di@NAc CaMKII-Cre@vHPC\CNO';
    elseif(nG==4)
        Group = 'vHPC-DMSO';
        Base_folder = 'G:\Project 5 - GP_Raw Data\Project #5-GP Rev\Experiments\Miniscope\hG6s_rghM4Di@NAc CaMKII-Cre@vHPC\DMSO';
    end

    cd(Base_folder)
    A = dir;
    A = A(~ismember({A.name},{'.','..'}));
    Folders = A([A.isdir]);

    load('Shelter_preference_ratio_over_time.mat')
    NSubj = max(unique(Conc_T.Subject));
    Avg_Conc_T = table();
    figure;
    x = [300:300:3900]/60;    
    k = 1;
    for(nSubj=1:NSubj)        
        for(t=-1:1)
            if(t==-1)
                col = [210 52 131]/255; %col_outside
            elseif(t==0)
                col = [128 128 128]/255; %col_nonspecific
            elseif(t==1)
                col = [59 197 124]/255; %col_shelter
            end               
            Avg_Conc_T.Subject{k,1} = Folders(nSubj).name;        
            Avg_Conc_T.Response_type(k,1) = t;            
            ROI_y = Conc_T(find(Conc_T.Subject==nSubj & Conc_T.Response_type==t),:);
            ROI_y.Shelter_preference_ratio(~isfinite(ROI_y.Shelter_preference_ratio)) = nan;
            Avg_Conc_T.Shelter_preference_ratio_x(k,1:numel(x)) = x;
            

            % y = ROI_y.Shelter_preference_ratio;
            % y(~isfinite(y)) = nan;
            % init_y = [];
            % for(i=1:size(ROI_y,1))
            %     init_y(i,1) = ROI_y.Shelter_preference_ratio(i,find(~isnan(ROI_y.Shelter_preference_ratio(i,:)),1));
            % end
            % dy = y - init_y;
            % % dy = y ;
            % Avg_Conc_T.Shelter_preference_ratio(k,1:numel(mean(dy,1,'omitnan'))) = mean(dy,1,'omitnan');
            % Avg_Conc_T.Shelter_preference_ratio_sem(k,1:numel(nansem(dy)))= nansem(dy);
            % k = k+1;      
            % shadedErrorBar(x(1:10), mean(dy(:,1:10),1,'omitnan'),nansem(dy(:,1:10),0,1),'lineProps',{'color',col,'markerfacecolor',col}); hold on;
            % % plot(x(1:10), mean(dy(:,1:10),1,'omitnan'),'color',col); hold on;

            y = ROI_y.Shelter_preference_ratio_change;
            k = k+1;      
            shadedErrorBar(x(1:10), mean(y(:,1:10),1,'omitnan'),nansem(y(:,1:10),0,1),'lineProps',{'color',col,'markerfacecolor',col}); hold on;
            % % plot(x(1:10), mean(dy(:,1:10),1,'omitnan'),'color',col); hold on;
        end
    end
    sgtitle(Group)
    title('SPI change')

    % figure;
    % x = [300:300:3900]/60;    
    % k = 1;
    % Avg_Conc_T = table();
    % for(t=-1:1)
    %     for(nSubj=1:NSubj)        
    %         if(t==-1)
    %             col = [210 52 131]/255; %col_outside
    %         elseif(t==0)
    %             col = [128 128 128]/255; %col_nonspecific
    %         elseif(t==1)
    %             col = [59 197 124]/255; %col_shelter
    %         end               
    %         Avg_Conc_T.Subject{k,1} = Folders(nSubj).name;        
    %         Avg_Conc_T.Response_type(k,1) = t;            
    %         ROI_y = Conc_T(find(Conc_T.Subject==nSubj & Conc_T.Response_type==t),:);
    %         ROI_y.Shelter_preference_ratio(~isfinite(ROI_y.Shelter_preference_ratio)) = nan;
    %         Avg_Conc_T.Shelter_preference_ratio_x(k,1:numel(x)) = x;
    %         y = ROI_y.Shelter_preference_ratio;
    %         y(~isfinite(y)) = nan;
    %         % dy = y - mean(y(:,1:2),2);
    %         dy = y ;
    % 
    %         Avg_Conc_T.Shelter_preference_ratio(k,1:numel(mean(dy,1,'omitnan'))) = mean(dy,1,'omitnan');
    %         Avg_Conc_T.Shelter_preference_ratio_sem(k,1:numel(nansem(dy)))= nansem(dy);
    %         k = k+1;                  
    %     end
    % end
    ylim([-2, 2])
    % % Avg_Conc_T = struct2table(Avg_Conc_T);


    Avg_SPI_Change = table();
    
    k = 1;
    for(nSubj=1:NSubj)                    
        Avg_SPI_Change.Subject{k,1} = Folders(nSubj).name;          
        ROI_y = Conc_T(find(Conc_T.Subject==nSubj),:);
        ROI_y.Shelter_preference_ratio(~isfinite(ROI_y.Shelter_preference_ratio)) = nan;            
        Avg_SPI_Change.All(k,1) = mean(mean(ROI_y.Shelter_preference_ratio(:,9:10),1,'omitnan'),2,'omitnan')-mean(mean(ROI_y.Shelter_preference_ratio(:,1:2),1,'omitnan'),2,'omitnan');                           
        k = k+1;    
    end
    k = 1;
    for(nSubj=1:NSubj)                    
        Avg_SPI_Change.Subject{k,1} = Folders(nSubj).name;        
        for(t=-1:1)
            ROI_y = Conc_T(find(Conc_T.Subject==nSubj & Conc_T.Response_type==t),:);
            % Avg_SPI_Change.Outside_NS_shelter(k,t+2) = mean(sum(ROI_y.Shelter_preference_ratio_change,2,'omitnan'));
            last_y = [];
            for(i=1:size(ROI_y,1))
                last_y(i,1) = ROI_y.Shelter_preference_ratio_change(i,find(~isnan(ROI_y.Shelter_preference_ratio_change(i,:)),1,'last'));
            end
            Avg_SPI_Change.Outside_NS_shelter(k,t+2) = mean(last_y);
        end
        k = k+1;
    end

    cd(Base_folder)
    save('Avg_Shelter_preference_ratio_over_time.mat', "Avg_Conc_T","Avg_SPI_Change")
    % Conc_Shelter_preference_Table.C_Shelter_preference_ratio = C_Shelter_preference_ratio;
end


% 
% 
% fig_1 = figure;
% set(gcf,'color','w','Position', [200, 200, 400, 400]);
% x = [300:300:3900]/60;
% y = Conc_T.Shelter_preference_ratio;
% y(~isfinite(y)) = nan;
% y = y - y(:,1);
% for(t=-1:1:1)
%     ROI_y = y(find(Conc_T.Response_type==t),:);
%     if(t==-1)
%         col = [210 52 131]/255; %col_outside
%     elseif(t==0)
%         col = [128 128 128]/255; %col_nonspecific
%     elseif(t==1)
%         col = [59 197 124]/255; %col_shelter
%     end    
% %     subplot(131)
% %     for(i=1:size(ROI_y,1))
% %         plot(x, ROI_y(i,:),'color',col); hold on;
% %     end
% %     xlim([0,55])
% %     set(gca,'box','off','xtick',[0:10:50],'tickdir','out')    
% %     subplot(132)    
%     shadedErrorBar(x(1:10), mean(ROI_y(:,1:10),'omitnan'),nansem(ROI_y(:,1:10)),'lineProps',{'color',col,'markerfacecolor',col}); hold on;
% %     shadedErrorBar(x,  'color',col); hold on;
%     ylim([-0.7 1.2])
%     xlim([0,55])
%     set(gca,'box','off','xtick',[0:10:50],'tickdir','out')
% end
% Fig_data.Shelter_preference_ratio_over_time.x = x(1:10);
% Fig_data.Shelter_preference_ratio_over_time.y = nanmean(ROI_y(:,1:10));
% Fig_data.Shelter_preference_ratio_over_time.y_sem = nansem(ROI_y(:,1:10));
% xlabel('Time (min)')
% ylabel('Shelter preference index')
% cd(Base_folder)
% saveas(fig_1,'Shelter_preference_ratio_over_time.pdf')
% 
% %% Stacked chart
% col_outside = [210 52 131]/255; 
% col_shelter = [59 197 124]/255; 
% col_nonspecific = [128 128 128]/255;
% 
% Shelter_Neurons_Population_T = [numel(find(Conc_T.Response_type==-1)) numel(find(Conc_T.Response_type==0)) numel(find(Conc_T.Response_type==1))]
% % figure
% % set(gcf,'color','w')
% % subplot(121)
% % b1 = bar(1,table2array(Shelter_Neurons_Population_T(1,:)),'stacked','Facecolor','flat');
% % b1(1).BarWidth = 0.25;
% % b1(1).CData = col_outside;
% % b1(2).CData = col_nonspecific;
% % b1(3).CData = col_shelter;
% % set(gca,'ydir','reverse','box','off')
% % subplot(122)
% 
% fig_2 = figure
% b2 = bar(2,Shelter_Neurons_Population_T/sum(Shelter_Neurons_Population_T)*100,'stacked','Facecolor','flat');
% b2(1).BarWidth = 0.25;
% b2(1).CData = col_outside;
% b2(2).CData = col_nonspecific;
% b2(3).CData = col_shelter;
% set(gca,'ydir','reverse','box','off')
% sgtitle(Group)
% 
% cd(Base_folder)
% saveas(fig_2,'Shelter_preference_population_bar.pdf')
% 
% fig_3 = figure
% set(gcf,'color','w','Position', [500, 200, 450, 450]);
% ax = gca();
% h = pie(ax,[Shelter_Neurons_Population_T]);
% ax.Colormap = [col_outside;col_nonspecific;col_shelter]; 
% % Create legend
% labels = {'Outside','Non-specific','Shelter'};
% lgd = legend(labels,'Location','southoutside')
% sgtitle(Group)
% cd(Base_folder)
% saveas(fig_3,'Shelter_preference_population_chart.pdf')
% 
% 
% fig_4 = figure; 
% set(gcf,'color','w','Position', [200, 200, 450, 450]);
% x = [300:300:3900];
% y(~isfinite(y)) = nan;
% 
% Data_cdf = [];
% % col_map = zeros(size(y,2));
% col_map(1:10,:) = viridis(10);
% % for(t=1:size(y,2))
% for(t=1:10)
%     if(t<=10)
%         ROI_y = y(:,t);
%         Data_cdf{t} = gen_cdf(ROI_y);
%         plot(Data_cdf{t}.x, Data_cdf{t}.y,'color',col_map(t,:)); hold on;
%     end
% end
% colormap(col_map)
% colorbar('location','southoutside');
% xlim([-2.2 3])
% sgtitle(Group)
% set(gca,'box','off','tickdir','out','ytick',[0:0.2:1])
% 
% cd(Base_folder)
% saveas(fig_4,'Shelter_preference_ratio_over_time_cdf.pdf')
% 
% Fig_data.Shelter_preference_ratio_cdf.Data_cdf = Data_cdf;
% Fig_data.Shelter_preference_ratio_cdf.Data_cdf = Data_cdf;
% 
% save('Fig_data.mat','Fig_data')
% end
% 
% 
% 
% %     db = [];
% %     display(pwd)
% %     [Heatmap_data] = Plot_trace_heatmap_func();
% %     
% %     Rot_Ref_Maze_XY_Coords = Heatmap_data.Rot_Ref_Maze_XY_Coords;
% %     x = Heatmap_data.x;
% %     y = Heatmap_data.y;
% %     Rot_Ref_Shelter_X0Y0 = Heatmap_data.Rot_Ref_Shelter_X0Y0;
% %     Rot_Ref_Shelter_XY_Coords = Heatmap_data.Rot_Ref_Shelter_XY_Coords;
% %     xbins = Heatmap_data.xbins;
% %     ybins = Heatmap_data.ybins;
% %     Hist_array = Heatmap_data.Hist_array; 
% %     Value_array = Heatmap_data.Value_array; 
% %     Avg_array = Heatmap_data.Avg_array; 
% %     norm_filtered_zg = Heatmap_data.norm_filtered_zg;
% %     
% %     G_Heatmap_data.norm_filtered_zg(:,:,nSubj) = Heatmap_data.norm_filtered_zg;    
% %     
% %     plot(Rot_Ref_Maze_XY_Coords(:,1), Rot_Ref_Maze_XY_Coords(:,2),'k'); hold on;
% %     plot(Rot_Ref_Maze_XY_Coords(:,1)/2, Rot_Ref_Maze_XY_Coords(:,2)/2,'k--'); hold on;
% % %     plot(x, y,'color',[0.5 0.5 0.5]); hold on;
% %     plot(Rot_Ref_Shelter_X0Y0(1), Rot_Ref_Shelter_X0Y0(2),'ro'); hold on;
% %     plot(Rot_Ref_Shelter_XY_Coords(:,1), Rot_Ref_Shelter_XY_Coords(:,2),'r'); hold on;    
% % end
% % G_Heatmap_data.Avg_norm_filtered_zg = nanmean(G_Heatmap_data.norm_filtered_zg,3);
