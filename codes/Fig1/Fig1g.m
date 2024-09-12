%% Code for Fig1g
% Created by Kanghoon Jung
% Contact: kanghoon.jung AT alleninstitute.org for bug reports, questions, or feature requests.

clear all
close all
data_folder = 'C:\Users\jungk\Google Drive\Github\Jung_NatNeuro_2024\codes\Fig1\';
cd(data_folder)

% Loop over different conditions (nC)
for(nC=1:3)
    % Load data based on the current condition
    if(nC==1)
        load("Fig1g_no_shelter_data.mat")
        nC_color = [0 0 0];
    elseif(nC==2)
        load("Fig1g_shelter_data.mat")
        nC_color = [33 64 154]/255;
    elseif(nC==3)
        load("Fig1g_shelter_removed_data.mat")
        nC_color = [91 0 173]/255;
    end

    % 1) Speed    
    bin_size = 0.1;
    edges = [-2:bin_size:6];    
    x_field = 'TimefromTrigger';
    y_field = 'Speed';    
    Conc_Avg.x = [];
    Conc_Avg.y = [];       
    for(i=1:size(Conc_Tcells_re,1))
        x = Conc_Tcells_re{i,1}.(x_field);
        y = Conc_Tcells_re{i,1}.(y_field);
        Avg = Y_averaging_wrt_X_edge(x,y,edges);
        Conc_Avg.x(i,:) = Avg.X;
        Conc_Avg.y(i,:) = Avg.Y;
    end
        
    figure(1), 
    set(gcf,'color','w','position',[600 200 400 400])
    shadedErrorBar(nanmean(Conc_Avg.x),nanmean(Conc_Avg.y),nansem(Conc_Avg.y),'lineProps', {'color', nC_color});
    xlim([edges(1) edges(end)])
    yrange = [0 50];
    
    ylim(yrange)
    xlabel(x_field);
    ylabel(y_field);
    ylabel('Speed');
    set(gca,'tickdir','out','box','off','FontSize',12, 'FontName', 'Helvetica')
    set(gca,'ytick',[0:10:yrange(2)],'xtick',[ceil(edges(1)):1:edges(end)]); 
    
    % 2) Dist2Shelter
    bin_size = 0.1;
    edges = [-2:bin_size:6];
    
    x_field = 'TimefromTrigger';
    y_field = 'Dist2Shelter';
    
    Conc_Avg.x = [];
    Conc_Avg.y = [];
    for(i=1:size(Conc_Tcells_re,1))
        x = Conc_Tcells_re{i,1}.(x_field);
        y = Conc_Tcells_re{i,1}.(y_field);
        Avg = Y_averaging_wrt_X_edge(x,y,edges);
        Conc_Avg.x(i,:) = Avg.X;
        Conc_Avg.y(i,:) = Avg.Y;
    end
        
    figure(2), 
    set(gcf,'color','w','position',[600 200 400 400])
    shadedErrorBar(nanmean(Conc_Avg.x),nanmean(Conc_Avg.y),nansem(Conc_Avg.y),'lineProps', {'color', nC_color});
    xlim([edges(1) edges(end)])
    yrange = [0 50];
    ylim(yrange)
    xlabel(x_field);
    ylabel(y_field);
    ylabel('Dist2Shelter');
    set(gca,'tickdir','out','box','off','FontSize',12, 'FontName', 'Helvetica')
    set(gca,'ytick',[0:10:yrange(2)],'xtick',[ceil(edges(1)):1:edges(end)]); 
end
