% Event_Analysis
close all
clear all
% Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_NpHR dPAG\m139_111322_Trial 81\Events\Threat';
% Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_NpHR dPAG\m139_111322_Trial 81\Events\589t +Threat (together)';
% Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_NpHR dPAG_Events_only\589t (first) - 473t (second)';
% Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_NpHR dPAG_Events_only\589t';
Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_NpHR dPAG_Events_only\589t +Threat (together)';
% Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_NpHR dPAG_Events_only\Threat';
% Flight_only = 1;

% Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_ChrimsonR Threats_Events_only\589t';
% Flight_only = 0;

% Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_ChrimsonR Threats_Events_only\589t +Threat (together)';
% Flight_only = 1;
% pl_color = [251 176 64]/255;

% Base_folder = 'H:\Project #5\Revision Experiments\Cal-light probe_ChrimsonR Threats_Events_only\Threat';
Flight_only = 0;
pl_color = [0.5 0.5 0.5];


foldername = strsplit(Base_folder,"\");
Event_type_name = cell2mat(foldername(end));

cd(Base_folder)
A = dir;
A = A(~ismember({A.name},{'.','..'}));
Folders = A([A.isdir]);
Num_Folders = size(Folders,1);

Avg_Head_Shelter_Angle = [];
Conc_Goal_Angle = [];
Conc_Speed = [];
Conc_ROI_Norm_Dist2Shelter_T = [];

% for nFolder = 1:Num_Folders
%     cd([Base_folder '\' Folders(nFolder).name])    
%     display(pwd)    
%     load('db_Event.mat')  
%     fig1 = figure(1);        
%     Body_Coords = db_Event.ROI_Track_T_Threat.Body_Coords;
%     plot(Body_Coords(:,1),Body_Coords(:,2),'k');
%     xlim([-30 30])
%     ylim([-30 30])
%     title(Folders(nFolder).name)    
%     saveas(fig1,'Trace.pdf')
%     close(gcf);    
% end
    
xrange = [-0.25 4];
set(gcf, 'renderer', 'painters'); % for vector
for nFolder = 1:Num_Folders
    cd([Base_folder '\' Folders(nFolder).name])    
    display(pwd)    
    load('db_Event.mat')  
    fig10 = figure(10);
    x = db_Event.ROI_Track_T.ROI_Timestamp_adjusted;
    y = db_Event.ROI_Track_T.Dist2Shelter;
    plot(x,y,'k');
    xlim(xrange)
    ylim([0 55])
    xlabel('Time from Stim onset (s)');
    ylabel('Distance to shelter (cm)');
    saveas(fig10,'XCorr_Dist.pdf')
    close(gcf);
    
    fig1 = figure(1), 
    set(gcf,'color','w','Position', [40, 660, 510, 300]);
    Plot_XCorr_Dist_Speed(db_Event,Flight_only)
    title(Event_type_name)
    
    fig2 = figure(2), 
    set(gcf,'color','w','Position', [570, 430, 300, 530]);
    Plot_Norm_Traces(db_Event)
    title(Event_type_name)
    
    fig3 = figure(3), 
    set(gcf,'color','w','Position', [780, 660, 180, 300]);
    Plot_XCorr_Norm_Dist(db_Event,Flight_only)
    title(Event_type_name)
    
    Event_onset_Frame = find(db_Event.ROI_Track_T.ROI_Timestamp_adjusted>=0,1);
    Event_offset_Frame = min([find(db_Event.ROI_Track_T.ROI_Timestamp_adjusted>=3,1), find(db_Event.ROI_Track_T.Dist2Shelter<=5,1)]);
%     ROI_Frames = [Event_onset_Frame:Event_offset_Frame];

    Conc_ROI_Norm_Dist2Shelter_T.Time(1:size(db_Event.ROI_Track_T,1),nFolder) = db_Event.ROI_Track_T.ROI_Timestamp_adjusted;
    Conc_ROI_Norm_Dist2Shelter_T.Norm_Dist2Shelter_T(1:size(db_Event.ROI_Track_T,1),nFolder) = db_Event.ROI_Track_T.Dist2Shelter/db_Event.ROI_Track_T.Dist2Shelter(Event_onset_Frame);
    
    
    
    if(Flight_only)
        if(find(db_Event.ROI_Track_T.Flight_Adjusted>0))    
            disp('Flight detected')
            ROI_Flight_Frames = db_Event.ROI_Flight_Frames;    
            ROI_Flight_Frames = [ROI_Flight_Frames(1):Event_offset_Frame];

            Conc_Goal_Angle = [Conc_Goal_Angle; db_Event.ROI_Track_T.Goal_Angle(ROI_Flight_Frames)];
            Conc_Speed = [Conc_Speed; db_Event.ROI_Track_T.Speed(ROI_Flight_Frames)];

            Avg_Head_Shelter_Angle(nFolder,1) = nanmean(abs(db_Event.ROI_Track_T.Goal_Angle(ROI_Flight_Frames)));
        else            
            Avg_Head_Shelter_Angle(nFolder,1) = nan;
        end
    else        
        ROI_Frames = [Event_onset_Frame:Event_offset_Frame];        
        Conc_Goal_Angle = [Conc_Goal_Angle; db_Event.ROI_Track_T.Goal_Angle(ROI_Frames)];
        Conc_Speed = [Conc_Speed; db_Event.ROI_Track_T.Speed(ROI_Frames)];

        Avg_Head_Shelter_Angle(nFolder,1) = nanmean(abs(db_Event.ROI_Track_T.Goal_Angle(ROI_Frames)));        
    end
%     subplot(133)
%     y1 = db_Event.ROI_Track_T.Speed;
%     col = [0.5 0.5 0.5];
%     plot(x,y1,'color',col); hold on;
%     plot(x(end),y1(end),'o','color',col); hold on;
%     set(gca,'box','off','tickdir','out')    
%     xlabel('Time from Stim onset (s)');
%     ylabel('');
end

pos = get(gcf, 'Position') %// gives x left, y bottom, width, height

%% Plot_Avg_XCorr_Norm_Dist
x = nanmean(Conc_ROI_Norm_Dist2Shelter_T.Time,2);      
y = nanmean(Conc_ROI_Norm_Dist2Shelter_T.Norm_Dist2Shelter_T,2);      
y_sem = std(Conc_ROI_Norm_Dist2Shelter_T.Norm_Dist2Shelter_T,0,2)/sqrt(size(Conc_ROI_Norm_Dist2Shelter_T.Norm_Dist2Shelter_T,2));      
p_color = [0.5 0.5 0.5];
fig4 = figure; 
set(gcf,'color','w','Position', [970, 660, 180, 300]);
shadedErrorBar(x,y,y_sem,'lineProps',{'color',p_color,'markerfacecolor',p_color});   
ylim([0 1.3])
xlim([0 4])
set(gca,'box','off','tickdir','out','xtick',[0:1:4],'ytick',[0:0.2:1.2])    
xlabel('Time from Stim onset (s)');
ylabel('Norm. Distance to Shelter (cm)');

%% Plot WindRose w/ Goal angle during events
D = Conc_Goal_Angle; % My reference is North = 0º, East = -90º.
V = Conc_Speed;  
[figure_handle,count,speeds,directions,Table] = WindRose_KJ(D,V,'anglenorth',0,'angleeast',-90,'labels',{'0�','-90�','180�','90�'},'freqlabelangle',45);
close(gcf)
edges = deg2rad([-5:10:355]);
fig5 = figure;
pax = polaraxes;
set(gcf,'color','w','Position', [40, 80, 520, 480]);
% polarhistogram(D,edges,'FaceColor','red','FaceAlpha',0.3);
polarhistogram('BinEdges',edges,'BinCounts',sum(count,2),'FaceColor',pl_color)
set(gca,'ThetaTickLabel',[]); 
sgtitle(Event_type_name)
pax.ThetaDir = 'counterclockwise';
pax.ThetaZeroLocation = 'top';
pax.FontSize = 18;
pax.RAxisLocation = 315;
thetaticks(0:30:360)
rticks([4:4:18])
rtickformat('percentage')
set(gcf, 'renderer', 'painters');
% title('Avg. goal angle dist')



%% Save result figures
cd(Base_folder)
saveas(fig1,'XCorr_Dist_Speed.pdf')
saveas(fig2,'Norm_Trace.pdf')
saveas(fig3,'XCorr_Norm_Dist.pdf')
saveas(fig4,'XCorr_Avg_Norm_Dist.pdf')
saveas(fig5,'WindRose.pdf')
saveas(fig5,'WindRose.pdf')

filename='Avg_Head_Shelter_Angle.xls';
writematrix(Avg_Head_Shelter_Angle, filename);

function Plot_XCorr_Dist_Speed(db_Event,Flight_only) % Plot Distance to Shelter and Speed
    xrange = [-0.25 4];
    if(Flight_only)
        if(find(db_Event.ROI_Track_T.Flight>0))
        ROI_Frame_start = 1;
        ROI_Frame_end = db_Event.ROI_Flight_Frames(end);
        ROI_Frames = [ROI_Frame_start:ROI_Frame_end];

    %     x = db_Event.ROI_Track_T.Timestamp(ROI_Flight_Frames) - db_Event.Event_time;  
        x = db_Event.ROI_Track_T.Timestamp(ROI_Frames) - db_Event.ROI_Track_T.Timestamp(db_Event.ROI_Flight_Frames(1));  
        subplot(121)
        y1 = db_Event.ROI_Track_T.Dist2Shelter(ROI_Frames);
        col = [0.5 0.5 0.5];
        plot(x,y1,'color',col); hold on;
%         plot(x(end),y1(end),'o','color',col,'MarkerFaceColor',col); hold on;
        xlim(xrange)
        ylim([0 53])    
        set(gca,'box','off','tickdir','out','ytick',[0:10:50])    
        xlabel('Time from Stim onset (s)');
        ylabel('Distance to Shelter (cm)');

        subplot(122)
        y2 = db_Event.ROI_Track_T.Speed(ROI_Frames);
        col = [0.5 0.5 0.5];
        plot(x,y2,'color',col); hold on;
%         plot(x(end),y2(end),'o','color',col,'MarkerFaceColor',col); hold on;
        xlim(xrange)
        ylim([0 55])
        set(gca,'box','off','tickdir','out')    
        xlabel('Time from Stim onset (s)');
        ylabel('Speed(cm/s)');
        end
    else        
        ROI_Frames = [1:size(db_Event.ROI_Track_T,1)];            
    %     x = db_Event.ROI_Track_T.Timestamp(ROI_Flight_Frames) - db_Event.Event_time;  
        x = db_Event.ROI_Track_T.ROI_Timestamp_adjusted(ROI_Frames);
        subplot(121)
        y1 = db_Event.ROI_Track_T.Dist2Shelter(ROI_Frames);
        col = [0.5 0.5 0.5];
        plot(x,y1,'color',col); hold on;
%         plot(x(end),y1(end),'o','color',col,'MarkerFaceColor',col); hold on;
        xlim(xrange)
        ylim([0 53])    
        set(gca,'box','off','tickdir','out','ytick',[0:10:50])    
        xlabel('Time from Stim onset (s)');
        ylabel('Distance to Shelter (cm)');

        subplot(122)
        y2 = db_Event.ROI_Track_T.Speed(ROI_Frames);
        col = [0.5 0.5 0.5];
        plot(x,y2,'color',col); hold on;
        plot(x(end),y2(end),'o','color',col,'MarkerFaceColor',col); hold on;
        xlim(xrange)
        ylim([0 55])
        set(gca,'box','off','tickdir','out')    
        xlabel('Time from Stim onset (s)');
        ylabel('Speed(cm/s)');
    end
end


function Plot_Norm_Traces(db_Event) % Plot Normalized Traces   
    col = [0.5 0.5 0.5];
    Event_onset_Frame = find(db_Event.ROI_Track_T.ROI_Timestamp_adjusted>=0,1);
    Event_offset_Frame = min(find(db_Event.ROI_Track_T.ROI_Timestamp_adjusted>=3,1), find(db_Event.ROI_Track_T.Shelter_inout_Adjusted==1,1));
    ROI_Frames = [Event_onset_Frame:Event_offset_Frame];

    
%     if(find(db_Event.ROI_Track_T.Flight>0))
%     ROI_Flight_Frame_start = min(union(find(db_Event.ROI_Track_T.Flight==1,1), find(diff(db_Event.ROI_Track_T.Flight)==1,1)));
%     ROI_Flight_Frame_end = max(union(find(diff(db_Event.ROI_Track_T.Flight)==-1,1,'last'), find(db_Event.ROI_Track_T.Flight==1,1,'last')));
%     ROI_Flight_Frames = [ROI_Flight_Frame_start:ROI_Flight_Frame_end];
% %     ROI_Frames = [ROI_Flight_Frame_start:size(db_Event.ROI_Track_T,2)];
%     ROI_Frames = ROI_Flight_Frames;
    
    x = db_Event.Norm_Trace.Rot_ROI_Norm_Vector(ROI_Frames,1);
    y = 1-db_Event.Norm_Trace.Rot_ROI_Norm_Vector(ROI_Frames,2);
    plot(x,y,'color',col); hold on;        
    
%     title(Event_Types_set{nT})
%     legend(Event_Types_set,'Location','NorthEastOutside')
    plot(0,0,'ks'); hold on; 
    % line([-1.0 -0.9], [-0.5 -0.5],'color','k'); hold on;
    set(gca,'box','off','tickdir','out','ytick',[0, 0.5, 1],'YDir','reverse')
    ylabel('Normalized Distance')
    xlim([-1.1, 1.1])    
    ylim([-0.3, 1.4])    
    %     set(gca,'box','off','tickdir','out','ytick',[0, 0.5, 1])
%     ylabel('Normalized Distance')
%     ylim([-0.4, 1.4])    
    
%     end
end

function Plot_XCorr_Norm_Dist(db_Event,Flight_only) % Plot Norm Distance to Shelter
    if(Flight_only)
        if(find(db_Event.ROI_Track_T.Flight>0))
            ROI_Frame_start = 1;
            ROI_Frame_end = db_Event.ROI_Flight_Frames(end);
            ROI_Frames = [ROI_Frame_start:ROI_Frame_end];            
            x = db_Event.ROI_Track_T.Timestamp(ROI_Frames) - db_Event.ROI_Track_T.Timestamp(db_Event.ROI_Flight_Frames(1));  
            y1 = db_Event.ROI_Track_T.Dist2Shelter(ROI_Frames)/db_Event.ROI_Track_T.Dist2Shelter(db_Event.ROI_Flight_Frames(1));

            col = [0.5 0.5 0.5];
            plot(x,y1,'color',col); hold on;
            plot(x(end),y1(end),'o','color',col,'MarkerFaceColor',col); hold on;
            ylim([0 1.3])
            xlim([0 4])
            set(gca,'box','off','tickdir','out','xtick',[0:1:4],'ytick',[0:0.2:1.2])    
            xlabel('Time from Stim onset (s)');
            ylabel('Norm. Distance to Shelter (cm)');
        end
    else        
        ROI_Frames = [1:size(db_Event.ROI_Track_T,1)];            
        x = db_Event.ROI_Track_T.ROI_Timestamp_adjusted;
        Stim_onset_Frame = find(db_Event.ROI_Track_T.ROI_Timestamp_adjusted>=0,1);
        y1 = db_Event.ROI_Track_T.Dist2Shelter(ROI_Frames)/db_Event.ROI_Track_T.Dist2Shelter(Stim_onset_Frame);
        
        col = [0.5 0.5 0.5];
        plot(x,y1,'color',col); hold on;
        plot(x(end),y1(end),'o','color',col,'MarkerFaceColor',col); hold on;
        ylim([0 1.3])
        xlim([0 4])
        set(gca,'box','off','tickdir','out','xtick',[0:1:4],'ytick',[0:0.2:1.2])    
        xlabel('Time from Stim onset (s)');
        ylabel('Norm. Distance to Shelter (cm)');
    end
    
    
end




%% Xcorr for Distance to Shelter aligned to Stim onset
% Plot_X_corr_x = Timestamp;   
% Plot_X_corr_Event_data = Events;
% Plot_type = 'raw';
% Plot_X_corr_ROI_X_corr_window = [-1 10];
% Plot_X_corr_sampling_freq = 1/mean(diff(Timestamp));

% for(i=1:3)
%     if(i==1)
%         Field_name = ['XCorr_Dist2Shelter'];
%         Plot_X_corr_y = Dist2Shelter;
%         fig_title = 'Dist2Shelter';
%     elseif(i==2)
%         Field_name = ['XCorr_Speed'];
%         Plot_X_corr_y = Speed;
%         fig_title = 'Speed';
%     elseif(i==3)
%         Field_name = ['XCorr_Goal_Angle'];
%         Plot_X_corr_y = Goal_Angle;
%         fig_title = 'Goal Angle';
%     end
%     [db] = Plot_X_corr_v2_basic(db, Plot_X_corr_x, Plot_X_corr_y, Plot_X_corr_Event_data, Plot_X_corr_ROI_X_corr_window, Plot_X_corr_sampling_freq, Field_name, Plot_type);
%     XCorr.x =  db.Plot.(Field_name).Data_x;
%     XCorr.y = db.Plot.(Field_name).Data_y.raw;
%     subplot(1,3,i)
% %     col = 'k';
%     for(j=1:size(XCorr.y,2))
%         col = Event_cols(Event_Types(j),:);
%         plot(XCorr.x,XCorr.y(:,j),'color',col); hold on;
%         plot(XCorr.x(end),XCorr.y(end,j),'o','color',col); hold on;
%     end
%     xlim(Plot_X_corr_ROI_X_corr_window)
%     xlabel('Time from Stim onset (s)');
%     ylabel(fig_title);
%     title(fig_title);
%     set(gca,'box','off','tickdir','out')     
% end