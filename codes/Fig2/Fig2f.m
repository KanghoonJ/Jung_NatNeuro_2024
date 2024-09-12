% Group data analysis
% Created by Kanghoon Jung(kjung AT jhmi.edu) in Jan. 2020
% Last modified: Jan. 27, 2020
% Please contact kjung AT jhmi.edu with any bug reports, questions or feature requests for this program. 
clear all
close all
%% Data base
gp_db = [];

%% Data Acquisition information 
%% Get subfolders
Base_folder = ['H:\Project #5\Experiments\Photometry VI_jR+dL_FI@NAc\Preprocessed Dataset (JHMI&MPFI)'];
% Base_folder = 'H:\Project #5\Experiments\Photometry VI_hG6f_FI_@NAc\Preprocessed Dataset (MPFI)'
cd(Base_folder) 
Sub_Folders = dir('*');
isub = [Sub_Folders(:).isdir];
Sub_Folders = Sub_Folders(isub);
Sub_Folders(ismember({Sub_Folders.name}',{'.','..'})) = [];
Num_Sub_Folders = size(Sub_Folders,1);
%% Load all db
for(nSub = 7:7)    
% for(nSub = 1:1)    
% for(nSub = 1:Num_Sub_Folders)    
    db = [];
    tic
    cd(Base_folder) 
    Sub.ID = Sub_Folders(nSub).name; 
    Subfolder = [Sub_Folders(nSub).name];
    cd(Subfolder)
    display(Subfolder)
    %% Load db
%     ROI_db_file = dir('*Analysized.mat');
%     load(ROI_db_file(end).name)        
%     Date = datetime({(ROI_db_file(:).date)}');
%     [Sorted_Date, I] = Sort_date(Date);    
    % Most recent one
%     ROI_db_file(I(end)).name        
    % Load most recent analysis datafile
    load('Dataset.mat');
%     display(db.DAQ.Subj_id)    
    toc
end



figure, 
ax1 = subplot(311)
plot(T.Timestamp, T.DFoF_R,'r'); hold on;
ax2 = subplot(312)
plot(T.Timestamp, T.DFoF_G,'g'); hold on;
ax3 = subplot(313)
plot(T.Timestamp, T.Noldus_Dist2Shelter,'k'); hold on;
linkaxes([ax1,ax2,ax3],'x');



%% Shelter info
Maze_XY_Coords = Info.Maze_XY_Coords;
Shelter_XY_Coords = Info.Shelter_XY_Coords;
Shelter_X0Y0 = mean(Shelter_XY_Coords);
% Shelter grid
binsize = 2;
xedges = floor(min(Maze_XY_Coords(:,1))/binsize)*binsize:binsize:ceil(max(Maze_XY_Coords(:,1))/binsize)*binsize;
yedges = floor(min(Maze_XY_Coords(:,2))/binsize)*binsize:binsize:ceil(max(Maze_XY_Coords(:,2))/binsize)*binsize;
[xq, yq] = meshgrid(xedges, yedges);
Grid_log_Dist2Shelter = log(sqrt((Shelter_X0Y0(1)-xq).^2+(Shelter_X0Y0(2)-yq).^2));
Grid_Dist2Shelter = sqrt((Shelter_X0Y0(1)-xq).^2+(Shelter_X0Y0(2)-yq).^2);

[in, on] = inpolygon(xq,yq,Maze_XY_Coords(:,1), Maze_XY_Coords(:,2));
Arena_mask = double(in);
Arena_mask(find(in==0))=nan;
Arena_mask(find(in==1))=0;
Arena_mask(find(Arena_mask==0)) = nan;



%% Individual trace plot
% Trajectory
%% Plot Trajectory 3D
Timestamp = T.Timestamp;
% T_lim = [T.Timestamp(1) T.Timestamp(end)];
T_lim = [T.Timestamp(1)+30 T.Timestamp(1)+450]; 
% T_lim = [T.Timestamp(1)+30 T.Timestamp(1)+1030];
% T_lim = [T.Timestamp(1)+100 T.Timestamp(1)+1100];




ROI_Frames = Find_Frames_in_TS(T_lim,Timestamp);
ROI_T = T(ROI_Frames,:);
x = ROI_T.Noldus_Body_Coord(:,1);
y = ROI_T.Noldus_Body_Coord(:,2);
z = ROI_T.Timestamp;
Nan_Frames = unique([find(isnan(x)),find(isnan(y))]);
ROI_T(Nan_Frames,:) = [];



%% Skewness map
all_x = T.Noldus_Body_Coord(:,1);
all_y = T.Noldus_Body_Coord(:,2);
[x_bins, x_edges] = discretize(all_x(:,1),xedges,'IncludedEdge','left'); 
y_bins = discretize(all_y,yedges,'IncludedEdge','left');                    
All_Trace_Grid_z = nan(size(xq));
k = [];
for i=1:size(All_Trace_Grid_z,2)
    k = find(x_bins==i);    
    for j=1:size(All_Trace_Grid_z,1)        
        n = find(y_bins(k)==j);
%         Trace_Grid_z(j,i) = numel(n);
        if numel(k(n))>0  
            All_Trace_Grid_z(j,i) = 1;
        end
    end    
end

T_coverage_mask = nan(size(xq));
T_coverage_mask(find(~isnan(Arena_mask)|~isnan(All_Trace_Grid_z))) = 1;

x = ROI_T.Noldus_Body_Coord(:,1);
y = ROI_T.Noldus_Body_Coord(:,2);
[x_bins, x_edges] = discretize(x(:,1),xedges,'IncludedEdge','left'); 
y_bins = discretize(y,yedges,'IncludedEdge','left');                    
ROI_Trace_Grid_z = nan(size(xq));
k = [];
for i=1:size(ROI_Trace_Grid_z,2)
    k = find(x_bins==i);    
    for j=1:size(ROI_Trace_Grid_z,1)        
        n = find(y_bins(k)==j);
%         Trace_Grid_z(j,i) = numel(n);
        if numel(k(n))>0  
            ROI_Trace_Grid_z(j,i) = 1;
        end
    end    
end


Coverage_rate = numel(find(~isnan(ROI_Trace_Grid_z)))/numel(find(~isnan(T_coverage_mask)));


% % Shelter
% r = 5;
% [X,Y,Z] = cylinder(r,50);
% Acclimation_duration = 420;
% Shelter_cylinder_X = X+Shelter_X0Y0(1);
% Shelter_cylinder_Y = Y+Shelter_X0Y0(2);
% Shelter_cylinder_Z = Acclimation_duration+Z*(Timestamp(end)-Acclimation_duration);
% Shelter_cylinder_Z2 = Acclimation_duration*Z;
% Shelter_Coords = [Shelter_cylinder_X, Shelter_cylinder_Y]; 
% 
% 
% Shelter_color = [130 208 255]/255;
% Shelter_color2 = [0.9 0.9 0.9];
% 
% 
% Arena_color = [220 220 220]/255;
% Maze_XY_Coords_z = zeros(size(Maze_XY_Coords,1),1);
% Maze_XYZ = [Maze_XY_Coords, Maze_XY_Coords_z];

dot_size = 5;

close all
%% Distribution of signal on movement traces on arena
figure, 
set(gcf,'color','w','Position', [2315, -136, 1600, 1100]);
% p(Maze_XYZ(:,1), Maze_XYZ(:,2), Maze_XYZ(:,3),Arena_color,'LineStyle','none'); hold on;
alpha2 = 0.5;
Maze_Height = z(1);
Maze_XYZ = [Maze_XY_Coords, zeros(size(Maze_XY_Coords,1),1)];
Shelter_XYZ = [Shelter_XY_Coords, zeros(size(Shelter_XY_Coords,1),1)];
x = ROI_T.Noldus_Body_Coord(:,1);
y = ROI_T.Noldus_Body_Coord(:,2);
for(i=1:4)
    subplot(2,2,i)   
%     delay = 10;
    if(i==1)
        ax1 = subplot(2,2,1);
        Col_Data = ROI_T.DFoF_R;
        title('R')    
    elseif(i==2)
        ax2 = subplot(2,2,2);
        Col_Data = ROI_T.DFoF_G;
        title('G')
    elseif(i==3)
        ax3 = subplot(2,2,3);
        Col_Data = ROI_T.DFoF_RG;
        title('RG')
    elseif(i==4)
        ax4 = subplot(2,2,4);
        Col_Data = ROI_T.Noldus_Vel;
        title('Vel')
    end
%     z = zscore(Col_Data);
    z = mat2gray(Col_Data);

   % caxis([0 0.15])
%     cmap = colormap(whitejet(50));
%     cmap = colormap(parula);
    cmap = colormap(flipud(gray));
    [C, EDGES] = Convert_value2color(z, cmap);

    dot_size = 10;
    h1 = scatter3(x,y,z,dot_size, z,'filled'); hold on; 
    set(h1,'MarkerEdgeColor','none','MarkerEdgeAlpha', alpha2, 'MarkerFaceAlpha', alpha2);
    % [caz,cel] = view()
%     caz = -45; cel = 52;
    caz = -110.968; cel = 44.8656;
    view([caz,cel])
    grid on;
    xlim([-30 30])
    ylim([-30 30])
%     zlim(T_lim);
    caxis([EDGES(1) EDGES(end)])
    

    % Tick control
    zt = zticks; 
    zticks([zt(1),zt(end)]);
    xt = xticks; 
    xticks([xt(1):10:xt(end)]);
    yt = yticks; 
    yticks([yt(1):10:yt(end)]);
    colormap(cmap);
%     caxis([0 1])
    cbh = colorbar;
    ct = cbh.Ticks;
    cbh.Ticks = [ct(1),ct(end)];

%     if(i==1)
%     Fig_title = ['T: ' num2str(j) ', t:' num2str(z(1)) '-' num2str(z(end))];
%     title(Fig_title)
%     end
    drawnow
    ax = gca;
    ax.GridColor = [0.8, 0.8, 0.8];  % [R, G, B]
    
    % caxis([0 EDGES(end)])
   plot3(Maze_XYZ(:,1), Maze_XYZ(:,2), Maze_XYZ(:,3),'k'); hold on;    
   plot3(Shelter_XYZ(:,1),Shelter_XYZ(:,2), Shelter_XYZ(:,3),'k'); hold on;      
end
Link = linkprop([ax1, ax2, ax3, ax4],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

set(gcf, 'renderer', 'painters');

    cmap = colormap(parula);

