
function Avg = Y_averaging_wrt_X_edge(x,y,edges)
%%Test
% x = G_DATA.Corr_Dist.Conc_Data{1}(:,4);
% y = G_DATA.Corr_Dist.Conc_Data{1}(:,3);  
% bin_size = 30;
% edges = [floor(min(x)/bin_size)*bin_size:50:ceil(max(x)/bin_size)*bin_size];

x_group = discretize(x,edges,'IncludedEdge','left');
% x_group = discretize(x,edges,'IncludedEdge','right');


Avg = [];
Avg.X = [edges(1:end-1) + (edges(2)-edges(1))/2]';

for(i=1:numel(edges)-1)    
    Avg.Y(i,1) = nanmean(y(find(x_group==i)));
    Avg.Y_sem(i,1) = nanstd(y(find(x_group==i)),0,1)/sqrt(numel(y(find(x_group==i))));    
end

% figure, 
% % plot(x,y,'.'); hold on;
% errorbar(Avg.X, Avg.Y,Avg.Y_sem);
%     