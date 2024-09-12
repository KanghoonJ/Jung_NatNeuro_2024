%% Plot shaded errorbar
% Created by Kanghoon Jung(kjung AT jhmi.edu) in 2019
% Last modified: Oct. 30, 2019
% Please contact kjung AT jhmi.edu with any bug reports, questions or feature requests for this program.


% % Data = [];
% x = Data(:,1);
% y1 = Data(:,2);
% y1_sem = Data(:,3);
% y2 = Data(:,4);
% y2_sem = Data(:,5);
% y3 = Data(:,6);
% y3_sem = Data(:,7);

% figure, 
% set(gcf,'color','w','position',[100 100 300 300])
% shadedErrorBar(x,y1,y1_sem,'lineProps', '-g'); hold on;
% shadedErrorBar(x,y2,y2_sem,'lineProps', '-r'); hold on;
% shadedErrorBar(x,y3,y3_sem,'lineProps', '-k'); hold on;
% set(gca,'tickdir','out','xlim',[0 55],'ylim',[-0.6 0.9],'box','off')



function Plot_shaded_errorbar(var)
x = var.x;
y = var.y;
y_sem = var.y_sem;
l_col = var.color;
% Timestamp = Track_T.Timestamp;
% Dist2Shelter = Track_T.Dist2Shelter;
% x = S{i,1}.ROI_Track_T.ROI_Timestamp_adjusted;
% y = S{i,1}.ROI_Track_T.Dist2Shelter;
% x_line = Shelter_in_time;

% figure
set(gcf,'color','w');
shadedErrorBar(x,y,y_sem,'lineProps',{'color',l_col,'markerfacecolor',l_col}); hold on;
% shadedplot(x,y); hold on;
% if(~isempty(Shelter_in_time))
%     y_line = Shelter_in_Dist2Shelter;
%     plot([x(1), x(end)], [y_line, y_line],'color',[0.3 0.3 0.3]); hold on;
%     yl = ylim;
%     plot([x_line, x_line], [yl(1), yl(end)],'color',[0.3 0.3 0.3]); hold on;
% end
set(gca,'box','off','tickdir','out')
xlabel('x')
ylabel('y')
% ylabel('Dist2Shelter')
