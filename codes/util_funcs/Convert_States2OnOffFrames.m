%% Convert State value to ON-OFF Frames
% Created by Kanghoon Jung
% Please contact kanghoon.jung AT alleninstitute.org with any bug reports, questions or feature requests for this program. 

function [OnOff_Frames, Adjusted_States] = Convert_States2OnOffFrames(States,minimum_frame_num)
%% On: 1
%% OFF: 0
% States = Threat_onoff;
% minimum_frame_num = 15;

Runlength = diff(find([1;diff(States);1]));
State_type = States(find([1;diff(States)]));
State = table(State_type,Runlength);
State.Cum_Runlength = cumsum(Runlength);
Adjusted_State = State;
Adjusted_State(find(State.Runlength<minimum_frame_num),:) = [];
if(size(Adjusted_State,1)>0)
    Adjusted_State(find(diff(Adjusted_State.State_type)==0),:) = [];
    Adjusted_State.Onset_frame = [1;Adjusted_State.Cum_Runlength(1:end-1)+1];
    Adjusted_State.Offset_frame = Adjusted_State.Cum_Runlength;
    % Adjusted_State.Onset_time = Timestamp(Adjusted_State.Onset_frame);
    % Adjusted_State.Offset_time = Timestamp(Adjusted_State.Offset_frame);
    OnOff_Frames = table2array([Adjusted_State(find(Adjusted_State.State_type>0),4:5)]);
    Adjusted_States = zeros(numel(States),1);
    for(i=1:size(OnOff_Frames,1))
        Adjusted_States(OnOff_Frames(i,1):OnOff_Frames(i,2),1) = 1;
    end
else
    OnOff_Frames = [];
    Adjusted_States = [];
end


% Old algorithm
% States = db.Behavior_data.Noldus.Track_data.In_zone;
% On_Frames = find(diff(States)==1)+1;
% Off_Frames = find(diff(States)==-1);
% if(States(1,1)==1)
%     On_Frames = [1; On_Frames];
% end
% if(States(end,1)==1)
%     Off_Frames = [Off_Frames; size(States,1)];
% end
% On_Frames(:,2) = ones(size(On_Frames,1),1);
% Off_Frames(:,2) = zeros(size(Off_Frames,1),1);
% 
% OnOff_Frames = On_Frames(:,1);
% for i=1:size(OnOff_Frames,1)
% %     if(i<=size(OnOff_Frames,1))
%     OnOff_Frames(i,2) = Off_Frames(find(Off_Frames(:,1)>On_Frames(i,1),1),1);
% %     end
% end
% Overlap = zeros(size(OnOff_Frames,1),1);
% for i=2:size(OnOff_Frames,1)
%     if(OnOff_Frames(i,2)==OnOff_Frames(i-1,2))
%         Overlap(i,1) = 1;
%     end
% end
% OnOff_Frames(find(Overlap==1),:) = [];
% 
% OnOff_Frames(diff(OnOff_Frames,1,2)<minimum_frame_num,:) = [];



