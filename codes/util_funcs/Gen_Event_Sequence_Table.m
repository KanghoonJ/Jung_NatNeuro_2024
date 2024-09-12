% Created by Kanghoon Jung
% Please contact kanghoon.jung AT alleninstitute.org with any bug reports, questions or feature requests for this program. 

function Sorted_Event_Sequence_Table = Gen_Event_Sequence_Table(db)
    %% Target Event
    Track_T = db.Behavior_data.Behavior_data_table;

    %% Shelter inout
    minimum_frame_num = 15;
    [OnOff_Frames, Adjusted_States] = Convert_States2OnOffFrames(Track_T.Shelter_inout_Adjusted,minimum_frame_num);

    ShelterIn_Event_Sequence_Table = table();
    ShelterIn_Event_Sequence_Table.Event_type = repmat({'Shelter_in'},size(OnOff_Frames,1),1);
    ShelterIn_Event_Sequence_Table.Event_time = Track_T.Timestamp(OnOff_Frames(:,1));

    ShelterOut_Event_Sequence_Table = table();
    ShelterOut_Event_Sequence_Table.Event_type = repmat({'Shelter_out'},size(OnOff_Frames,1),1);
    ShelterOut_Event_Sequence_Table.Event_time = Track_T.Timestamp(OnOff_Frames(:,2));

    %% Threat onoff
    [OnOff_Frames, Adjusted_States] = Convert_States2OnOffFrames(Track_T.Threat_onoff,minimum_frame_num);

    ThreatOn_Event_Sequence_Table = table();
    ThreatOn_Event_Sequence_Table.Event_type = repmat({'Threat_on'},size(OnOff_Frames,1),1);
    ThreatOn_Event_Sequence_Table.Event_time = Track_T.Timestamp(OnOff_Frames(:,1));

    ThreatOff_Event_Sequence_Table = table();
    ThreatOff_Event_Sequence_Table.Event_type = repmat({'Threat_off'},size(OnOff_Frames,1),1);
    ThreatOff_Event_Sequence_Table.Event_time = Track_T.Timestamp(OnOff_Frames(:,2));

    %% Movement onoff
    [OnOff_Frames, Adjusted_States] = Convert_States2OnOffFrames(Track_T.Movement_Adjusted,minimum_frame_num);

    MovementOn_Event_Sequence_Table = table();
    MovementOn_Event_Sequence_Table.Event_type = repmat({'Movement_on'},size(OnOff_Frames,1),1);
    MovementOn_Event_Sequence_Table.Event_time = Track_T.Timestamp(OnOff_Frames(:,1));

    MovementOff_Event_Sequence_Table = table();
    MovementOff_Event_Sequence_Table.Event_type = repmat({'Movement_off'},size(OnOff_Frames,1),1);
    MovementOff_Event_Sequence_Table.Event_time = Track_T.Timestamp(OnOff_Frames(:,2));

    %% Flight onoff
    [OnOff_Frames, Adjusted_States] = Convert_States2OnOffFrames(Track_T.Flight_Adjusted,minimum_frame_num);

    FlightOn_Event_Sequence_Table = table();
    FlightOn_Event_Sequence_Table.Event_type = repmat({'Flight_on'},size(OnOff_Frames,1),1);
    FlightOn_Event_Sequence_Table.Event_time = Track_T.Timestamp(OnOff_Frames(:,1));

    FlightOff_Event_Sequence_Table = table();
    FlightOff_Event_Sequence_Table.Event_type = repmat({'Flight_off'},size(OnOff_Frames,1),1);
    FlightOff_Event_Sequence_Table.Event_time = Track_T.Timestamp(OnOff_Frames(:,2));


    Event_Sequence_Table = [ShelterIn_Event_Sequence_Table; ShelterOut_Event_Sequence_Table; ThreatOn_Event_Sequence_Table; ThreatOff_Event_Sequence_Table; MovementOn_Event_Sequence_Table; MovementOff_Event_Sequence_Table; FlightOn_Event_Sequence_Table; FlightOff_Event_Sequence_Table];

    Sorted_Event_Sequence_Table = sortrows(Event_Sequence_Table,2)
end