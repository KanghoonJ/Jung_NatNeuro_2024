% Find_Frames in Timestamp
function ROI_Frames = Find_Frames_in_TS(t_range, TS)
ROI_Frames = find(TS>=t_range(1) & TS<=t_range(2));
