%% EDF1m - Analysis Script
% This script analyzes speed data during exploration and escape flight conditions.

% Close all open figures and clear the workspace to prevent conflicts from previous sessions
% close all;
% clear all;
% 
% %% Initialization
% % Define the folder where the data is stored
% data_folder = 'C:\Users\jungk\Google Drive\Github\Jung_NatNeuro_2024\codes\EDF1\';
% 
% % Change the current directory to the data folder where the .mat file is located
% cd(data_folder);
% 
% % Load the data from the .mat file (make sure "Conc_T" table is part of the data)
% load("EDF1m_data.mat");

%% Speed distribution analysis
% Extract speed values from the table "Conc_T"
Speed = Conc_T.Speed;

% Define the frames for different conditions based on movement, threat, and shelter status

% 1. Exploration frames: Animal is moving, there is no threat, and it's outside the shelter
Exploration_Frames = find(Conc_T.Movement_Adjusted == 1 & Conc_T.Threat_onoff == 0 & Conc_T.Shelter_inout_Adjusted == 0);

% 2. Escape flight frames during threat: Animal is moving, threat is present, and it's outside the shelter
Escape_Flight_Frames_Threat = find(Conc_T.Movement_Adjusted == 1 & Conc_T.Threat_onoff > 0 & Conc_T.Shelter_inout_Adjusted == 0);

% 3. Escape flight frames during non-threat: Similar to exploration but added for completeness
Escape_Flight_Frames_Non_Threat = find(Conc_T.Movement_Adjusted == 1 & Conc_T.Threat_onoff == 0 & Conc_T.Shelter_inout_Adjusted == 0);

% Extract speed values for the identified frames
Speed_Exploration = Speed(Exploration_Frames);           % Speed during exploration
Speed_Escape_Flight_Threat = Speed(Escape_Flight_Frames_Threat);  % Speed during escape flight with threat

%% Plot Histograms of Speed Distributions
% Create a figure with two subplots for histograms

figure;

% Subplot 1: Histogram for exploration speed
subplot(1, 2, 1);  % Create a 1x2 grid of subplots, position at (1)
hist(Speed_Exploration);  % Plot histogram of exploration speed
title('Exploration Speed');  % Add a title

% Subplot 2: Histogram for escape flight (threat) speed
subplot(1, 2, 2);  % Position at (2)
hist(Speed_Escape_Flight_Threat);  % Plot histogram of escape flight (threat) speed
title('Escape Flight Threat Speed');  % Add a title

%% Generate and Plot CDF (Cumulative Distribution Function) for comparison
% Generate CDFs for exploration and threat escape speeds using the "gen_cdf" function

% Generate CDF for exploration speed
Results_1.cdf = gen_cdf(Speed_Exploration);

% Generate CDF for escape flight speed (threat)
Results_2.cdf = gen_cdf(Speed_Escape_Flight_Threat);

% Create a new figure for the CDF plots
figure;
set(gcf, 'color', 'w', 'position', [200 200 300 300]);  % Set the figure's background color and size

% Plot the CDF for exploration speed in black
plot(Results_1.cdf.x, Results_1.cdf.y, 'k'); 
hold on;  % Keep the figure open to add the next plot

% Plot the CDF for escape flight speed (threat) in red
plot(Results_2.cdf.x, Results_2.cdf.y, 'r');

% Customize axes properties
set(gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 8, 'FontName', 'Helvetica', ...
         'xtick', 0:10:100, 'ytick', 0:0.2:1);  % Customize ticks and fonts

title('Speed Distribution CDF');  % Add a title

%% Statistical analysis of the speed distributions
% Calculate and display the mean and standard error of the mean (SEM) for both conditions

exploration_mean = nanmean(Speed_Exploration);  % Mean speed during exploration
exploration_sem = nansem(Speed_Exploration);    % SEM during exploration

threat_mean = nanmean(Speed_Escape_Flight_Threat);  % Mean speed during escape flight (threat)
threat_sem = nansem(Speed_Escape_Flight_Threat);    % SEM during escape flight (threat)

disp(['Exploration Mean: ', num2str(exploration_mean), ' � ', num2str(exploration_sem)]);
disp(['Escape Flight Threat Mean: ', num2str(threat_mean), ' � ', num2str(threat_sem)]);

