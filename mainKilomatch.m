% Set the path to Kilomatch and settings
path_kilomatch = '.\Kilomatch'; % The path where Kilomatch is installed
path_settings = '.\settings.json'; % Please make sure the settings in the file are accurate

addpath(path_kilomatch);
addpath(genpath(fullfile(path_kilomatch, 'Functions')));

user_settings = jsonc.jsoncDecode(fileread(path_settings)); % Read the settings
tic;

%% Run Kilomatch
Step1_Preprocess;
Step2_MotionEstimation;
Step3_ComputeWaveformFeatures;
Step4_ComputeAllSimilarity;
Step5_IterativeClustering;
Step6_AutoCuration;
Step7_VisualizeClusters;





