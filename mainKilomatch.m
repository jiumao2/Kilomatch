% set the path to Kilomatch and settings
path_kilomatch = '.\Kilomatch';
path_settings = '.\settings.json';

addpath(path_kilomatch);
addpath(genpath(fullfile(path_kilomatch, 'Functions')));

user_settings = jsonc.jsoncDecode(fileread(path_settings));
tic;

%% Kilomatch
Step1_Preprocess;
Step2_MotionEstimation;
Step3_ComputeWaveformFeatures;
Step4_ComputeAllSimilarity;
Step5_Clustering;
Step6_AutoCuration;
Step7_VisualizeClusters;





