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
Step3_ComputePC_Features;
Step4_ComputeWaveformFeatures;
Step5_ComputeAllSimilarity;
Step6_Clustering;
Step7_AutoCuration;
Step8_VisualizeClusters;





