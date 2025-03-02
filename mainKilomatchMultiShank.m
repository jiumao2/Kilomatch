% set the path to Kilomatch and settings
path_kilomatch = '.\Kilomatch';
path_settings = '.\settings.json';

addpath(path_kilomatch);
addpath(genpath(fullfile(path_kilomatch, 'Functions')));

user_settings = jsonc.jsoncDecode(fileread(path_settings));
tic;

%% Kilomatch
Step1_Preprocess;

output_folder = user_settings.output_folder;
path_to_data = user_settings.path_to_data;

shanks_data = arrayfun(@(x)x.Kcoords(x.Channel), spikeInfo);
shankIDs = unique(spikeInfo(1).Kcoords);

for i_shank = 1:length(shankID)
    load(path_to_data);
    shankID = shankIDs(i_shank);
    spikeInfo = spikeInfo(shanks_data == shankID);

    user_settings.output_folder = fullfile(output_folder, ['Shank', num2str(shankID)]);

    Step2_MotionEstimation;
    Step3_ComputeWaveformFeatures;
    Step4_ComputeAllSimilarity;
    Step5_IterativeClustering;
    Step6_AutoCuration;
    Step7_VisualizeClusters;
end





