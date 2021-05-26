%                                *Main*
%
%   Sets up paths and display options used in the  other scripts.
%
%--------------------------------------------------------------------------
%
%   *Author :* Cl\'ement Toupoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Paths

paths.root = '../';
paths.ext_functions = [paths.root '/ext_functions'];

paths.data = [paths.root 'data/'];
paths.output = [paths.root 'output/'];

paths.conv_on = [paths.data 'conv_on/'];
paths.conv_off = [paths.data 'conv_off/'];
paths.misc = [paths.data 'misc/'];
paths.dyes = [paths.data 'dyes/'];
paths.temperature_probes = [paths.data 'temperature_probes/'];

paths.spreadsheet = [paths.data 'spreadsheet_conv_on.tsv'];
paths.spreadsheet_conv_off = [paths.data 'spreadsheet_conv_off.tsv'];

paths.output_processed_images = [paths.output 'temperature_fields/'];
paths.output_data = [paths.output 'data/'];

paths.mask_common_roi = [paths.misc 'mask_common_roi.mat'];
paths.Kalman_filter_variance = [paths.misc 'Kalman_filter_variance.mat'];
paths.tform_1to2 = [paths.misc 'tform_1to2.mat'];
paths.tform_2to1 = [paths.misc 'tform_2to1.mat'];
paths.tform_1to2_d = [paths.misc 'tform_1to2_d.mat'];
paths.tform_2to1_d = [paths.misc 'tform_2to1_d.mat'];
paths.stereo_params = [paths.misc 'stereo_params.mat'];
paths.params_cam1 = [paths.misc 'params_cam1.mat'];
paths.params_cam2 = [paths.misc 'params_cam2.mat'];
paths.params_cam_stereo = [paths.misc 'params_cam_stereo.mat'];
paths.intensity_correction = [paths.misc 'intensity_correction.mat'];
paths.intensity_correction_unpadded = [paths.misc 'intensity_correction_unpadded.mat'];
paths.temperature_calibration = [paths.misc 'temperature_calibration.mat'];
paths.transfer_function_run_averaged = [paths.misc 'transfer_function_run_averaged.mat'];
paths.crenel_mask = [paths.misc 'crenel_mask.mat'];

paths.temp_dep_Sakakibara = [paths.misc 'Sakakibara1999_temp_dependence.csv'];


addpath(genpath(paths.ext_functions));


%% Display otpions

warning('off','MATLAB:legend:IgnoringExtraEntries');
set(0,'defaulttextinterpreter','latex')

col_order = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840  0.5560;...
    0.9451, 0.5490, 0.1333;...
    0.3010, 0.7450, 0.9330;...
    0.6350, 0.0780, 0.1840;...
    0.3255, 0.3176, 0.3294;...
    0.6706, 0.4078, 0.3412];

marker_order = {'o','s','d','^','<','>','v','p','h'};

% Get info for display purposes
ScrSize=get(0,'MonitorPositions');

% Setting default parameters for figures
set(groot,'DefaultAxesFontSize',24);
set(groot,'DefaultAxesposition',[0.1300 0.1100 0.7750 0.8150]);
set(groot,'DefaultAxesBox','on');
set(groot,'defaultAxesColorOrder',col_order);
set(groot,'DefaultTextInterpreter', 'latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextFontName','Times')
set(groot,'defaultfigurecolor',[1 1 1])

LWidth = 2;
MSize = 15;
