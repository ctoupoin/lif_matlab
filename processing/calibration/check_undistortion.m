%                          Check undistortion
%
%   This scripts applieds to un-distortion process to calibration target
%   images in order to check that it is correctly applied.
%
%--------------------------------------------------------------------------
%
%   Author : Cl\'ement Toupoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all

if ~exist('paths','var')
    run('../main_double.m')
end

format_imgs = 'tif';


%% Paths

addpath(genpath(paths.ext_functions));

paths.target = [paths.data 'calibration_target/'];
paths.params_cam0 = [paths.misc 'params_cam1.mat'];
paths.params_cam1 = [paths.misc 'params_cam2.mat'];

data_files = dir(paths.target);
count = 0;
for ii = 1:length(data_files)
    if length(data_files(ii).name) > 10 && strcmp(data_files(ii).name(1:11),'target_plan')
        count = count + 1;
        paths.target_plan(count,:) = [data_files(ii).folder '/' data_files(ii).name];
    end
end
paths.no_target = [data_files(ii).folder '/no_target/'];

load(paths.params_cam0);
load(paths.params_cam1);


%% Load  images

imgs = dir([paths.target_plan(1,:) '/' '*.' format_imgs]);

ind_cam0 = ismember(vertcat(imgs.name),'cam0');
ind_cam0 = ind_cam0(:,4) == true;
ind_cam0 = ind_cam0';
ind_cam1 = ~ind_cam0;

if sum(ind_cam0) == 0
    ind_cam0 = ismember(vertcat(imgs.name),'C1');
    ind_cam0 = ind_cam0(:,17) == true;
    ind_cam0 = ind_cam0';
    ind_cam1 = ~ind_cam0;
end

ind_cam0 = find(ind_cam0);
ind_cam1 = find(ind_cam1);

num_img = length(ind_cam0);

char_length = length([imgs(end).folder '/' imgs(end).name]);


img_name0 = [imgs(1).folder '/' imgs(1).name];
img0 = imrotate(double(imread(img_name0)),-90);

img_name1 = [imgs(1).folder '/' imgs(1).name];
img1 = imrotate(double(imread(img_name1)),90);


%% Undistortion

alpha_c = 0;

KK_new_0 = [params_cam0.fc(1) alpha_c*params_cam0.fc(1) params_cam0.cc(1);...
    0 params_cam0.fc(2) params_cam0.cc(2) ;...
    0 0 1];
KK_new_1 = [params_cam1.fc(1) alpha_c*params_cam1.fc(1) params_cam1.cc(1);...
    0 params_cam1.fc(2) params_cam1.cc(2) ;...
    0 0 1];


img0_ud1 = rect(img0,eye(3),params_cam0.fc,params_cam0.cc,params_cam0.kc,KK_new_0);
img1_ud1 = rect(img1,eye(3),params_cam1.fc,params_cam1.cc,params_cam1.kc,KK_new_1);

img0_ud = rect(img0,eye(3),params_cam0.fc,params_cam0.cc,params_cam0.kc,params_cam0.KK_new);
img1_ud = rect(img1,eye(3),params_cam1.fc,params_cam1.cc,params_cam1.kc,params_cam1.KK_new);


%% Plot

img0_adj = imadjust(img0);
img0_ud_adj = imadjust(img0_ud);