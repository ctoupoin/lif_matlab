
clearvars
close all

if ~exist('paths','var')
    run('main_double.m')
end


%% User inputs

% Run to process
exp_date = 20210426;
run_num = 1;
convection = 'off';
img_to_display = 50;


%% Paths \& data reading

% Paths
[paths,spreadsheet_data,imgs_name1,imgs_name2] = read_data_paths(paths,convection,exp_date,run_num);

% Spreadsheet
fps = spreadsheet_data.fps;
misc = spreadsheet_data.misc;

% Miscellaneaous
load(paths.tform_1to2);
load(paths.tform_2to1);
load(paths.tform_1to2_d);
load(paths.tform_2to1_d);
load(paths.intensity_correction);
load(paths.temperature_calibration);
load(paths.params_cam1);
load(paths.params_cam2);


%% Read images

img1 = imrotate(double(imread(imgs_name1(img_to_display,:))),-90);
% img1_ud = img1;

img2 = imrotate(double(imread(imgs_name2(img_to_display,:))),90);
% img2_ud = img2;

img1_ud = rect(img1,eye(3),...
    params_cam1.fc,params_cam1.cc,params_cam1.kc,params_cam1.KK_new);
img2_ud = rect(img2,eye(3),...
    params_cam2.fc,params_cam2.cc,params_cam2.kc,params_cam2.KK_new);


%% Acquire coordinates

user_input = 0;
while user_input == 0
    figure,
    imshow(img1_ud,[0 40000])
    disp('Click on the upper left corner, upper right corner, bottom right corner and bottom left corner');
    [x,y] = ginput(4);
    
    mask_common_roi = round([min(y(1:2)) min(y(3:4)) min(x([1 4])) max(x(2:3))]);
    close

    figure,
    imshow(img1_ud(mask_common_roi(1):mask_common_roi(2),mask_common_roi(3):mask_common_roi(4)),[0 40000])
    
    user_input = input('Is result acceptable? [Y=1,N=else]');
end

save(paths.mask_common_roi,'mask_common_roi');