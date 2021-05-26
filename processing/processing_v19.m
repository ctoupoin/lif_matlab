%                  *Processing preliminary experiments*
%
%   Script used to process the experiments of visualization of a thermal
%   plume using LIF technique.
%
%   Inputs:
%       - raw images
%       - camera calibration data and transfer matrix to go from the fov of
%      one camera to that of the other
%       - transfer matrix to max the histogram of one image with that of
%       the other one
%       - temperature calibration linking the mean temeprature of the cell
%       with light intensity
%
%   Output:
%       - temperature field inside the convection cell
%
%--------------------------------------------------------------------------
%
%   *Author:* Cl\'ement Toupoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changelog:
%   - v19 : simplified the script and removed now useless parts
%   - v18 : separated the computation of the Kalman filter variance to a
%   function
%   - v17 : ~~added the ability to process sevral runs with the same
%   experimental conditions~~ change not kept
%   - v16 : the Kalman filter is first implemented on a portion of the
%   images, before its noise value is kept for the rest
%   - v15 : ~~building of spatio-temporal diagrams~~ change not kept
%   - v14 : implemented a Kalman filter to remove noise
%   - v13 : changed the intensity correction to remove FFTs
%   - v12 : denoising tests
%   - v11 : changed camera names from 0/1 to 1/2
%   - v10 : adapted the code to work with a degree 2 polynomial fit for the
%   temperature calibration, and reversed the change of v9 so that it only
%   uses the transfer_function 1over0, which was shown to be better
%   - v9 : compares images obtained when using the 0over1 or the 1over0
%   transfer function
%   - v8 : modification of the ffts in the code to include 0 padding
%   - v7 : standardized geometric & intensity image transformations by
%   creating subfunctions
%   - v6 : allowed the script to run on runs with a uniform temperature
%   field
%   - v5 : changed the routine so that the geometric transform is done from
%   camera 1 to camera 0
%   - v4 :adapted the script to work with PCO cameras
%   - v3 : adapted the script for the use of camera calibration parameters
%   as computed by JYB's calibration toolbox, instead of MATLAB built-in
%   solution -which requires a dedicated toolbox-
%   - v2 : added call to main.m script where paths are dfined, and display
%   options are set

clearvars
close all

if ~exist('paths','var')
    run('main_double.m')
end


%% User inputs

% Run to process
exp_date = 20210517;
run_num = 1;

% Other user controls
convection = 'on';
disp_imgs = 'off';
save_data = 'on';
img_start = 5; % always start at 5, the first images are not well synchronized
img_end = 100;
imgs_to_process = (img_end - img_start) + 1;

% Display parameters
num_bins_colorscale = 500;
cmap = redblue(num_bins_colorscale);
T_scale_span = 2;
filter_size = 4;

% Filter parameters
num_imgs_filter = 500;
noise_var = 5e-2;
Kalman_filter_gain = 0.5;

% Physical parameters
g = 9.81;           % gravity
h = 41.5e-2;        % cell width
l = 10.5e-2;        % cell depht
q = 400 / (h * l);  % heat flux


%% Paths \& data reading

% Generate paths for the current run, and get name of the image files
[paths,spreadsheet_data,imgs_name1,imgs_name2] = read_data_paths(paths,convection,exp_date,run_num);

if strcmp(save_data,'on')
    mkdir(paths.output_processed_images_run)
end

% Get run info from the spreadsheet
fps = spreadsheet_data.fps;
dt = 1 / fps;
misc = spreadsheet_data.misc;
c_ratio = spreadsheet_data.c_rB / spreadsheet_data.c_r110;
db_filename = [paths.temperature_probes num2str(spreadsheet_data.db_name) '.tsv'];
db_start_acquisition = spreadsheet_data.db_time;

% Miscellaneaous
load(paths.mask_common_roi);
load(paths.tform_2to1);
load(paths.intensity_correction);
load(paths.temperature_calibration);
load(paths.params_cam1);
load(paths.params_cam2);
% load(paths.tform_1to2_d);
% load(paths.tform_2to1_d);
% load(paths.tform_1to2);

transfer_function_2over1 = intensity_correction.transfer_function_2over1;
gamma_px = params_cam1.gamma * 10^3;
% transfer_function_1over2 = intensity_correction.transfer_function_1over2;
% roi1_ref = intensity_correction.roi1;
% roi2_ref = intensity_correction.roi2;

disp([num2str(exp_date) ' run' num2str(run_num,'%02.f') ', convection ' convection])

if abs(temp_calib.c_ratio-c_ratio) > 0.1
    error('The concentration ratios differ between the current run and the temperature calibration');
end


%% Reading temperature

% Read temperature from the temperature probes. Accessing the database
% directly necessitates a toolbox in MATLAB, so the data needs to be
% extracted by the user in .tsv
T_int = spreadsheet_data.T_int; 
temperature_data = read_temprature_data(db_filename,db_start_acquisition,imgs_to_process,dt,T_int);

T_int = temperature_data.T_internal;
delta_T = temperature_data.Delta_T;
T_disp = T_int;


%% Physical characteristics of water \& dimensionless number

water_props = compute_water_properties_v2(T_int);

Ra = (g * water_props.thermal_expansion * temperature_data.Delta_T * h^3) /...
    (water_props.kinematic_viscosity * water_props.thermal_diffusivity);
Pr = water_props.kinematic_viscosity / water_props.thermal_diffusivity;
Nu = q * h / (water_props.thermal_conductivity * temperature_data.Delta_T);


%% Background creation w/o laser

% Images of the cell without laser illumination are recorded and averaged.
% They will serve as a background later on.
[bkgrd_no_laser1, bkgrd_no_laser2] = background_creation_v2(paths,...
    params_cam1,params_cam2,mask_common_roi,tform_2to1);
[size_x,size_y] = size(bkgrd_no_laser1);
size_px = size(bkgrd_no_laser1);
size_cm = (size_px ./ gamma_px) * 100;


%% Kalman filtering - getting variance

% A Kalman filter is a temporal filter used to reduce noise in the image
% series. Its variance converges iteratively. We first run the filter on a
% sufficient number of images (500 is enough) to get an acceptable value of
% the variance. This value is subsequently used to process the whole run,
% so that there is no difference in processing between the beginning and
% the end ofthe run.
% The data is stored in a file so that this operation does not need to be
% repeated on subsequent processing.

% see for the orignal file, slightly modified in the current code:
% https://www.mathworks.com/matlabcentral/fileexchange/26334-kalman-filter-for-noisy-movies
% https://imagej.nih.gov/ij/plugins/kalman.html

if ~exist(paths.Kalman_filter_variance,'file')
    [Kalman_var_final,mean_var] = get_kalman_variance(num_imgs_filter,noise_var,...
        Kalman_filter_gain,imgs_name1,imgs_name2,params_cam1,params_cam2,mask_common_roi,tform_2to1,...
        transfer_function_2over1,bkgrd_no_laser1,bkgrd_no_laser2,filter_size,img_start);
    
    Kalman_filter.variance = Kalman_var_final;
    Kalman_filter.gain = Kalman_filter_gain;
    Kalman_filter.mean_var = mean_var;
    Kalman_filter.num_imgs_filter = num_imgs_filter;
    Kalman_filter.noise_var = noise_var;
    
    save(paths.Kalman_filter_variance,'Kalman_filter')
    
    if strcmp(disp_imgs,'on')
        figure('Visible',disp_imgs,'position',[100 100 1000 1000]);
        plot(mean_var)
        xlabel('frames')
        ylabel('$\sigma$')
    end
else
    load(paths.Kalman_filter_variance)
    
    disp(['Variance of the Kalman filter loaded from data with ' num2str(Kalman_filter.num_imgs_filter) ' images used.'])
    
    if strcmp(disp_imgs,'on')
        figure('Visible',disp_imgs,'position',[100 100 1000 1000]);
        plot(processed_data.mean_var)
        xlabel('frames')
        ylabel('$\sigma$')
    end
end


%% Image processing + Kalman filtering

% Main processing step, transforming raw images into temperature field.
% The temperature fields are saved as .mat if save_imgs is 'on'.
[processing_time,mean_T,sigma_T] = processing_filtering_v2(imgs_name1,imgs_name2,...
    paths,params_cam1,params_cam2,mask_common_roi,tform_2to1,transfer_function_2over1,...
    bkgrd_no_laser1,bkgrd_no_laser2,filter_size,img_start,img_end,size_px,...
    save_data,Kalman_filter,temp_calib);

disp([num2str(imgs_to_process) ' images procesed in ' num2str(processing_time.min)...
    'min ' num2str(processing_time.sec,'%.0f') 'sec'])