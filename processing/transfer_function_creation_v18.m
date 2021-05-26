%                      *Transfer function creation*
%
%   Script used to create an otpical trasnfer function, which matches the
%   histogram of the images of one camera so that it matches that of the
%   other camera.
%
%--------------------------------------------------------------------------
%
%   *Author :* Cl\'ement Toupoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changelog:
%   - v18: attempt to allow the transfer function to evolve with time
%   - v17 : removed FFTs altogether
%   - v16 : opted to compute the trasnfer function based on a median image,
%   instead of on an average image
%   - v15 : chaned the routine so that only one transfer function is
%   created, instead of the previous two, now that the most effective one is
%   known (2over1)
%   - v14 : extended MATLAB's regular padding for the FFTs
%   - v13 : went back to regular MATLAB padding for the FFTs
%   - v12c : cleaned up code
%   - v12 : changed the 0 padding
%   - v11 : averaged image first, THEN create transfer function, isntead
%   of averaging transfer functions
%   - v10 : modification of the ffts in the code, with 0 padding and perhaps
%   Hamming windowing
%   - v9 : standardized geometric & intensity image transformations by
%   creating subfunctions
%   - v8 : change the script so that it computes the transfer function
%   based on a single run, specified by the user with date and run number
%   - v7 : basically reverses change done in v4, so that we create transfer
%   function for each pair of images, then average them
%   - v6 : changed the routine so that the geometric transform is done from
%   camera 1 to camera 0
%   - v5 : adapted the routine for PCO cameras, along with various QOL
%   changes
%   - v4 : changed to script so that it averages images, THEN creates the
%   transfer function instead of averaging transfer functions
%   - v3 : Created transfer functions for several concentration ratios and
%   laser powers
%   - v2 : implemented new background creation function

clearvars
close all

if ~exist('paths','var')
    run('main_double.m')
end


%% User inputs

exp_date = 20210524;
run_num = 1;
img_start = 5;
img_end = 2600;
imgs_to_average = (img_end - img_start) + 1;


%% Data reading
% Load previously-computed data -geometric transformation for image
% merging, run data and so on...-

[paths,spreadsheet_data,imgs_name1,imgs_name2] = read_data_paths(paths,'off',exp_date,run_num);

c_ratio = spreadsheet_data.c_rB / spreadsheet_data.c_r110;
T_int = spreadsheet_data.T_int;

load(paths.mask_common_roi);
load(paths.tform_1to2);
load(paths.tform_2to1);
load(paths.tform_1to2_d);
load(paths.tform_2to1_d);
load(paths.params_cam1);
load(paths.params_cam2);


%% Reading database

disp([num2str(exp_date) ' run' num2str(run_num,'%02.f') ', T=' num2str(T_int) '°C, P='...
    num2str(spreadsheet_data.laser_pwr) 'W, c_B/c_110=' num2str(c_ratio)]);


%% Background creation
% Images of the cell without laser illumination are recorded and averaged.
% They will serve as a background later on.

[bkgrd_no_laser1, bkgrd_no_laser2] = background_creation_v2(paths,...
            params_cam1,params_cam2,mask_common_roi,tform_2to1);

        
%% Transfer function creation
% Create the transfer function defined as $G = I^{110}_u ./ I^{B}_u,
% where $I^{110}_u \& I^{B}_u$ are images of the fluorescence of
% rhodamine 110 and B respectively, with a uniform temperature field.

[size_x,size_y] = size(bkgrd_no_laser1);

% Average images in the chosen run
num_img = min(imgs_to_average,size(imgs_name1,1));
I1 = zeros(1,num_img);
I2 = zeros(1,num_img);
% roi_stack1 = zeros(size_x,size_y,imgs_to_average);
% roi_stack2 = zeros(size_x,size_y,imgs_to_average);
count_loop = 0;
for jj = img_start:img_end
    count_loop = count_loop + 1;
    
    [roi1,roi2] = geometric_transforms(imgs_name1,imgs_name2,params_cam1,...
        params_cam2,mask_common_roi,tform_2to1,jj);
    
    if any((size(bkgrd_no_laser1) == size(roi1)) == false)
        disp('Background size incorrect')
        
        delete(paths.bkgrd_no_laser_averaged_cam1);
        delete(paths.bkgrd_no_laser_averaged_cam2);
        
        [bkgrd_no_laser1, bkgrd_no_laser2] = background_creation_v2(paths,...
            params_cam1,params_cam2,mask_common_roi,tform_2to1);
    end
    
    roi1 = imsubtract(roi1,bkgrd_no_laser1);
    roi2 = imsubtract(roi2,bkgrd_no_laser2);
    
%     roi_stack1(:,:,count_loop) = roi1;
%     roi_stack2(:,:,count_loop) = roi2;
    
    I1(count_loop) = mean(mean(roi1));
    I2(count_loop) = mean(mean(roi2));
end

% roi1_median = median(roi_stack1,3);
% roi2_median = median(roi_stack2,3);
%         
% % Computes the trasnfer function
% transfer_function_2over1 = roi2_median ./ roi1_median;
% transfer_function_1over2 = roi1_median ./ roi2_median;


%% Saving data

intensity_correction.I1 = I1;
intensity_correction.I2 = I2;

intensity_correction.transfer_function_2over1 = transfer_function_2over1;
intensity_correction.transfer_function_1over2 = transfer_function_1over2;
intensity_correction.roi1 = roi1_median;
intensity_correction.roi2 = roi2_median;
intensity_correction.temperature = T_int;
intensity_correction.laser_pwr = spreadsheet_data.laser_pwr;
intensity_correction.c_ratio = c_ratio;
intensity_correction.optical_params = [spreadsheet_data.exposure_cam1 spreadsheet_data.exposure_cam2];
intensity_correction.exp_date = exp_date;
intensity_correction.run_num = run_num;
intensity_correction.num_imgs = imgs_to_average;

% save(paths.intensity_correction,'intensity_correction','-v7.3')