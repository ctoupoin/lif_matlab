function [processing_time,mean_T,sigma_T] = processing_filtering_v2(imgs_name1,imgs_name2,...
    paths,params_cam1,params_cam2,mask_common_roi,tform_2to1,transfer_function_2over1,...
    bkgrd_no_laser1,bkgrd_no_laser2,filter_size,img_start,img_end,size_px,...
    save_data,Kalman_filter,temp_calib)

% Changelog:
%   - v2: removed the option to make a video

imgs_to_process = (img_end - img_start) + 1;

% Keep track of the mean temperature, and of the std over a horizontal
% profile. Useful to evaluate the accuracy of the temperature over a
% uniform temperature field.
mean_T = zeros(1,imgs_to_process);
sigma_T = zeros(1,imgs_to_process);

tic;
% Initialization for the Kalman filter
im = image_division_v1(imgs_name1,imgs_name2,params_cam1,...
    params_cam2,mask_common_roi,tform_2to1,transfer_function_2over1,...
    bkgrd_no_laser1,bkgrd_no_laser2,filter_size,img_start);
predicted_im = im;
count_loop = 0;
for ii = img_start+1:img_end
    count_loop = count_loop + 1;
    
    im = image_division_v1(imgs_name1,imgs_name2,params_cam1,...
        params_cam2,mask_common_roi,tform_2to1,transfer_function_2over1,...
        bkgrd_no_laser1,bkgrd_no_laser2,filter_size,ii);
    
    observed_im = im;
    Kalman = Kalman_filter.variance ./ (Kalman_filter.variance + Kalman_filter.noise_var);
    corrected_im = Kalman_filter.gain * predicted_im + (1.0 - Kalman_filter.gain) *...
        observed_im + Kalman .* (observed_im - predicted_im);
    predicted_im = corrected_im;
    im = corrected_im;
    
    % Temperature conversion
    im_T = polyval(temp_calib.fit,im);
    im_T = flip(im_T,1);
    
    mean_T(count_loop) = mean(mean(im_T));
    sigma_T(count_loop) = std(im_T(round(size_px(1)/2),5:end-5));
    
    if strcmp(save_data,'on')
        if ~exist(paths.output_processed_images_run,'dir')
            mkdir(paths.output_processed_images_run)
        end
       image_name = [paths.output_processed_images_run 'im_T_' num2str(count_loop,'%04.f') '.mat'];
       save(image_name,'im_T') 
    end
end
run_time = toc;

processing_time.min = floor(run_time/60);
processing_time.sec = rem(run_time,60);
