%                       *Temperature calibration*
%
%   Temperature calibration for the LIF image processing. Images of the
%   cell with an uniforme temperature field and the laser sheet are taken.
%   The whole processing used in processing_vX.m is applied to these
%   images. The mean intensity value of the resulting image is plotted
%   versus the temperature of the cell.
%
%   Input : 
%       - runs with uniform temperature
%
%   Output :
%       - polynomial fit of the light intensity versus temperature
%
%--------------------------------------------------------------------------
%
%   *Author :* Cl\'ement Toupoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changelog:
%   - v11 : used the median image instead of the mean one
%   - v10 : clean-up comments
%   - v9 : converts camera notation from 0/1 to 1/2
%   - v8 : reverses change made in v7 so now only one temperature
%   calibration is computed
%   - v7 : plots 2 temperature calibrations, depending on which transfer
%   function is used
%   - v6b : clean-up of the code
%   - v6 : standardized geometric & intensity image transformations by
%   creating subfunctions
%   - v5 : adapted the routine to work with the new transfer function
%   syntax
%   - v4 : changed the routine so that the geometric transform is done from
%   camera 1 to camera 0
%   - v3 : adpated code for using PCO cameras
%   - v2 : added comparisons between several concentration ratios

clearvars
close all

if ~exist('paths','var')
    run('main_double.m')
end


%% User inputs

% Desired experimental parameters
laser_pwr = 4;
desired_c_ratio = 6;
imgs_to_average = 30;
imgs_reject = 3;
deg_poly = 1;
filter_size = 4;


%% Reading data

load(paths.mask_common_roi);

load(paths.tform_2to1);
load(paths.intensity_correction);
load(paths.params_cam1);
load(paths.params_cam2);
% load(paths.tform_1to2);
% load(paths.tform_1to2_d);
% load(paths.tform_2to1_d);

transfer_function_2over1 = intensity_correction.transfer_function_2over1;
% transfer_function_1over2 = intensity_correction.transfer_function_1over2;
% roi1_ref = intensity_correction.roi1;
% roi2_ref = intensity_correction.roi2;

fid = fopen(paths.spreadsheet_conv_off);
spreadsheet = textscan(fid,repmat('%s',1,17),'delimiter','\t');
fclose(fid);

ind_laser_pwr = str2double(spreadsheet{1,4}(:)) == laser_pwr;
if isstring(desired_c_ratio)
    ind_c_ratio = true(1,length(spreadsheet{1,1}));
else
    ind_c_ratio = abs( (str2double(spreadsheet{1,10}(:)) ./...
        str2double(spreadsheet{1,11}(:))) - desired_c_ratio) < 0.1;
end

ind_OK = str2double(spreadsheet{1,17}(:));

ind_calib_temp = find(ind_laser_pwr & ind_c_ratio & ind_OK)';

available_temperatures = str2double(spreadsheet{1,5}(ind_calib_temp));
[~,ind_temp_min] = min(available_temperatures);
[~,ind_temp_max] = max(available_temperatures);
ind_temp_min = ind_calib_temp(ind_temp_min);
ind_temp_max = ind_calib_temp(ind_temp_max);


%%  Compute light intensity

% Initialization
T_bulk = zeros(1,length(ind_calib_temp));
I1 = zeros(1,length(ind_calib_temp));
I2 = zeros(1,length(ind_calib_temp));
I1_final = zeros(1,length(ind_calib_temp));
I2_final = zeros(1,length(ind_calib_temp));
I = zeros(1,length(ind_calib_temp));
T_int = zeros(1,length(ind_calib_temp));
c_ratio = zeros(1,length(ind_calib_temp));
laser_pwr = zeros(1,length(ind_calib_temp));
int_gamma = zeros(1,length(ind_calib_temp));
exp_date_list = repmat('a',length(ind_calib_temp),8);

count_loop = 0;
for ii = ind_calib_temp
    count_loop = count_loop + 1;
    
    exp_date = spreadsheet{1,1}{ii};
    run_num = spreadsheet{1,2}{ii};
    
    [paths,~,imgs_name1,imgs_name2] = read_data_paths(paths,'off',exp_date,run_num);
    
    laser_pwr(count_loop) = str2double(spreadsheet{1,4}(ii));
    c_ratio(count_loop) = str2double(spreadsheet{1,10}(ii)) ./ str2double(spreadsheet{1,11}(ii));
    exp_date_list(count_loop,:) = exp_date;
    
    db_filename = [paths.temperature_probes spreadsheet{1,15}{ii,1} '.tsv'];
    db_start_acquisition = str2double(spreadsheet{1,16}{ii});

    fps = str2double(spreadsheet{1,15}{ii,1});
    dt = 1 / fps;
    T_int(count_loop) = str2double(spreadsheet{1,6}(ii));
    
    temperature_data = read_temprature_data(db_filename,db_start_acquisition,imgs_to_average,dt,T_int(count_loop));
    T_int(count_loop) = temperature_data.T_internal;
    
    
    % Background creation
    % Images of the cell without laser illumination are recorded and averaged.
    % They will serve as a background later on.
    [bkgrd_no_laser1, bkgrd_no_laser2] = background_creation_v2(paths,...
        params_cam1,params_cam2,mask_common_roi,tform_2to1_d);
    [size_x,size_y] = size(bkgrd_no_laser1);
    
    
    % Intensity computation
    num_img = min(imgs_to_average,size(imgs_name1,1))+1;
    
    int1 = zeros(1,num_img);
    int2 = zeros(1,num_img);
    int1_final = zeros(1,num_img);
    int2_final = zeros(1,num_img);
    int_gamma = zeros(1,num_img);
    im_stack = zeros(size_x,size_y,num_img);
    count_loop2 = 0;
    for jj = imgs_reject+1:num_img+imgs_reject
        count_loop2 = count_loop2 + 1;
        
        [roi1,roi2] = geometric_transforms(imgs_name1,imgs_name2,params_cam1,params_cam2,mask_common_roi,...
            tform_2to1_d,jj);
        
        if any((size(bkgrd_no_laser1) == size(roi1)) == false)
            disp('Background size incorrect')
            
            delete(paths.bkgrd_no_laser_averaged_cam1);
            delete(paths.bkgrd_no_laser_averaged_cam2);
            
            [bkgrd_no_laser1, bkgrd_no_laser2] = background_creation_v2(paths,...
                params_cam1,params_cam2,mask_common_roi,tform_2to1_d);
            [size_x,size_y] = size(bkgrd_no_laser1);
        end
        
        roi1 = imsubtract(roi1,bkgrd_no_laser1);
        roi2 = imsubtract(roi2,bkgrd_no_laser2);

        roi1_c = roi1 .* transfer_function_2over1;
%         roi2_c = roi2 .* transfer_function_1over2;
        
        im = imdivide(roi1_c,roi2);
%         im = imdivide(roi1,roi2_c);
        im = wiener2(im,[4 4]);

        im_stack(:,:,count_loop2) = im;
        
        int1(count_loop2) = mean(mean(roi1));
        int2(count_loop2) = mean(mean(roi2));
        int_gamma(count_loop2) = mean(mean(im));
    end
    
    im_med = median(im_stack,3);

    I1(count_loop) = mean(int1);
    I2(count_loop) = mean(int2);

    I(count_loop) = mean(mean(im_med));
end


%% Save data

temp_calib.c_ratio = desired_c_ratio;

[T_sorted, ind_sort] = sort(T_int);
I1_sorted = I1(ind_sort);
I2_sorted = I2(ind_sort);
I_sorted = I(ind_sort);
exp_date_sorted = str2num(exp_date_list(ind_sort,:));
ind_date = exp_date_sorted > 20210519;

if length(unique(T_sorted)) > 1
    temp_calib.fit = polyfit(I_sorted,T_sorted,deg_poly);
    temp_calib.fit1 = polyfit(I1_sorted,T_sorted,1);
    temp_calib.fit2 = polyfit(I2_sorted,T_sorted,1);
else
    temp_calib.fit = NaN;
    temp_calib.fit1 = NaN;
    temp_calib.fit2 = NaN;
end

temp_calib.transfer_function_date = intensity_correction.exp_date;
temp_calib.transfer_function_run_num = intensity_correction.run_num;
temp_calib.exp_date = exp_date_sorted;
temp_calib.T_int = T_sorted;
temp_calib.I1 = I1_sorted;
temp_calib.I2 = I2_sorted;
temp_calib.I = I_sorted;
temp_calib.t_min = min(available_temperatures);
temp_calib.t_max = max(available_temperatures);

save(paths.temperature_calibration,'temp_calib');


%%

ind_date = exp_date_sorted > 20210519;

figure('position',[0 0 1000 1000]),
hold on,

plot(T_sorted,I_sorted,marker_order{1},...
    'Color',col_order(1,:),'LineWidth',2,'MarkerSize',15);
plot(T_sorted(ind_date),I_sorted(ind_date),marker_order{1},...
    'Color',col_order(2,:),'LineWidth',2,'MarkerSize',15);

xlim([min(T_sorted)-5 max(T_sorted)+5])
ylim([min(I_sorted)-0.2 max(I_sorted)+0.2])

y_fit = linspace(0.5,1.5,100);
x_fit = polyval(temp_calib.fit,y_fit);
plot(x_fit,y_fit,'Color',col_order(1,:),'LineWidth',2,'MarkerSize',15);
xlabel('T ($^{\circ}C$)')
ylabel('$\langle I_U \rangle$')

export_fig([paths.output 'temperature_calibration_date.pdf'],'-transparent','-pdf')


%%

figure('Position',[0 0 1000 1000]),
hold on,
plot(T_sorted,I2_sorted,marker_order{2},...
    'Color',col_order(2,:),'LineWidth',2,'MarkerSize',15);
plot(T_sorted,I1_sorted,marker_order{1},...
    'Color',col_order(1,:),'LineWidth',2,'MarkerSize',15);
yl = ylim;
xl = xlim;
y_fit = linspace(min(yl),max(yl),100);
% x_fit = polyval(temp_calib.fit1,y_fit);
% x_fit2 = polyval(temp_calib.fit2,y_fit);

% plot(x_fit,y_fit,'Color',col_order(1,:),'LineWidth',2,'MarkerSize',15);
% plot(x_fit2,y_fit,'Color',col_order(2,:),'LineWidth',2,'MarkerSize',15);
axis([min(xl) max(xl) 0 max(yl)])
legend('Rh. 110','Rh. B','location','best')
xlabel('T [$^{\circ}$C]')
ylabel('$I_{B, 110}$')

% export_fig([paths.output 'light_intensity_rhodamines_T.pdf'],'-transparent','-pdf')


%%
% figure('Position',[0 0 1000 1000]),
% hold on,
% plot(T_sorted,I1_sorted./I2_sorted,...
%     marker_order{1},'Color',col_order(1,:),'LineWidth',2,'MarkerSize',15);
% % ylim([0.9 1.4])
% xlabel('T [$^{\circ}$C]')
% ylabel('$I_{B} / I_{110}$')

% export_fig([paths.output 'light_intensity_ratio_T.pdf'],'-transparent','-pdf')


