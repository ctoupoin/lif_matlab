function [paths,spreadsheet_data,imgs_name1,imgs_name2] = read_data_paths(paths,convection,exp_date,run_num)

% if nargin < 5
%     run_num_end = NaN;
% end

%% Paths

if strcmp(convection,'on')
    paths.expdate = [paths.conv_on num2str(exp_date) '/'];
    paths.output_processed_images_run = [paths.output_processed_images num2str(exp_date)...
        '/run' num2str(run_num,'%02.f') '/'];
    paths.movie = [paths.output_processed_images_run num2str(exp_date) '_run' num2str(run_num,'%02.f') '.avi'];
    filename_spreadsheet = paths.spreadsheet;
else
    paths.expdate = [paths.conv_off num2str(exp_date) '/'];
    paths.output_processed_images_run = [paths.output_processed_images 'conv_off/' num2str(exp_date)...
        '/run' num2str(run_num,'%02.f') '/'];
    paths.movie = [paths.output_processed_images_run 'conv_off/' num2str(exp_date) '_run' num2str(run_num,'%02.f') '.avi'];
    filename_spreadsheet = paths.spreadsheet_conv_off;
end

paths.run = [paths.expdate 'run' num2str(run_num,'%02.f') '/'];
paths.output_file = [paths.output_data num2str(exp_date) '_run' num2str(run_num,'%02.f') '.mat'];
paths.bkgrd_no_laser = [paths.conv_on num2str(exp_date) '/background_no_laser/'];
paths.bkgrd_no_laser_averaged_cam1 = [paths.bkgrd_no_laser 'avgd_cam1.tif'];
paths.bkgrd_no_laser_averaged_cam2 = [paths.bkgrd_no_laser 'avgd_cam2.tif'];
paths.db = [paths.run 'session.db'];
paths.db_float = [paths.run 'db_float.tsv'];


%% Spreadsheet reading

disp(['Reading ' filename_spreadsheet '.'])
fid = fopen(filename_spreadsheet);
spreadsheet = textscan(fid,repmat('%s',1,17),'delimiter','\t');
fclose(fid);

ind_date = strcmp({spreadsheet{1,1}{:}},num2str(exp_date));
ind_run = strcmp({spreadsheet{1,2}{:}},num2str(run_num,'%02.f'));

spreadsheet_data.num_imgs_th = str2double(spreadsheet{3}{ind_date & ind_run});
spreadsheet_data.laser_pwr = str2double(spreadsheet{4}{ind_date & ind_run});
spreadsheet_data.T_bulk = str2double(spreadsheet{5}{ind_date & ind_run});
spreadsheet_data.T_int = str2double(spreadsheet{6}{ind_date & ind_run});
spreadsheet_data.T_up = str2double(spreadsheet{7}{ind_date & ind_run});
spreadsheet_data.T_down = str2double(spreadsheet{8}{ind_date & ind_run});
spreadsheet_data.T_room = str2double(spreadsheet{9}{ind_date & ind_run});
spreadsheet_data.c_rB = str2double(spreadsheet{10}{ind_date & ind_run});
spreadsheet_data.c_r110 = str2double(spreadsheet{11}{ind_date & ind_run});
spreadsheet_data.exposure_cam1 = str2double(spreadsheet{12}{ind_date & ind_run});
spreadsheet_data.exposure_cam2 = str2double(spreadsheet{13}{ind_date & ind_run});
spreadsheet_data.fps = str2double(spreadsheet{14}{ind_date & ind_run});
spreadsheet_data.db_name = str2double(spreadsheet{15}{ind_date & ind_run});
spreadsheet_data.db_time = str2double(spreadsheet{16}{ind_date & ind_run});
spreadsheet_data.misc = spreadsheet{17}{ind_date & ind_run};


%% Gather images name

% if isnan(run_num_end)
imgs = dir([paths.run]);
format = imgs(3).name(end-2:end);
imgs = dir([paths.run '*.' format]);
imgs_names = vertcat(imgs.name);

ind_cam1 = imgs_names(:,17) == '1';
ind_cam1 = ind_cam1';
ind_cam2 = ~ind_cam1;

ind_cam1 = find(ind_cam1);
ind_cam2 = find(ind_cam2);

lgth_filename = length([imgs(end).folder '/' imgs(end).name]);
imgs_name1 = repmat('a',length(ind_cam1),lgth_filename);
imgs_name2 = repmat('a',length(ind_cam2),lgth_filename);
for ii = ind_cam1
    imgs_name1(ii,:) = [imgs(ii).folder '/' imgs(ii).name];
end
count = 0;
for ii = ind_cam2
    count = count + 1;
    imgs_name2(count,:) = [imgs(ii).folder '/' imgs(ii).name];
end
% else
%     num_imgs = zeros(1,run_num_end-run_num+1);
%
%     count_loop = 0;
%     for ii = run_num:run_num_end
%         count_loop = count_loop + 1;
%
%         path_run = [paths.expdate 'run' num2str(ii,'%02.f') '/'];
%
%         imgs = dir(path_run);
%         format = imgs(3).name(end-2:end);
%         imgs = dir([paths.run '*.' format]);
%         imgs_names = vertcat(imgs.name);
%
%         ind_cam1 = imgs_names(:,17) == '1';
%         ind_cam1 = find(ind_cam1);
%
%         num_imgs(count_loop) = length(ind_cam1);
%         lgth_filename = length([imgs(end).folder '/' imgs(end).name]);
%     end
%     num_imgs = sum(num_imgs);
%
%     imgs_name1 = repmat('a',num_imgs,lgth_filename);
%     imgs_name2 = repmat('a',num_imgs,lgth_filename);
%
%     count_loop = 0;
%     for ii = run_num:run_num_end
%         count_loop = count_loop + 1;
%
%         path_run = [paths.expdate 'run' num2str(ii,'%02.f') '/'];
%
%         imgs = dir(path_run);
%         format = imgs(3).name(end-2:end);
%         imgs = dir([paths.run '*.' format]);
%         imgs_names = vertcat(imgs.name);
%
%         ind_cam1 = imgs_names(:,17) == '1';
%         ind_cam1 = ind_cam1';
%         ind_cam2 = ~ind_cam1;
%
%         ind_cam1 = find(ind_cam1);
%         ind_cam2 = find(ind_cam2);
%
%         for jj = ind_cam1
%             imgs_name1(count_loop,:) = [imgs(jj).folder '/' imgs(jj).name];
%         end
% %         count = 0;
%         for jj = ind_cam2
% %             count = count + 1;
%             imgs_name2(count_loop,:) = [imgs(jj).folder '/' imgs(jj).name];
%         end
%     end
% end


end