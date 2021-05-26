function temperature_data = read_temprature_data(db_filename,db_start_acquisition,imgs_to_process,dt,T_int)

fid = fopen(db_filename);
temperature_data_raw = textscan(fid,'%f%s%f','delimiter','\t');
fclose(fid);

field_names = unique(temperature_data_raw{1,2});
available_timesteps = unique(temperature_data_raw{1,1});

temperature_data = struct();
count_loop = 0;
for ii = 1:length(available_timesteps)
    count_loop = count_loop + 1;
    
    ind = find(temperature_data_raw{1,1} == available_timesteps(ii));
    
    for jj = 1:length(field_names)
        temperature_data.(field_names{jj})(count_loop,1) = temperature_data_raw{1,3}(ind(jj));
    end
    temperature_data.t(count_loop) = available_timesteps(ii) - available_timesteps(1);
end
temperature_data.t_h = temperature_data.t ./ 3600;


dt_temp = mean(diff(temperature_data.t));

acq_time = imgs_to_process .* dt;

[~,ind_acq_start] =  min(abs([temperature_data.t_h] - db_start_acquisition));

ind_acq_max = min(ind_acq_start + (round(acq_time ./ dt_temp)+1) + round(120 ./ dt_temp),length([temperature_data.t_h]));
ind_acq = ind_acq_start:ind_acq_max;
ind_temp =  (abs([temperature_data.T_internal] - T_int) < 0.1)';

ind_acq_ext = false(size(ind_temp));
ind_acq_ext(max(1,ind_acq_start - round(600 ./ dt_temp)):min(length(ind_temp),ind_acq_start + round(600 ./ dt_temp))) = true;

temperature_data_temp = temperature_data;

temperature_data.T_bas1 = temperature_data_temp.T_bain;
temperature_data.T_bas2 = temperature_data_temp.T_bas1;
temperature_data.T_bas3 = temperature_data_temp.T_bas2;
temperature_data.T_haut1 = temperature_data_temp.T_bas3;
temperature_data.T_haut2 = temperature_data_temp.T_haut1;
temperature_data.T_bain = temperature_data_temp.T_haut2;

T_down(1,:) = temperature_data.T_bas1;
T_down(2,:) = temperature_data.T_bas2;
T_down(3,:) = temperature_data.T_bas3;

T_up(1,:) = temperature_data.T_haut1;
T_up(2,:) = temperature_data.T_haut2;

temperature_data.T_down = mean(T_down,1);
temperature_data.T_up = mean(T_up,1);

temperature_data.T_down = mean(temperature_data.T_down(ind_acq));
temperature_data.T_up = mean(temperature_data.T_up(ind_acq));
temperature_data.Delta_T = temperature_data.T_down - temperature_data.T_up;

% [~,ind_T_int] = min(abs(ind_temp - ind_acq_start));
% abs(diff(ind_temp - ind_temp(ind_T_int)))

temperature_data.T_internal = mean(temperature_data.T_internal(ind_acq_ext & ind_temp));
temperature_data.ind_acq = ind_acq;