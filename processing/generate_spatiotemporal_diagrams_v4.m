%                   *Generate spatio-temporal diagrams*
%
%   Generate spatio-temporal diagrams from temperature fields
%
%   Inputs :
%       - processed temperature fields
%
%   Output :
%       - spatio-temporal diagrams
%
%--------------------------------------------------------------------------
%
%   *Author :* Cl\'ement Toupoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changelog:
%   - v4: added the optionn to choose at which imae to start processing,
%   for each run
%   - v3: added the ability to process multiple runs at once
%   - v2: changed functionnalities


clearvars
close all

if ~exist('paths','var')
    run('main_double.m')
end


%% User input

exp_date_to_process = [20210421, 20210507, 20210517, 20210519];
run_num_to_process = [1, 1, 1, 1];
convection = 'on';
img_start = [800, 800, 1, 1];
img_end = [2595, 2595, 2595, 2595];
% imgs_to_process = img_end - img_start;
disp_imgs = 'off';
save_imgs = 'off';

path_output = 'D:/Post_doc/LateX/slides_20210521/figures/';

sub_roi = [5 1524 6 1576];
filter_size = 4;

% Color display
num_bins_colorscale = 500;
cmap = redblue(num_bins_colorscale);
T_scale_span = 2;
T_rms_display = 0.5;

% Physical parameters
width_rugosity_cm = 0.5;
heigth_rugosity_cm = 0.2;

g = 9.81;
H = 40e-2;
Q = 400;


%% Load data for all runs

% Miscellaneaous
load(paths.params_cam1);
load(paths.params_cam2);
load(paths.temperature_calibration)

gamma_px = params_cam1.gamma * 10^3;

% Measurement positions
T_roi_size = 10e-2*gamma_px;

height_rugosity_px = round(heigth_rugosity_cm*10^(-2) * gamma_px);
half_rugosity = round(0.5*height_rugosity_px);

position_horizontal_profiles_cm = [0.5 1:39 39.5];
position_vertical_profiles_cm = [0.52 1.04:1.04:40.56 41.08];

position_BL_z_cm = [0.03 0.06 0.08 0.2 0.3 1];
position_notch_x_px = 783;
position_plot_x_px = 765;

position_single_point_x_px = [17 760 1498];
position_single_point_z_px = [114 780 1435];


%% Processing loop

% all_processed_data = struct();
count_main_loop = 0;
for ii = 1:length(exp_date_to_process)
    count_main_loop = count_main_loop + 1;
    
    exp_date = exp_date_to_process(ii);
    run_num = run_num_to_process(ii);
    
    processed_data.exp_date = exp_date;
    processed_data.run_num = run_num;
    
    imgs_to_process = img_end(ii)-img_start(ii)+1;
    
    
    %% Load data
    
    [paths,spreadsheet_data,imgs_name1,imgs_name2] = read_data_paths(paths,convection,exp_date,run_num);
    
    fps = spreadsheet_data.fps;
    misc = spreadsheet_data.misc;
    db_filename = [paths.temperature_probes num2str(spreadsheet_data.db_name) '.tsv'];
    db_start_acquisition = spreadsheet_data.db_time;
    
    imgs_processed = dir([paths.output_processed_images_run '*.mat']);

    dt = 1 / fps;
    c_ratio = spreadsheet_data.c_rB / spreadsheet_data.c_r110;

    
    %% Reading temperature
    
    T_int = spreadsheet_data.T_int;
    T_room = spreadsheet_data.T_room;
    temperature_data = read_temprature_data(db_filename,db_start_acquisition,imgs_to_process,dt,T_int);
    
    T_int = temperature_data.T_internal;
    delta_T = temperature_data.Delta_T;
    T_disp = T_int;
    
    water_props = compute_water_properties_v2(T_int);
    
    Ra = (g * water_props.thermal_expansion * temperature_data.Delta_T * (H)^3) /...
        (water_props.kinematic_viscosity * water_props.thermal_diffusivity);
    Nu = Q * H  / (water_props.thermal_conductivity * temperature_data.Delta_T);
    
    
    Ra_smooth = (2 * g * water_props.thermal_expansion * (temperature_data.T_internal - temperature_data.T_up) *...
        (40e-2)^3) / (water_props.kinematic_viscosity * water_props.thermal_diffusivity);
    Ra_rough = (2 * g * water_props.thermal_expansion * (temperature_data.T_down - temperature_data.T_internal) *...
        (40e-2)^3) / (water_props.kinematic_viscosity * water_props.thermal_diffusivity);
    
    Nu_smooth = Q * H  / (2 * water_props.thermal_conductivity * (temperature_data.T_internal - temperature_data.T_up));
    Nu_rough = Q * H  / (2 * water_props.thermal_conductivity * (temperature_data.T_down - temperature_data.T_internal));
    
    ind_acq = temperature_data.ind_acq;
    
    processed_data.t_probes = temperature_data.t(ind_acq) - temperature_data.t(ind_acq(1));
    
    processed_data.T_bas1 = temperature_data.T_bas1(ind_acq);
    processed_data.T_bas2 = temperature_data.T_bas2(ind_acq);
    processed_data.T_bas3 = temperature_data.T_bas3(ind_acq);
    
    processed_data.T_haut1 = temperature_data.T_haut1(ind_acq);
    processed_data.T_haut2 = temperature_data.T_haut2(ind_acq);
    
    
    processed_data.T_int = T_int;
    processed_data.T_room = T_room;
    processed_data.T_hot = temperature_data.T_down;
    processed_data.T_cold = temperature_data.T_up;
    
    processed_data.Ra = Ra;
    processed_data.Nu = Nu;
    
    processed_data.Ra_smooth = Ra_smooth;
    processed_data.Ra_rough = Ra_rough;
    
    processed_data.Nu_smooth = Nu_smooth;
    processed_data.Nu_rough = Nu_rough;
    processed_data.Ra = Ra;
    
    
    %% Process data
    
    if ~exist(paths.output_file,'file')
        img_name = [imgs_processed(1).folder '/' imgs_processed(1).name];
        load(img_name)

        im_T = im_T(sub_roi(1):sub_roi(2),sub_roi(3):sub_roi(4));
        size_px = size(im_T);
        size_cm = 100*(size_px ./ gamma_px);
        
        
        processed_data.position_measures_x_px = position_single_point_x_px;
        processed_data.position_measures_z_px = position_single_point_z_px;
        
        processed_data.position_x_cm = 100*(processed_data.position_measures_x_px ./ gamma_px);
        processed_data.position_z_cm = 100*(processed_data.position_measures_z_px ./ gamma_px);
        
        
        processed_data.middle_position_z_px_notch = [round((position_BL_z_cm ./ 100) .* gamma_px)...
            round(((size_cm(1) - position_BL_z_cm) ./ 100) .* gamma_px)];
        processed_data.middle_position_z_px_plot = [processed_data.middle_position_z_px_notch(1:length(processed_data.middle_position_z_px_notch)/2)+...
            height_rugosity_px processed_data.middle_position_z_px_notch(length(processed_data.middle_position_z_px_notch)/2:end)];
        
        processed_data.middle_position_x_px_notch = position_notch_x_px;
        processed_data.middle_position_x_px_plot = position_plot_x_px;
        
        processed_data.middle_position_x_cm_notch = 100*(processed_data.middle_position_x_px_notch ./ gamma_px);
        processed_data.middle_position_x_cm_plot = 100*(processed_data.middle_position_x_px_plot ./ gamma_px);
        
        
        % Read temperature field, build the mean image
        count_loop = 0;
        T_mean = zeros(1,imgs_to_process);
        for jj = img_start(ii):img_end(ii)
            count_loop = count_loop + 1;
            
            img_name = [imgs_processed(jj).folder '/' imgs_processed(jj).name];
            load(img_name)
            
            im_T = im_T(sub_roi(1):sub_roi(2),sub_roi(3):sub_roi(4));
            
            T_mean(count_loop) = mean(mean(im_T(round(0.5*(size_px(1)-T_roi_size):0.5*(size_px(1)+T_roi_size)),...
                round(0.5*(size_px(2)-T_roi_size):0.5*(size_px(2)+T_roi_size)))));
%             T_mean(count_loop) = mean(mean(im_T));
            
            for kk = 1:length(position_BL_z_cm)
                position_notch_z_px = round(position_BL_z_cm(kk)*10^(-2) * gamma_px) + 3;
                position_plot_z_px = position_notch_z_px + height_rugosity_px + 3;

                processed_data.(['middle_position_notch_T_' num2str(kk)])(count_loop) = im_T(position_notch_z_px,position_notch_x_px);
                processed_data.(['middle_position_plot_T_' num2str(kk)])(count_loop) = im_T(position_plot_z_px,position_plot_x_px);
            end
            
            
            if count_loop == 1
                im_T_mean = im_T;
            else
                im_T_mean = im_T_mean + im_T;
            end
        end
        im_T_mean = im_T_mean ./ imgs_to_process;
        
        
        % Make most measures with subtracted mean field
        count_loop = 0;
        for jj = img_start(ii):img_end(ii)
            count_loop = count_loop + 1;
            
            img_name = [imgs_processed(jj).folder '/' imgs_processed(jj).name];
            load(img_name)
            
            im_T = im_T(sub_roi(1):sub_roi(2),sub_roi(3):sub_roi(4));
            im_T = im_T - im_T_mean;
            
            if count_loop == 1
                positions_over_rugosity = 7:2*height_rugosity_px:size_px(1);
                positions_over_trench = 7+height_rugosity_px:2*height_rugosity_px:size_px(1);
            end
            
            
            % Horizontal and vertical temperature profiles
            for kk = 1:length(position_horizontal_profiles_cm)
                y_px = round(position_horizontal_profiles_cm(kk)*10^(-2) * gamma_px);
                
                field_name = ['y_' num2str(kk)];
                processed_data.(field_name)(count_loop,:) = im_T(y_px,:);
            end
            for kk = 1:length(position_vertical_profiles_cm)
                x_px = round(position_vertical_profiles_cm(kk)*10^(-2) * gamma_px);
                
                field_name = ['x_' num2str(kk)];
                processed_data.(field_name)(:,count_loop) = im_T(:,x_px);
            end
            
            for kk = 1:length(position_single_point_x_px)
                for ll = 1:length(position_single_point_z_px)
                    data_T = im_T(position_single_point_x_px(kk)-half_rugosity:position_single_point_x_px(kk)+half_rugosity,...
                        position_single_point_z_px(ll)-half_rugosity:position_single_point_z_px(ll)+half_rugosity);
                    
%                     processed_data.(['single_point_T_' num2str(kk) '_' num2str(ll)])(count_loop) = im_T(position_single_point_x_px(kk),...
%                         position_single_point_z_px(ll));
                    processed_data.(['single_point_T_' num2str(kk) '_' num2str(ll)])(count_loop) = im_T(position_single_point_x_px(kk),...
                        position_single_point_z_px(ll));
                    processed_data.(['single_point_hist_' num2str(kk) '_' num2str(ll)])(count_loop,:) = reshape(data_T,1,81);
                end
            end
            
            
            if count_loop == 1
                im_T_rms = im_T.^2;
            else
                im_T_rms = im_T_rms + im_T.^2;
            end
        end
        im_T_rms = sqrt( im_T_rms ./ imgs_to_process);
        
        processed_data.horizontal_profiles_cm = position_horizontal_profiles_cm;
        processed_data.vertical_profiles_cm = position_vertical_profiles_cm;
        
        processed_data.T_mean = T_mean;
        
        processed_data.im_T_mean = im_T_mean;
        processed_data.im_T_rms = im_T_rms;
        
        save(paths.output_file,'processed_data','-v7.3')
    else
        load(paths.output_file)
        
        size_px = size(processed_data.im_T_mean);
        size_cm = 100*(size_px ./ gamma_px);
    end
    
    all_processed_data(count_main_loop) = processed_data;
end

%%
for ii = 1:length(all_processed_data)
    field_names = fieldnames(all_processed_data(ii));
    
   for jj = 1:length(field_names)
       if length(all_processed_data(ii).(field_names{jj})) == 2595
           ind_zero = all_processed_data(ii).(field_names{jj}) == 0;
           
           all_processed_data(ii).(field_names{jj})(ind_zero) = [];
       end
   end
end

%%
for ii = 1:length(all_processed_data)
    im_rms = all_processed_data(ii).im_T_rms;
    
    im_rms = im_rms.^2;
    im_rms = (im_rms .* length(imgs_processed)) ./ (length(imgs_processed) - 800);
    im_rms = sqrt(im_rms);
    
    all_processed_data(ii).im_T_rms = im_rms;
%     all_processed_data(ii).im_T_rms = sqrt((all_processed_data(ii).im_T_rms.^2 .* length(imgs_processed)) ./ (length(imgs_processed) - 800));

end


%% Temperature calibration

figure('Position',[100 100 1000 1000]),
hold on,
plot([temp_calib.T_int],[temp_calib.I],'Color',col_order(1,:),'Marker',marker_order{1},'LineWidth',2,'MarkerSize',15,'LineStyle','none')
y_fit = linspace(0.8,1.5,100);
x_fit = polyval(temp_calib.fit,y_fit);
plot(x_fit,y_fit,'Color',col_order(1,:),'LineWidth',2)
xlim([min([temp_calib.T_int])-3 max([temp_calib.T_int])+3])
xlabel('$T_{int}$ ($^{\circ}$C)')
ylabel('$I$')

if strcmp(save_imgs,'on')
    export_fig([path_output 'temperature_calibration.pdf'],'-transparent','-pdf')
end


%% Mean temperature evolution

figure('Position',[100 100 1000 1000]),
hold on,
for ii = 1:length(all_processed_data)
    t = 0:dt:(length(all_processed_data(ii).T_mean)-1)*dt;
    T_probe = all_processed_data(ii).T_int;
    
    plot(t,all_processed_data(ii).T_mean,'LineWidth',2,'Color',col_order(ii,:))
%     plot([min(t) max(t)],[T_probe T_probe],'k','LineWidth',2)
end
% axis([0 260 32 34])
legend('run01','run02','run03','run04','location','SouthEast')
xlabel('$t$ (s)')
ylabel('$<T>$ ($^{\circ}$C)')

if strcmp(save_imgs,'on')
    export_fig([path_output 'mean_T_evolution.pdf'],'-transparent','-pdf')
end


%% Temperature difference

figure('Position',[100 100 1000 1000]),
hold on,
t = 0:dt:(length(all_processed_data(1).T_mean)-1)*dt;
plot(t,all_processed_data(2).T_mean-all_processed_data(1).T_mean,'LineWidth',2,'Color',col_order(1,:))

text(100,0.3,['$\sigma_{T}=' num2str(std(all_processed_data(2).T_mean-all_processed_data(1).T_mean),'%.02f')...
        '^{\circ}$C'],'FontSize',20)
    
axis([0 260 0 1])
xlabel('$t$ (s)')
ylabel('$<T_2> - <T_1>$ ($^{\circ}$C)')

if strcmp(save_imgs,'on')
    export_fig([path_output 'T_difference.pdf'],'-transparent','-pdf')
end


%% Temperature probe data

figure('Position',[100 100 1500 1000]),

[ha,~] = tight_subplot(1,length(all_processed_data),[.05 0.01],[.1 .05],[.07 .01]);

for ii = 1:length(all_processed_data)
    axes(ha(ii))
    t_max = max(all_processed_data(ii).t_probes);
    
    hold on,
    plot(all_processed_data(ii).t_probes,all_processed_data(ii).T_bas1,'LineWidth',2,'Color',col_order(1,:))
    plot(all_processed_data(ii).t_probes,all_processed_data(ii).T_bas2,'LineWidth',2,'Color',col_order(2,:))
    plot(all_processed_data(ii).t_probes,all_processed_data(ii).T_bas3,'LineWidth',2,'Color',col_order(3,:))
    plot(all_processed_data(ii).t_probes,all_processed_data(ii).T_haut1,'LineWidth',2,'Color',col_order(4,:))
    plot(all_processed_data(ii).t_probes,all_processed_data(ii).T_haut2,'LineWidth',2,'Color',col_order(5,:))
    
    plot([0 t_max],[all_processed_data(ii).T_int all_processed_data(ii).T_int],...
        'k','LineWidth',2)
    
    text(0.3*t_max,all_processed_data(ii).T_int-5,['$T_{int}=' num2str(all_processed_data(ii).T_int,'%.02f')...
        '^{\circ}$C'],'FontSize',20)
    text(0.3*t_max,all_processed_data(ii).T_int-7,['$T_{room}=' num2str(all_processed_data(ii).T_room,'%.02f')...
        '^{\circ}$C'],...
        'FontSize',20)
      
    
    title(['run' num2str(ii,'%02.f')])
    
    ax=gca;
    x_ticks = 0:60:t_max;
    y_ticks = 15:5:45;
    
    ax.XTick = x_ticks;
    ax.XTickLabel = num2str(x_ticks');
    
    axis([0 t_max 15 45])
    
    xlabel('$t$ (s)')
    if ii == 1
        ax.YTick = y_ticks;
        ax.YTickLabel = num2str(y_ticks');
        ylabel('$T$ ($^{\circ}$C)')
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'temperature_probes_signal.pdf'],'-transparent','-pdf')
end


%% Dye spectra

rhB_em_file = [paths.dyes 'rhB_em.csv'];
rhB_abs_file = [paths.dyes 'rhB_abs.csv'];
rh110_em_file = [paths.dyes 'rh110_em.csv'];
rh110_abs_file = [paths.dyes 'rh110_abs.csv'];

fid = fopen(rhB_abs_file);
rhB_abs = textscan(fid,repmat('%f',1,2),'delimiter',';');
fclose(fid);

fid = fopen(rh110_abs_file);
rh110_abs = textscan(fid,repmat('%f',1,2),'delimiter',';');
fclose(fid);

fid = fopen(rhB_em_file);
rhB_em = textscan(fid,repmat('%f',1,2),'delimiter',';');
fclose(fid);

fid = fopen(rh110_em_file);
rh110_em = textscan(fid,repmat('%f',1,2),'delimiter',';');
fclose(fid);


figure('Position',[100 100 1000 1000]),
hold on,
plot(rhB_abs{1},rhB_abs{2},'LineWidth',2,'LineStyle','--','Color',col_order(4,:))
plot(rh110_em{1},rh110_em{2},'LineWidth',2,'Color',col_order(5,:))
axis([400 650 0 1])
% yl = ylim;
% plot([488 488],[min(yl) max(yl)],'b','LineWidth',2)
legend('Rh B, abs','Rh 110, em')
ylabel('Normalized intensity')
xlabel('$\lambda$ (nm)')


if strcmp(save_imgs,'on')
    export_fig([path_output 'spectral_em110_absB.pdf'],'-transparent','-pdf')
end


%% Profiles position

if strcmp(disp_imgs,'on')
    im_disp = scaled2rgb(processed_data.im_T_rms,cmap,[0 T_rms_display]);
    
    figure,
    imshow(im_disp,[]),
    hold on,
    ax=gca;
    ax.YDir = 'normal';
    
    size_img_px = size(im_disp);
    size_img_cm = 100 .* (size_img_px ./ gamma_px);
    
    xticks_desired = 0:10:size_img_cm(2);
    yticks_desired = 0:10:size_img_cm(1);
    
    xticks = (xticks_desired./100) .* gamma_px;
    yticks = (yticks_desired./100) .* gamma_px;
    
    xticks_label = num2str(xticks_desired');
    yticks_label = num2str(yticks_desired');
    
    ax = gca;
    ax.YDir = 'normal';
    ax.Visible = 'On';
    ax.YTickLabel = yticks_label;
    ax.YTick = yticks;
    ax.XTickLabel = xticks_label;
    ax.XTick = xticks;
    
    xlabel('$x \, $(cm)')
    ylabel('$z \, $(cm)')
    
    colormap(cmap),
    h2 = colorbar;
    ylabel(h2, 'RMS($T-\overline{T}$) (K)','interpreter','latex')
    caxis([0 T_rms_display])
    
    
    for jj = 1:length(temperature_profiles.horizontal_px)
        y_px = temperature_profiles.horizontal_px(jj);
        plot([1 size_img_px(2)],[y_px y_px],'g','LineWidth',2),
    end
    for jj = 1:length(temperature_profiles.vertical_px)
        x_px = temperature_profiles.vertical_px(jj);
        plot([x_px x_px],[1 size_img_px(1)],'y','LineWidth',2),
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_position.pdf'],'-transparent','-pdf')
end


%% Mean image

ind = 4;

im_disp = scaled2rgb(all_processed_data(ind).im_T_mean - all_processed_data(ind).T_int,cmap,[-4 4]);

figure,
imshow(im_disp,[]),
hold on,
ax=gca;
ax.YDir = 'normal';

size_img_px = size(im_disp);
size_img_cm = 100 .* (size_img_px ./ gamma_px);

xticks_desired = 0:10:size_img_cm(2);
yticks_desired = 0:10:size_img_cm(1);

xticks = (xticks_desired./100) .* gamma_px;
yticks = (yticks_desired./100) .* gamma_px;

xticks_label = num2str(xticks_desired');
yticks_label = num2str(yticks_desired');

ax = gca;
ax.YDir = 'normal';
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;

xlabel('$x \, $(cm)')
ylabel('$z \, $(cm)')

colormap(cmap),
h2 = colorbar;
ylabel(h2,'$\overline{T}-T_i$ (K)','interpreter','latex')

caxis([-T_scale_span T_scale_span])


if strcmp(save_imgs,'on')
    export_fig([path_output 'mean_image.pdf'],'-pdf')
end



%% RMS image

im_disp = scaled2rgb(all_processed_data(3).im_T_rms,cmap,[0 T_rms_display]);

figure,
imshow(im_disp,[]),
hold on,
ax=gca;
ax.YDir = 'normal';

size_img_px = size(im_disp);
size_img_cm = 100 .* (size_img_px ./ gamma_px);

xticks_desired = 0:10:size_img_cm(2);
yticks_desired = 0:10:size_img_cm(1);

xticks = (xticks_desired./100) .* gamma_px;
yticks = (yticks_desired./100) .* gamma_px;

xticks_label = num2str(xticks_desired');
yticks_label = num2str(yticks_desired');

ax = gca;
ax.YDir = 'normal';
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;

xlabel('$x \, $(cm)')
ylabel('$z \, $(cm)')

colormap(cmap),
h2 = colorbar;
ylabel(h2, 'RMS($T-\overline{T}$) (K)','interpreter','latex')

caxis([0 T_rms_display])


if strcmp(save_imgs,'on')
    export_fig([path_output 'rms_image_normal.pdf'],'-transparent','-pdf')
end


%% RMS data

position_rms_profiles_cm = [2 10 20 30 40];
position_rms_profiles_px = round((position_rms_profiles_cm*10^(-2)) .* gamma_px);

for ii = 1:length(position_rms_profiles_px)
    
    field_name = ['profile_' num2str(ii)];
    
    for jj = 1:length(all_processed_data)
        im_RMS = all_processed_data(jj).im_T_rms;
        profiles_rms.(field_name)(:,jj) = im_RMS(:,position_rms_profiles_px(ii));
    end
end


figure('Position',[100 100 2000 1000]),

[ha,~] = tight_subplot(1,length(position_rms_profiles_px),[.05 0.01],[.1 .05],[.07 .01]);

for ii = 1:length(position_rms_profiles_px)
    
    axes(ha(ii))
    hold on,

    title(['x=' num2str(position_rms_profiles_cm(ii),'%02.f') 'cm'])
    field_name = ['profile_' num2str(ii)];
    
    for jj = 1:length(all_processed_data)
        z = ((1:length(profiles_rms.(field_name)(:,jj))) ./ gamma_px) * 100;
        plot(z,profiles_rms.(field_name)(:,jj),'LineWidth',2,'Color',col_order(jj,:))
    end
    
    x_ticks = 0:10:max(z);
    y_ticks = 0:0.1:2;
    ax = gca;
    ax.XTick = x_ticks;
    ax.XTickLabel = num2str(x_ticks');

    xlabel('$z$ (cm)')
    if ii == 1
        ax.YTick = y_ticks;
        ax.YTickLabel = num2str(y_ticks');
        ylabel('$RMS(T)$ ($^{\circ}$C)')
        
        legend('run01','run02','run03','run04')
    end
    
    axis([0 max(z) 0 0.8])
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'rms_profiles_zoomed.pdf'],'-transparent','-pdf')
end


%% RMS data

ind_run = 4;

size_area_cm = 10;
size_area_px =  round(size_area_cm ./ 100 .* gamma_px);

im_RMS = all_processed_data(ind_run).im_T_rms;

[X,Y] = meshgrid(1:size_area_px,1:size_area_px);
Z = im_RMS(1:size_area_px,1:size_area_px);
% Z = im_RMS(size_px(1)-size_area_px+1:size_px(1),size_px(2)-size_area_px+1:size_px(2));
% Z = im_RMS(size_px(1)-size_area_px+1:size_px(1),1:size_area_px);
Z = im_RMS(1:size_area_px,size_px(2)-size_area_px+1:size_px(2));

figure('Position',[100 100 2000 1000]),
surf(X,Y,Z)

x_ticks = (1:size_area_cm) ./ 100 .* gamma_px;
x_ticks_label = num2str((1:size_area_cm)');

y_ticks = x_ticks;
y_ticks_label = x_ticks_label;

ax=gca;
ax.XTick = x_ticks;
ax.XTickLabel = x_ticks_label;
ax.YTick = y_ticks;
ax.YTickLabel = y_ticks_label;

xlabel('x (cm)')
ylabel('z (cm)')
zlabel('RMS(T) ($^{\circ}$C)')

view(130,45)

if strcmp(save_imgs,'on')
    export_fig([path_output 'rms_3dview_' num2str(ind_run) '.png'],'-transparent','-png')
end


%% Single points measurement positions

T_rms_display = 0.5;

if strcmp(disp_imgs,'on')
    im_disp = scaled2rgb(processed_data.im_T_rms,cmap,[0 T_rms_display]);
    
    for jj = 1:length(processed_data.position_measures_x_px)
        for kk = 1:length(processed_data.position_measures_z_px)
            position_measures_x = processed_data.position_measures_x_px(jj);
            position_measures_z = processed_data.position_measures_z_px(kk);
            
            im_disp(position_measures_x-half_rugosity:position_measures_x+half_rugosity,...
                position_measures_z-half_rugosity:position_measures_z+half_rugosity,:) = 0;
        end
    end
    
    figure,
    imshow(im_disp,[]),
    hold on,
    ax=gca;
    ax.YDir = 'normal';
    
    size_img_px = size(im_disp);
    size_img_cm = 100 .* (size_img_px ./ gamma_px);
    
    xticks_desired = 0:10:size_img_cm(2);
    yticks_desired = 0:10:size_img_cm(1);
    
    xticks = (xticks_desired./100) .* gamma_px;
    yticks = (yticks_desired./100) .* gamma_px;
    
    xticks_label = num2str(xticks_desired');
    yticks_label = num2str(yticks_desired');
    
    ax = gca;
    ax.YDir = 'normal';
    ax.Visible = 'On';
    ax.YTickLabel = yticks_label;
    ax.YTick = yticks;
    ax.XTickLabel = xticks_label;
    ax.XTick = xticks;
    
    xlabel('$x \, $(cm)')
    ylabel('$z \, $(cm)')
    
    colormap(cmap),
    h2 = colorbar;
    ylabel(h2, 'RMS($T-\overline{T}$) (K)','interpreter','latex')
    
    caxis([0 T_rms_display])
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'single_point_positions.jpg'],'-jpg')
end


%% Plots - horizontal bot
clear profiles

ind_run = 4;

ind = 1:3;

for ii = ind
    profile = all_processed_data(ind_run).(['y_' num2str(ii)]);
    profile = scaled2rgb(profile,cmap,[-T_scale_span T_scale_span]);
    
    profiles.(['hor_' num2str(ii)]) = profile;
end


figure('position',[100 100 2000 1000]);

[ha,~] = tight_subplot(1,3,[.05 -0.1],[.1 .05],[.01 .01]);
xticks_desired = 0:10:size_cm(2);
yticks_desired = 0:10:round((imgs_to_process-1)/10);

xticks = (xticks_desired./100) .* gamma_px;
yticks = yticks_desired.*10 + 1;

xticks_label = num2str(xticks_desired');
yticks_label = num2str(yticks_desired');

axes(ha(1))
imshow(profiles.hor_1,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),

xlabel('$x \, $(cm)')
ylabel('$t \, $(s)')
title(['z=' num2str(processed_data.horizontal_profiles_cm(1)) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])


axes(ha(2))
imshow(profiles.hor_2,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),

xlabel('$x \, $(cm)')
title(['z=' num2str(processed_data.horizontal_profiles_cm(2)) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = [];
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])

axes(ha(3))
imshow(profiles.hor_3,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),
h2 = colorbar;
ylabel(h2, '$T-\overline{T}$ (K)','interpreter','latex')

xlabel('$x \, $(cm)')
title(['z=' num2str(num2str(processed_data.horizontal_profiles_cm(3))) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = [];
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_horizontal_bot_4.jpg'],'-jpg')
end


%% Plots - horizontal top
clear profiles

ind_run = 4;

ind = length(all_processed_data(ind_run).horizontal_profiles_cm)-2:length(all_processed_data(ind_run).horizontal_profiles_cm);

count_loop = 0;
for ii = ind
    count_loop = count_loop + 1;
    
    profile = all_processed_data(ind_run).(['y_' num2str(ii)]);
    profile = scaled2rgb(profile,cmap,[-T_scale_span T_scale_span]);
    
    profiles.(['hor_' num2str(count_loop)]) = profile;
end


figure('position',[100 100 2000 1000]);

[ha,~] = tight_subplot(1,3,[.05 -0.1],[.1 .05],[.01 .01]);
xticks_desired = 0:10:size_cm(2);
yticks_desired = 0:10:round((imgs_to_process-1)/10);

xticks = (xticks_desired./100) .* gamma_px;
yticks = yticks_desired.*10 + 1;

xticks_label = num2str(xticks_desired');
yticks_label = num2str(yticks_desired');

axes(ha(1))
imshow(profiles.hor_1,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),

xlabel('$x \, $(cm)')
ylabel('$t \, $(s)')
title(['z=' num2str(processed_data.horizontal_profiles_cm(ind(1))) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])


axes(ha(2))
imshow(profiles.hor_2,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),

xlabel('$x \, $(cm)')
title(['z=' num2str(processed_data.horizontal_profiles_cm(ind(2))) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = [];
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])

axes(ha(3))
imshow(profiles.hor_3,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),
h2 = colorbar;
ylabel(h2, '$T-\overline{T}$ (K)','interpreter','latex')

xlabel('$x \, $(cm)')
title(['z=' num2str(num2str(processed_data.horizontal_profiles_cm(ind(3)))) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = [];
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_horizontal_top_4.jpg'],'-jpg')
end


%% Plots - vertical left
clear profiles

ind = 1:3;

count_loop = 0;
for ii = ind
    count_loop = count_loop + 1;
    
    profile = processed_data.(['x_' num2str(ii)]);
    profile = scaled2rgb(profile,cmap,[-T_scale_span T_scale_span]);
    
    profiles.(['ver_' num2str(count_loop )]) = profile;
end


figure('position',[100 100 1000 1200]);

[ha,~] = tight_subplot(3,1,[.05 -0.1],[.1 .05],[.01 .01]);
yticks_desired = 0:10:size_cm(2);
xticks_desired = 0:30:round((imgs_to_process-1)/10);

yticks = (yticks_desired./100) .* gamma_px;
xticks = xticks_desired.*10 + 1;

yticks_label = num2str(yticks_desired');
xticks_label = num2str(xticks_desired');

axes(ha(1))
imshow(profiles.ver_1,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),
h2 = colorbar;
ylabel(h2, '$T-\overline{T}$ (K)','interpreter','latex')

ylabel('$z \, $(cm)')
title(['x=' num2str(processed_data.vertical_profiles_cm(ind(1))) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = [];
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])


axes(ha(2))
imshow(profiles.ver_2,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),
h2 = colorbar;
ylabel(h2, '$T-\overline{T}$ (K)','interpreter','latex')

ylabel('$z \, $(cm)')
title(['x=' num2str(processed_data.vertical_profiles_cm(ind(2))) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = [];
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])

axes(ha(3))
imshow(profiles.ver_3,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),
h2 = colorbar;
ylabel(h2, '$T-\overline{T}$ (K)','interpreter','latex')

xlabel('$t \, $(s)')
ylabel('$z \, $(cm)')
title(['x=' num2str(num2str(processed_data.vertical_profiles_cm(ind(3)))) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_vertical_left.jpg'],'-jpg')
end


%% Plots - vertical right

ind = length(processed_data.horizontal_profiles_cm)-2:length(processed_data.horizontal_profiles_cm);
for ii = ind
    profile = processed_data.(['x_' num2str(ii)]);
    profile = scaled2rgb(profile,cmap,[-T_scale_span T_scale_span]);
    
    profiles.(['ver_' num2str(ii)]) = profile;
end


figure('position',[100 100 1000 1200]);

[ha,~] = tight_subplot(3,1,[.05 -0.1],[.1 .05],[.01 .01]);
yticks_desired = 0:10:size_cm(2);
xticks_desired = 0:30:round((imgs_to_process-1)/10);

yticks = (yticks_desired./100) .* gamma_px;
xticks = xticks_desired.*10 + 1;

yticks_label = num2str(yticks_desired');
xticks_label = num2str(xticks_desired');

axes(ha(1))
imshow(profiles.ver_41,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),
h2 = colorbar;
ylabel(h2, '$T-\overline{T}$ (K)','interpreter','latex')

ylabel('$z \, $(cm)')
title(['x=' num2str(processed_data.vertical_profiles_cm(41)) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = [];
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])


axes(ha(2))
imshow(profiles.ver_40,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),
h2 = colorbar;
ylabel(h2, '$T-\overline{T}$ (K)','interpreter','latex')

ylabel('$z \, $(cm)')
title(['x=' num2str(processed_data.vertical_profiles_cm(40)) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = [];
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])

axes(ha(3))
imshow(profiles.ver_39,[])
ax = gca;
ax.YDir = 'normal';
colormap(cmap),
h2 = colorbar;
ylabel(h2, '$T-\overline{T}$ (K)','interpreter','latex')

xlabel('$t \, $(s)')
ylabel('$z \, $(cm)')
title(['x=' num2str(num2str(processed_data.vertical_profiles_cm(39))) 'cm'])

ax = gca;
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;
caxis([-T_scale_span T_scale_span])


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_vertical_right.jpg'],'-jpg')
end


%% Single point histograms

ind_run = 4;

figure('position',[100 100 1200 1200]);
[ha,~] = tight_subplot(3,3,[.03 .03],[.08 .01],[.08 .01]);
count_loop = 0;
for ii = [3 2 1]
    for jj = [1 2 3]
        count_loop = count_loop + 1;
        
        axes(ha(count_loop))
        data = all_processed_data(ind_run).(['single_point_hist_' num2str(ii) '_' num2str(jj)]);
        histogram(data,'BinWidth',0.1)
        axis([-5 5 1 1e5])
        ax=gca;
        ax.YScale = 'log';
        
        if count_loop > 6
            xlabel('$T-\overline{T}$ ($^{\circ}$ C)')
            ax.XTick = -5:1:5;
        else
            ax.XTick = [];
        end
        
        if count_loop  == 1 || count_loop  == 4 || count_loop  == 7
            ax.YTick = [1 10 100 1e3 1e4 1e5];
            ylabel('count')
        else
            ax.YTick = [];
        end
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'single_point_hist_data' num2str(ind_run) '.pdf'],'-pdf')
end


%% Single point histograms - comparison

b_width = 0.1;
x_values = 0:b_width:5;

figure('position',[100 100 1500 1000]);
[ha,~] = tight_subplot(1,3,[.03 .03],[.08 .01],[.08 .01]);

axes(ha(1));
hold on,

data_1 = abs(processed_data.('single_point_hist_1_1'));
data_2 = abs(processed_data.('single_point_hist_3_3'));

size_data = size(data);

data_1 = reshape(data_1,size_data(1)*size_data(2),1);
data_2 = reshape(data_2,size_data(1)*size_data(2),1);

pd_1 = fitdist(data_1,'kernel','Kernel','normal','width',b_width);
pd_2 = fitdist(data_2,'kernel','Kernel','normal','width',b_width);

y_1 = pdf(pd_1,x_values);
y_2 = pdf(pd_2,x_values);

hold on
plot(x_values,y_1,'LineWidth',2,'Color',col_order(1,:))
plot(x_values,y_2,'LineWidth',2,'Color',col_order(2,:))
ax=gca;
ax.YScale='log';
axis([0 5 1e-5 10])
xlabel('$|T-\overline{T}|$ (K)')
ylabel('p.d.f.')
legend('Rough','Smooth')


axes(ha(2));
hold on,

data_1 = abs(processed_data.('single_point_hist_1_2'));
data_2 = abs(processed_data.('single_point_hist_3_2'));

size_data = size(data);

data_1 = reshape(data_1,size_data(1)*size_data(2),1);
data_2 = reshape(data_2,size_data(1)*size_data(2),1);

pd_1 = fitdist(data_1,'kernel','Kernel','normal','width',b_width);
pd_2 = fitdist(data_2,'kernel','Kernel','normal','width',b_width);

y_1 = pdf(pd_1,x_values);
y_2 = pdf(pd_2,x_values);

hold on
plot(x_values,y_1,'LineWidth',2,'Color',col_order(1,:))
plot(x_values,y_2,'LineWidth',2,'Color',col_order(2,:))
ax=gca;
ax.YScale='log';
ax.YTickLabel = [];
axis([0 5 1e-5 10])
xlabel('$|T-\overline{T}|$ (K)')
legend('Rough','Smooth')


axes(ha(3));
hold on,

data_1 = abs(processed_data.('single_point_hist_1_3'));
data_2 = abs(processed_data.('single_point_hist_3_1'));

size_data = size(data);

data_1 = reshape(data_1,size_data(1)*size_data(2),1);
data_2 = reshape(data_2,size_data(1)*size_data(2),1);

pd_1 = fitdist(data_1,'kernel','Kernel','normal','width',b_width);
pd_2 = fitdist(data_2,'kernel','Kernel','normal','width',b_width);

y_1 = pdf(pd_1,x_values);
y_2 = pdf(pd_2,x_values);

hold on
plot(x_values,y_1,'LineWidth',2,'Color',col_order(1,:))
plot(x_values,y_2,'LineWidth',2,'Color',col_order(2,:))
ax=gca;
ax.YScale='log';
ax.YTickLabel = [];
axis([0 5 1e-5 10])
xlabel('$|T-\overline{T}|$ (K)')
legend('Rough','Smooth')


if strcmp(save_imgs,'on')
    export_fig([path_output 'single_point_hist_comparison.pdf'],'-pdf')
end


%% Single point temporal T

ind_run = 1;

figure('position',[100 100 1200 1200]);
[ha,~] = tight_subplot(3,3,[.03 .03],[.08 .01],[.08 .01]);
count_loop = 0;
for ii = [3 2 1]
    for jj = [1 2 3]
        count_loop = count_loop + 1;
        
        axes(ha(count_loop))
        data = all_processed_data(ind_run).(['single_point_T_' num2str(ii) '_' num2str(jj)]);
        t = 0.1 .* (1:length(data));
        plot(t,data)
        axis([0 180 -4 4])
        ax=gca;
        
        if count_loop > 6
            xlabel('t (s)')
        else
            ax.XTick = [];
        end
        if count_loop  == 1 || count_loop  == 4 || count_loop  == 7
            ylabel('$T-\overline{T}$ ($^{\circ}$C)')
        else
            ax.YTick = [];
        end
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'single_point_temp_data_' num2str(ind_run) '.pdf'],'-pdf')
end


%% Plume detection

ii = 1;
jj = 3;
position_measures_x = processed_data.position_measures_x_px(ii);
position_measures_z = processed_data.position_measures_z_px(jj);
data = processed_data.(['single_point_T_' num2str(ii) '_' num2str(jj)]);
data_grad = gradient(data);
t = (1:length(data));

T_rms = processed_data.im_T_rms(position_measures_x,position_measures_z);

ind_t0 = diff(data) > 0.35 .* T_rms;
ind_t0 = [false ind_t0];

ind_grad = abs(data_grad) > 0.1;

ind_t0 = ind_t0 & ind_grad;


figure('position',[100 100 1000 1000]);
hold on,
plot(t,data)
ylabel('$T-\overline{T}$ ($^{\circ}$C)')
xlabel('t (s)')
% plot(t(ind_t0),data_grad(ind_t0),'^r')
axis([0 2600 -4 4])

% export_fig([path_output 'single_point_temperature.pdf'],'-pdf')
%%

ind_t0_f = find(ind_t0);
ind_grad_0 = abs(data_grad) < 3e-3;

ind_t1 = false(1,length(data));
ind_t2 = false(1,length(data));
for ii = 1:length(ind_t0_f)
    ind_min = max(ind_t0_f(ii)-20,0):ind_t0_f(ii);
    local_min = find(islocalmin(data_grad(ind_min)),1,'last');
    
    ind_max = ind_t0_f(ii):min(length(ind_t0_f),ind_t0_f(ii)+20);
    local_max = find(islocalmax(data_grad(ind_max)),1,'first');
    
    ind_t1(ind_t0_f(ii) - local_min - 1) = true;
    ind_t2(ind_t0_f(ii) + local_max - 1) = true;
end

%%

figure('position',[100 100 1000 1000]);
hold on,
plot(t,data)
plot(t(pks),data(pks),'bo')
% plot(t(ind_t2),data(ind_t2),'ro')
axis([0 2600 -5 5])



%%
mean_T = zeros(1,length(data)-1);
count_loop = 0;
for ii = 2:length(data)
    count_loop = count_loop + 1;
    
    mean_T(count_loop) =  0.5 .* (data(ii-1) + data(ii));
end

ind_T = mean_T > 0;
ind_T = [false ind_T];
%
ind_plumes = ind_T & ind_t0;

figure('position',[100 100 1000 1000]);
hold on,
plot(t,data)
plot(t(ind_t0),data(ind_t0),'^r')
axis([0 260 -4 4])


%% Profile histograms - horizontal bot

ind = 1:5;
num_profiles = length(ind);

figure('position',[100 100 500 1000]);
[ha,~] = tight_subplot(num_profiles,1,[.01 .01],[.1 .01],[.1 .01]);
count_loop = 0;
for ii = ind
    count_loop = count_loop + 1;
    
    axes(ha(count_loop))
    profile = processed_data.(['y_' num2str(ii)]);
    profile_position = processed_data.horizontal_profiles_cm(ii);
    histogram(profile,'BinWidth',0.1)
    axis([-5 5 1 1e6])
    ax=gca;
    ax.YScale = 'log';
    
    legend(['z=' num2str(profile_position) 'cm']);
    
    ylabel('count')
    if ii == ind(end)
        xlabel('$T - \overline{T}$ ($^{\circ}$C)')
    else
        ax.XTick = [];
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_histogram_hor_bot.png'],'-png')
end


%% Profile histograms - horizontal mid

ind = 18:22;
num_profiles = length(ind);

figure('position',[100 100 500 1000]);
[ha,~] = tight_subplot(num_profiles,1,[.01 .01],[.1 .01],[.1 .01]);
count_loop = 0;
for ii = ind
    count_loop = count_loop + 1;
    
    axes(ha(count_loop))
    profile = processed_data.(['y_' num2str(ii)]);
    profile_position = processed_data.horizontal_profiles_cm(ii);
    histogram(profile,'BinWidth',0.1)
    axis([-5 5 1 1e6])
    ax=gca;
    ax.YScale = 'log';
    
    legend(['z=' num2str(profile_position) 'cm']);
    
    ylabel('count')
    if ii == ind(end)
        xlabel('$T - \overline{T}$ ($^{\circ}$C)')
    else
        ax.XTick = [];
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_histogram_hor_mid.png'],'-png')
end


%% Profile histograms - horizontal top

ind = 37:41;
num_profiles = length(ind);

figure('position',[100 100 500 1000]);
[ha,~] = tight_subplot(num_profiles,1,[.01 .01],[.1 .01],[.1 .05]);
count_loop = 0;
for ii = ind
    count_loop = count_loop + 1;
    
    axes(ha(count_loop))
    profile = processed_data.(['y_' num2str(ii)]);
    profile_position = processed_data.horizontal_profiles_cm(ii);
    histogram(profile,'BinWidth',0.1)
    axis([-5 5 1 1e6])
    ax=gca;
    ax.YScale = 'log';
    
    legend(['z=' num2str(profile_position) 'cm']);
    
    ylabel('count')
    if ii == ind(end)
        xlabel('$T - \overline{T}$ ($^{\circ}$C)')
    else
        ax.XTick = [];
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_histogram_hor_top.png'],'-png')
end


%% Profile histograms - vertical left

ind = 1:5;
num_profiles = length(ind);

figure('position',[100 100 500 1000]);
[ha,~] = tight_subplot(num_profiles,1,[.01 .01],[.1 .01],[.1 .05]);
count_loop = 0;
for ii = ind
    count_loop = count_loop + 1;
    
    axes(ha(count_loop))
    profile = processed_data.(['x_' num2str(ii)]);
    profile_position = processed_data.vertical_profiles_cm(ii);
    histogram(profile,'BinWidth',0.1)
    axis([-7 7 1 1e6])
    ax=gca;
    ax.YScale = 'log';
    legend(['x=' num2str(profile_position) 'cm']);
    
    ylabel('count')
    if ii == ind(end)
        xlabel('$T - \overline{T}$ ($^{\circ}$C)')
    else
        ax.XTick = [];
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_histogram_ver_left.png'],'-png')
end


%% Profile histograms - vertical middle

ind = 18:22;
num_profiles = length(ind);

figure('position',[100 100 500 1000]);
[ha,~] = tight_subplot(num_profiles,1,[.01 .01],[.1 .01],[.1 .05]);
count_loop = 0;
for ii = ind
    count_loop = count_loop + 1;
    
    axes(ha(count_loop))
    profile = processed_data.(['x_' num2str(ii)]);
    profile_position = processed_data.vertical_profiles_cm(ii);
    histogram(profile,'BinWidth',0.1)
    axis([-7 7 1 1e6])
    ax=gca;
    ax.YScale = 'log';
    legend(['x=' num2str(profile_position) 'cm']);
    
    ylabel('count')
    if ii == ind(end)
        xlabel('$T - \overline{T}$ ($^{\circ}$C)')
    else
        ax.XTick = [];
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_histogram_ver_mid.png'],'-png')
end



%% Profile histograms - vertical right

ind = 37:41;
num_profiles = length(ind);

figure('position',[100 100 500 1000]);
[ha,~] = tight_subplot(num_profiles,1,[.01 .01],[.1 .01],[.1 .05]);
count_loop = 0;
for ii = ind
    count_loop = count_loop + 1;
    
    axes(ha(count_loop))
    profile = processed_data.(['x_' num2str(ii)]);
    profile_position = processed_data.vertical_profiles_cm(ii);
    histogram(profile,'BinWidth',0.1)
    axis([-7 7 1 1e6])
    ax=gca;
    ax.YScale = 'log';
    legend(['x=' num2str(profile_position) 'cm']);
    
    ylabel('count')
    if ii == ind(end)
        xlabel('$T - \overline{T}$ ($^{\circ}$C)')
    else
        ax.XTick = [];
    end
end


if strcmp(save_imgs,'on')
    export_fig([path_output 'profiles_histogram_ver_right.png'],'-png')
end


%% Profile histograms - middle position

b_width = 0.1;
x_values = 10:b_width:50;
x_ticks = min(x_values):5:max(x_values);
x_ticks_label = num2str(x_ticks');

num_pts = (length(processed_data.middle_position_z_px_plot)-1)/2;

figure('position',[100 100 1400 700]);
[ha,~] = tight_subplot(1,2,[.01 .02],[.13 .08],[.08 .05]);

axes(ha(1));
hold on,
legend_info = cell(1,num_pts);
for ii = 1:num_pts
    fieldname = ['middle_position_notch_T_' num2str(ii)];
    legend_info{ii} = ['$z=' num2str(10.*position_BL_z_cm(ii),'%.02f') 'mm$'];
    
    data = processed_data.(fieldname);
    
    pd = fitdist(data','kernel','Kernel','normal','width',b_width);
    
    y = pdf(pd,x_values);
    
    plot(x_values,y,'LineWidth',2)
end
title('Rough notch')
legend(legend_info,'location','NorthWest')
axis([25 40 1e-4 2])
ax=gca;
ax.YScale = 'log';
ax.YTickLabel = num2str([1e-4 1e-3 1e-2 1e-1 1 10]');
ax.XTick = x_ticks;
ax.XTickLabel = x_ticks_label;
xlabel('T ($^{\circ}$C)')
ylabel('p.d.f.')

axes(ha(2));
hold on,
for ii = 1:num_pts
    fieldname = ['middle_position_plot_T_' num2str(ii)];
    
    data = processed_data.(fieldname);
    
    pd = fitdist(data','kernel','Kernel','normal','width',b_width);
    
    y = pdf(pd,x_values);
    
    plot(x_values,y,'LineWidth',2)
end
title('Smooth')
axis([25 40 1e-4 2])
ax=gca;
ax.YScale = 'log';
ax.XTick = x_ticks;
ax.XTickLabel = x_ticks_label;
xlabel('T ($^{\circ}$C)')

if strcmp(save_imgs,'on')
    export_fig([path_output 'histograms_close_plates.png'],'-png')
end


%% Profile histograms - middle position

b_width = 0.1;
x_values = 10:b_width:50;
x_ticks = min(x_values):5:max(x_values);
x_ticks_label = num2str(x_ticks');

num_pts = length(processed_data.middle_position_z_px_notch)/2;

figure('position',[100 100 1400 700]);
[ha,~] = tight_subplot(1,2,[.01 .02],[.13 .08],[.08 .05]);

axes(ha(1));
hold on,
legend_info = cell(1,num_pts);
for ii = 1:num_pts
    fieldname = ['middle_position_notch_T_' num2str(ii)];
    legend_info{ii} = ['$z=' num2str(10.*position_BL_z_cm(ii),'%.02f') 'mm$'];
    
    data = processed_data.(fieldname);
    
    pd = fitdist(data','kernel','Kernel','normal','width',b_width);
    
    y = pdf(pd,x_values);
    
    plot(x_values,y,'LineWidth',2)
end
title('Rough notch')
legend(legend_info,'location','NorthWest')
axis([25 40 1e-4 2])
ax=gca;
ax.YScale = 'log';
ax.YTickLabel = num2str([1e-4 1e-3 1e-2 1e-1 1 10]');
ax.XTick = x_ticks;
ax.XTickLabel = x_ticks_label;
xlabel('T ($^{\circ}$C)')
ylabel('p.d.f.')


axes(ha(2));
hold on,
legend_info = cell(1,num_pts);
for ii = 1:num_pts
    fieldname = ['middle_position_plot_T_'  num2str(ii)];
    
    data = processed_data.(fieldname);
    
    pd = fitdist(data','kernel','Kernel','normal','width',b_width);
    
    y = pdf(pd,x_values);
    
    plot(x_values,y,'LineWidth',2)
end
title('Rough plot')
axis([25 40 1e-4 2])
ax=gca;
ax.YScale = 'log';
ax.XTick = x_ticks;
ax.XTickLabel = x_ticks_label;
xlabel('T ($^{\circ}$C)')

if strcmp(save_imgs,'on')
    export_fig([path_output 'histograms_close_plates_plot.png'],'-png')
end


%% Profile - middle position

im_disp = scaled2rgb(processed_data.im_T_mean - T_int,cmap,[-T_scale_span T_scale_span]);

figure,
imshow(im_disp,[]),
hold on,
ax=gca;
ax.YDir = 'normal';

size_img_px = size(im_disp);
size_img_cm = 100 .* (size_img_px ./ gamma_px);

xticks_desired = 0:10:size_img_cm(2);
yticks_desired = 0:10:size_img_cm(1);

xticks = (xticks_desired./100) .* gamma_px;
yticks = (yticks_desired./100) .* gamma_px;

xticks_label = num2str(xticks_desired');
yticks_label = num2str(yticks_desired');

ax = gca;
ax.YDir = 'normal';
ax.Visible = 'On';
ax.YTickLabel = yticks_label;
ax.YTick = yticks;
ax.XTickLabel = xticks_label;
ax.XTick = xticks;

xlabel('$x \, $(cm)')
ylabel('$z \, $(cm)')

colormap(cmap),
h2 = colorbar;
ylabel(h2,'$\overline{T}-T_i$ (K)','interpreter','latex')

caxis([-T_scale_span T_scale_span])

for ii = 1:length(processed_data.middle_position_z_px_notch)
    plot([765 765],[1 size_px(2)],'k','LineWidth',2,'MarkerSize',15)
    plot([783 783],[1 size_px(2)],'k--','LineWidth',2,'MarkerSize',15)
end
legend('Plot','Notch')

if strcmp(save_imgs,'on')
    export_fig([path_output 'vertical_profile_postion.png'],'-png')
end


%% Profile - middle position zoom

im_disp = scaled2rgb(processed_data.im_T_mean - T_int,cmap,[-T_scale_span T_scale_span]);

size_roi_cm = 5;
size_roi_px = ceil((size_roi_cm/100) * gamma_px);

roi = [1:size_roi_px; 774-0.5*size_roi_px:773+0.5*size_roi_px];
im_disp = im_disp(roi(1,:),roi(2,:),:);

fig = figure('Position',[100 100 1000 1000]);
imshow(im_disp,[]),
hold on,
ax=gca;
ax.YDir = 'normal';

size_img_px = size(im_disp);
size_img_cm = 100 .* (size_img_px ./ gamma_px);

colormap(cmap),
h2 = colorbar;
ylabel(h2,'$\overline{T}-T_i$ (K)','interpreter','latex')

caxis([-T_scale_span T_scale_span])

for ii = 1:(length(processed_data.middle_position_z_px_plot)-1)/2
    plot(processed_data.middle_position_x_px_plot - (774-0.5*size_roi_px),processed_data.middle_position_z_px_plot(ii),'kx','LineWidth',2,'MarkerSize',15)
    plot(processed_data.middle_position_x_px_notch - (774-0.5*size_roi_px),processed_data.middle_position_z_px_notch(ii),'kx','LineWidth',2,'MarkerSize',15)
    %     plot([765-(774-0.5*size_roi_px) 765-(774-0.5*size_roi_px)],[1 size_img_px(2)],'k','LineWidth',2,'MarkerSize',15)
    %     plot([783-(774-0.5*size_roi_px) 783-(774-0.5*size_roi_px)],[1 size_img_px(2)],'k--','LineWidth',2,'MarkerSize',15)
end
fig.Position = [100 100 1000 1000];
% legend('Plot','Notch')

if strcmp(save_imgs,'on')
    export_fig([path_output 'vertical_profile_postion_zoom_measures.png'],'-png')
end


%% Profile - middle position zoom

im_disp = scaled2rgb(processed_data.im_T_mean - T_int,cmap,[-T_scale_span T_scale_span]);

size_roi_cm = 5;
size_roi_px = ceil((size_roi_cm/100) * gamma_px);

roi = [1:size_roi_px; 774-0.5*size_roi_px:773+0.5*size_roi_px];
im_disp = im_disp(roi(1,:),roi(2,:),:);

fig = figure('Position',[100 100 1000 1000]);
imshow(im_disp,[]),
hold on,
ax=gca;
ax.YDir = 'normal';

size_img_px = size(im_disp);
size_img_cm = 100 .* (size_img_px ./ gamma_px);

colormap(cmap),
h2 = colorbar;
ylabel(h2,'$\overline{T}-T_i$ (K)','interpreter','latex')

caxis([-T_scale_span T_scale_span])

for ii = 1:length(processed_data.middle_position_z_px)
    %     plot(processed_data.middle_position_x_px,processed_data.middle_position_z_px(ii),'kx','LineWidth',2,'MarkerSize',15)
    plot([765-(774-0.5*size_roi_px) 765-(774-0.5*size_roi_px)],[1 size_img_px(2)],'k','LineWidth',2,'MarkerSize',15)
    plot([783-(774-0.5*size_roi_px) 783-(774-0.5*size_roi_px)],[1 size_img_px(2)],'k--','LineWidth',2,'MarkerSize',15)
end
fig.Position = [100 100 1000 1000];
% legend('Plot','Notch')

if strcmp(save_imgs,'on')
    export_fig([path_output 'vertical_profile_postion_zoom_measures.png'],'-png')
end



%% Profile - vertical plots

% z_vector = linspace(0,size_cm(1),size_px(1)-vertical_offset_px_top-vertical_offset_px_bot+1);
z_vector = linspace(0,size_cm(1),size_px(1));

figure('Position',[100 100 1000 1000]);
hold on,
% plot(z_vector,processed_data.im_T_mean(vertical_offset_px_bot:size_px(1)-vertical_offset_px_top,middle_position_x_px(1)),'LineWidth',2);
% plot(z_vector,processed_data.im_T_mean(vertical_offset_px_bot:size_px(1)-vertical_offset_px_top,middle_position_x_px(2)),'--','LineWidth',2);
plot(z_vector,processed_data.im_T_mean(:,middle_position_x_px(1)),'LineWidth',2);
plot(z_vector,processed_data.im_T_mean(:,middle_position_x_px(2)),'--','LineWidth',2);
axis([0 40.7 20 50])
ylabel('$\overline{T}$ ($^{\circ}$C)')
xlabel('z (cm)')
legend('Plot','Notch')

if strcmp(save_imgs,'on')
    export_fig([path_output 'vertical_profile_postion_plot_full.png'],'-png')
end


%% Profile - vertical plots

z_vector = linspace(0,size_cm(1),size_px(1));

figure('Position',[100 100 600 1200]);
[ha,~] = tight_subplot(2,1,[.1 .03],[.15 .08],[.18 .08]);

axes(ha(1));
hold on,
plot(processed_data.im_T_mean(:,middle_position_x_px(1)),'LineWidth',2);
plot(processed_data.im_T_mean(:,middle_position_x_px(2)),'--','LineWidth',2);
axis([0 20 18 50])
ax = gca;
ax.XTick = 0:5:20;
ax.XTickLabel = num2str((0:5:20)');
ax.YTick = 20:3:50;
ax.YTickLabel = num2str((20:3:50)');
text(5,23,'Rough','FontSize',25)
ylabel('$\overline{T}$ ($^{\circ}$C)')
xlabel('z (px)')
legend('Plot','Notch','location','North')

axes(ha(2));
hold on,
plot(processed_data.im_T_mean(:,middle_position_x_px(1)),'LineWidth',2);
% plot(processed_data.im_T_mean(:,middle_position_x_px(2)),'--','LineWidth',2);
axis([size_px(1)-20 size_px(1) 18 50])
ax = gca;
ax.XTick = size_px(1)-20:5:size_px(1);
ax.XTickLabel = num2str((size_px(1)-20:5:size_px(1))');
ax.YTick = 20:3:50;
ax.YTickLabel = num2str((20:3:50)');
text(size_px(1)-15,23,'Smooth','FontSize',25)
xlabel('z (px)')
ylabel('$\overline{T}$ ($^{\circ}$C)')


if strcmp(save_imgs,'on')
    export_fig([path_output 'vertical_profile_postion_plot_BL.png'],'-png')
end


%}

%% Generate movies

% inf_run = 4;

if strcmp(disp_imgs,'on')
    num_imgs = 500;
    imgs_processed = dir([paths.output_processed_images_run '*.mat']);
    
    for ii = 1:num_imgs%length(imgs_processed)
        img_name = [imgs_processed(ii).folder '/' imgs_processed(ii).name];
        
        load(img_name)
        
        im_T = im_T(sub_roi(1):sub_roi(2),sub_roi(3):sub_roi(4));
        im_T = im_T - processed_data.im_T_mean;
        
        im_disp = scaled2rgb(im_T,cmap,[-T_scale_span T_scale_span]);
        if ii == 1
            fig = figure('Visible','off','Position',[20 20 1300 1300]);
        end
        
        imshow(im_disp);
        ax=gca;
        ax.YDir = 'normal';
        
        size_img_px = size(im_disp);
        size_img_cm = 100 .* (size_img_px ./ gamma_px);
        
        xticks_desired = 0:10:size_img_cm(2);
        yticks_desired = 0:10:size_img_cm(1);
        
        xticks = (xticks_desired./100) .* gamma_px;
        yticks = (yticks_desired./100) .* gamma_px;
        
        xticks_label = num2str(xticks_desired');
        yticks_label = num2str(yticks_desired');
        
        ax = gca;
        ax.YDir = 'normal';
        ax.Visible = 'On';
        ax.YTickLabel = yticks_label;
        ax.YTick = yticks;
        ax.XTickLabel = xticks_label;
        ax.XTick = xticks;
        
        xlabel('$x \, $(cm)')
        ylabel('$z \, $(cm)')
        
        colormap(cmap),
        h2 = colorbar;
        caxis([-T_scale_span T_scale_span])
        ylabel(h2, '$T-<T> \,$(K)','interpreter','latex')
        
        if strcmp(save_imgs,'on')
            export_fig([path_output 'movie/' 'im_' num2str(ii,'%03.f') '.jpg'],'-transparent','-jpg')
        end
    end
end


%% Generate zoomed in movies

if strcmp(disp_imgs,'on')
    num_imgs = 500;
    size_profile_cm = 10;
    size_profile_px = ceil((size_profile_cm./100) .* gamma_px);
    
    
    imgs_processed = dir([paths.output_processed_images_run '*.mat']);
    
    for ii = 1:num_imgs
        img_name = [imgs_processed(ii).folder '/' imgs_processed(ii).name];
        load(img_name)
        
        im_T = im_T(sub_roi(1):sub_roi(2),sub_roi(3):sub_roi(4));
        im_T = im_T - processed_data.im_T_mean;
        
        im_disp = scaled2rgb(im_T,cmap,[-T_scale_span T_scale_span]);
        if ii == 1
            
            fig = figure('Visible','off','Position',[20 20 1300 1300]);
            [ha,~] = tight_subplot(2,2,[-.1 .02],[.05 .01],[.08 .1]);
            
            size_img_px = size(im_disp);
            size_img_cm = 100 .* (size_img_px ./ gamma_px);
            
            yticks_desired_bot = 0:1:size_profile_cm;
            yticks_bot = (yticks_desired_bot./100) .* gamma_px + 1;
            yticks_label_bot = num2str(yticks_desired_bot');
            
            yticks_desired_top = floor(size_img_cm(1) - size_profile_cm:1:size_img_cm(1));
            yticks_top = flip(((size_img_cm(1) - yticks_desired_top)./100) .* gamma_px + 1);
            yticks_label_top = num2str(yticks_desired_top');
            
            xticks_desired_left = 0:1:size_profile_cm;
            xticks_left = (xticks_desired_left./100) .* gamma_px + 1;
            xticks_label_left = num2str(xticks_desired_left');
            
            xticks_desired_right = floor(size_img_cm(2) - size_profile_cm:1:size_img_cm(2));
            xticks_right = ((xticks_desired_left - (size_img_cm(2) - max(xticks_desired_right)))./100) .* gamma_px + 1;
            xticks_label_right = num2str(xticks_desired_right');
        end
        
        axes(ha(1))
        imshow(im_disp(size_img_px(1)-size_profile_px-8:size_img_px(1)-8,1:1+size_profile_px,:));
        fig.Visible = 'off';
        
        ax = gca;
        ax.Visible = 'On';
        ax.YDir = 'normal';
        ax.XTickLabel = [];
        ax.XTick = xticks_left;
        ax.YTickLabel = yticks_label_top;
        ax.YTick = yticks_top;
        ylabel('$z \, $(cm)')
        
        colormap(cmap),
        caxis([-T_scale_span T_scale_span])
        
        
        axes(ha(2))
        imshow(im_disp(size_img_px(1)-size_profile_px-8:size_img_px(1)-8,size_img_px(2)-size_profile_px:size_img_px(2),:));
        originalSize = get(gca, 'Position');
        
        ax = gca;
        ax.Visible = 'On';
        ax.YDir = 'normal';
        ax.XTickLabel = [];
        ax.XTick = xticks_right;
        ax.YTickLabel = [];
        ax.YTick = yticks;
        
        colormap(cmap),
        h2 = colorbar;
        ylabel(h2, '$T-<T> \,$(K)','interpreter','latex')
        caxis([-T_scale_span T_scale_span])
        set(ha(2), 'Position', originalSize);
        
        
        axes(ha(3))
        imshow(im_disp(1:1+size_profile_px,1:1+size_profile_px,:));
        
        ax = gca;
        ax.Visible = 'On';
        ax.YDir = 'normal';
        ax.XTickLabel = xticks_label_left;
        ax.XTick = xticks_left;
        ax.YTickLabel = yticks_label_bot;
        ax.YTick = yticks_bot;
        xlabel('$x \, $(cm)')
        ylabel('$z \, $(cm)')
        
        colormap(cmap),
        caxis([-T_scale_span T_scale_span])
        
        
        axes(ha(4))
        imshow(im_disp(1:1+size_profile_px,size_img_px(2)-size_profile_px:size_img_px(2),:));
        originalSize = get(gca, 'Position');
        
        ax = gca;
        ax.Visible = 'On';
        ax.YDir = 'normal';
        ax.XTickLabel = xticks_label_right;
        ax.XTick = xticks_right;
        ax.YTickLabel = [];
        ax.YTick = yticks_bot;
        xlabel('$x \, $(cm)')
        
        colormap(cmap),
        h2 = colorbar;
        ylabel(h2, '$T-<T> \,$(K)','interpreter','latex')
        caxis([-T_scale_span T_scale_span])
        set(ha(4), 'Position', originalSize);
        
        if strcmp(save_imgs,'on')
            export_fig([path_output 'movie_zoom/' 'im_zoom_' num2str(ii,'%03.f') '.jpg'],'-transparent','-jpg')
        end
    end
end