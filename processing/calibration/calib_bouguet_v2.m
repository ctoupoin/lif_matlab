function [params_cam1,params_cam2] = calib_bouguet_v2(num_target_plans,imageSize,...
    calib_pattern_u,calib_pattern_d,c1_matched_u,c1_matched_d,c2_matched_u,c2_matched_d)

n_ima = 2.*num_target_plans;
nx = imageSize(2);
ny = imageSize(1);
% check_cond = 0;
two_focals_init = 0;
est_aspect_ratio = 0;
center_optim = 4;
est_dist = [1;1;1;1;0];


c1_matched = [c1_matched_d c1_matched_u];
c2_matched = [c2_matched_d c2_matched_u];

%% Calibration - cam0

% We are using the syntax of Jea-Yves Bouguet's calibration toolbox so that
% we can shamelessly use his functions as they are.

% imageSize = size(img0);

for ii = 1:num_target_plans
    p = c1_matched_d{ii};
    m = calib_pattern_d(:,:);
    
    eval(sprintf('x_%d = p(:,:)'';',ii));
    eval(sprintf('X_%d = m(:,:)'';',ii));
    eval(sprintf('X_%d(3,:) = 0;',ii));
    
    eval(sprintf('x_left_%d = x_%d'';',ii,ii));
    eval(sprintf('X_left_%d = X_%d'';',ii,ii));
end
for ii = 1:num_target_plans
    p = c1_matched_u{ii};
    m = calib_pattern_u(:,:);
    
    eval(sprintf('x_%d = p(:,:)'';',ii+num_target_plans));
    eval(sprintf('X_%d = m(:,:)'';',ii+num_target_plans));
    eval(sprintf('X_%d(3,:) = 0;',ii+num_target_plans));
    
    eval(sprintf('x_left_%d = x_%d'';',ii+num_target_plans,ii+num_target_plans));
    eval(sprintf('X_left_%d = X_%d'';',ii+num_target_plans,ii+num_target_plans));
end


go_calib_optim;

r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
dist_amount = 1;
fc_new = dist_amount * fc;

KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];

params_cam1.fc = fc;
params_cam1.cc = cc;
params_cam1.kc = kc;
params_cam1.KK_new = KK_new;

for ii = 1:n_ima
    eval(sprintf('[omc_left_%d,Tc_left_%d,~,~] = compute_extrinsic(x_%d,X_%d,fc,cc,kc,alpha_c);',ii,ii,ii,ii));
end

ext_calib
saving_calib_left


%% Calibration - cam1

% We are using the syntac of Jea-Yves Bouguet's calibration toolbox so that
% we can shamelessly use his functions as they are.

% imageSize = size(img1);

for ii = 1:num_target_plans
    p = c2_matched_u{ii};
    m = calib_pattern_u(:,:);
    
    eval(sprintf('x_%d = p(:,:)'';',ii));
    eval(sprintf('X_%d = m(:,:)'';',ii));
    eval(sprintf('X_%d(3,:) = 0;',ii));
    
    eval(sprintf('x_right_%d = x_%d'';',ii,ii));
    eval(sprintf('X_right_%d = X_%d'';',ii,ii));
end
for ii = 1:num_target_plans
    p = c2_matched_d{ii};
    m = calib_pattern_d(:,:);
    
    eval(sprintf('x_%d = p(:,:)'';',ii+num_target_plans));
    eval(sprintf('X_%d = m(:,:)'';',ii+num_target_plans));
    eval(sprintf('X_%d(3,:) = 0;',ii+num_target_plans));
    
    eval(sprintf('x_right_%d = x_%d'';',ii+num_target_plans,ii+num_target_plans));
    eval(sprintf('X_right_%d = X_%d'';',ii+num_target_plans,ii+num_target_plans));
end

go_calib_optim;

r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
dist_amount = 1;
fc_new = dist_amount * fc;

KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];

params_cam2.fc = fc;
params_cam2.cc = cc;
params_cam2.kc = kc;
params_cam2.KK_new = KK_new;


for ii = 1:n_ima
    eval(sprintf('[omc_right_%d,Tc_right_%d,~,~] = compute_extrinsic(x_%d,X_%d,fc,cc,kc,alpha_c);',ii,ii,ii,ii));
end

ext_calib
saving_calib_right


%% Calibration - stereo

est_dist_left = [1;0;0;0;0];
est_aspect_ratio_left = 0;
est_alpha_left = 0;
est_fc_left = 0;
alpha_c_left = 0;


fc_left = params_cam1.fc;
cc_left = params_cam1.cc;
kc_left = params_cam1.kc;
KK_new_left = params_cam1.KK_new;


for ii = 1:n_ima
    p = c2_matched{ii};
    m = calib_pattern(:,:);
    eval(sprintf('x_%d = p(:,:)'';',ii));
    eval(sprintf('X_%d = m(:,:)'';',ii));
    eval(sprintf('X_%d(3,:) = 0;',ii));
end

center_optim_right = 1;
est_dist_right = [1;0;0;0;0];
est_aspect_ratio_right = 0;
est_alpha_right = 0;
est_fc_right = 0;
alpha_c_right = 0;


fc_right = params_cam1.fc;
cc_right = params_cam1.cc;
kc_right = params_cam1.kc;
KK_new_right = params_cam1.KK_new;

active_images_left = ones(1,10);
active_images_right = ones(1,10);

calib_stereo


params_cam1.fc = fc_left;
params_cam1.cc = cc_left;
params_cam1.kc = kc_left;
params_cam1.KK_new = KK_left;

params_cam2.fc = fc_right;
params_cam2.cc = cc_right;
params_cam2.kc = kc_right;
params_cam2.KK_new = KK_right;


%% Calibration - individual cameras again
% We restart the calibration for the individual cameras in case the stereo
% calibration discarded some bad views. This will improve the general
% calibration result.

for ii = 1:n_ima
    p = c1_matched{ii};
    m = calib_pattern(:,:);
    
    eval(sprintf('x_%d = p(:,:)'';',ii));
    eval(sprintf('X_%d = m(:,:)'';',ii));
    eval(sprintf('X_%d(3,:) = 0;',ii));
    
    eval(sprintf('x_left_%d = x_%d'';',ii,ii));
    eval(sprintf('X_left_%d = X_%d'';',ii,ii));
end

go_calib_optim;

r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
dist_amount = 1;
fc_new = dist_amount * fc;

KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];

% params_cam1.fc = fc;
% params_cam1.cc = cc;
% params_cam1.kc = kc;
% params_cam1.KK_new = KK_new;

for ii = 1:n_ima
    eval(sprintf('[omc_left_%d,Tc_left_%d,~,~] = compute_extrinsic(x_%d,X_%d,fc,cc,kc,alpha_c);',ii,ii,ii,ii));
end

ext_calib
saving_calib_left


for ii = 1:n_ima
    p = c2_matched{ii};
    m = calib_pattern(:,:);
    
    eval(sprintf('x_%d = p(:,:)'';',ii));
    eval(sprintf('X_%d = m(:,:)'';',ii));
    eval(sprintf('X_%d(3,:) = 0;',ii));
    
    eval(sprintf('x_right_%d = x_%d'';',ii,ii));
    eval(sprintf('X_right_%d = X_%d'';',ii,ii));
end

go_calib_optim;

r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
dist_amount = 1;
fc_new = dist_amount * fc;

KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];

% params_cam2.fc = fc;
% params_cam2.cc = cc;
% params_cam2.kc = kc;
% params_cam2.KK_new = KK_new;


for ii = 1:n_ima
    eval(sprintf('[omc_right_%d,Tc_right_%d,~,~] = compute_extrinsic(x_%d,X_%d,fc,cc,kc,alpha_c);',ii,ii,ii,ii));
end

ext_calib
saving_calib_right
