%                          Camera calibration
%
% Calibrates the cameras, and uses the calibration target to generate the
% transfer matrix allowing to go for the fov of one camera to the other, as
% well as the one allowing to go from the image plane to world coordinates.
%
%--------------------------------------------------------------------------
%
%   *Author :* Cl\'ement Toupoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changelog :
%   - v12 : updated the script to work with new calibration images, where
%   two calibration grids are visible per plane
%   - v11 : changed notations from cam0/cam1 to cam1/cam2
%   - v8 : subtracted irregular background from calibration images to
%   improve detection
%   - v7 : after un-distrotion of the images, we now only match grid points
%   for the plane corresponding to the laser sheet (i.e. the first one), as
%   it's the only one needed. Added the computation of the geometric
%   transform to go from cam 0 to 1 and inversely.
%   - v6 : adaptation of the code for wirlong with PCO
%   - v5 : added computation of the transform matrix without removing
%   camera distortion
%   - v4 : finally, changed to using JYB's calibration toolbox
%   - v3 : revamp / simplification of the code
%   - v2 : implemented a calibration method using MATLAB's own tools, from
% the computer vision toolbox

clearvars
close all

% Here we just get rid of some pesky -and useless- warning messages
warning off images:initSize:adjustingMag;
warning off images:imfindcircles:warnForLargeRadiusRange

if ~exist('paths','var')
    run('../main_double.m')
end


%% User inputs

% Calibration target : cicles  3mm in diameter, big circles 6 mm in
% diameter, total lenght 390 mm
% The calibration target is made of two different grids of points. The
% calibration routine can only process grid of evenly spaced points.
% Therefore the two grid have to be processed independently, even though
% they are visible on a same image. In the following, fields related to the
% bottom grid are appended a _d (down) and fields related to the upper
% grid, _u (up).

format_imgs = 'tif';
d_small = 2;                    % diameter of the small dots
d_big = 6;                      % diameter of the big dots -unused-
d_pts = 5;                      % distance between 2 dots

rows_d = 21;                    % number of rows of dots to be processed
rows_u = 21;                    % number of rows of dots to be processed
columns = 55;                   % number of columns of rows to be processed

dist_mm = d_pts * ( columns - 1);  % distance between the two dots at top left and top right corners of the grid
d_px = [5 15];                  % estimation of the pixel diameters of the dots. Used in the imfindcircles functions, and needs to be adjusted when the calibration pattern is changed or the cameras are moved.


%% Paths

% Path used to read and save information
paths.target = [paths.data 'calibration_target/'];
paths.mask1 = [paths.misc 'mask1.mat'];
paths.mask2 = [paths.misc 'mask2.mat'];
paths.params_cam1 = [paths.misc 'params_cam1.mat'];
paths.params_cam2 = [paths.misc 'params_cam2.mat'];
paths.c1_matched = [paths.misc 'c1_matched.mat'];
paths.r1_matched = [paths.misc 'r1_matched.mat'];
paths.c2_matched = [paths.misc 'c2_matched.mat'];
paths.r2_matched = [paths.misc 'r2_matched.mat'];
paths.c1_matched_ud = [paths.misc 'c1_matched_ud.mat'];
paths.r1_matched_ud = [paths.misc 'r1_matched_ud.mat'];
paths.c2_matched_ud = [paths.misc 'c2_matched_ud.mat'];
paths.r2_matched_ud = [paths.misc 'r2_matched_ud.mat'];
paths.target_infos = [paths.misc 'target_infos.mat'];

data_files = dir(paths.target);
count = 0;
for ii = 1:length(data_files)
    if length(data_files(ii).name) > 10 && strcmp(data_files(ii).name(1:11),'target_plan')
        count = count + 1;
        paths.target_plan(count,:) = [data_files(ii).folder '/' data_files(ii).name];
    end
end
paths.no_target = [data_files(ii).folder '/no_target/'];

addpath(genpath(paths.ext_functions));


%% Load  and average images

% When several images of a static calibration pattern are available, they
% are averaged in order to reduce noise.

num_target_plans = size(paths.target_plan,1);

pattern1 = cell(1,num_target_plans);
pattern2 = cell(1,num_target_plans);

for ii = 1:size(paths.target_plan,1)
    imgs = dir([paths.target_plan(ii,:) '/' '*.' format_imgs]);
    
    ind_cam1 = ismember(vertcat(imgs.name),'cam1');
    ind_cam1 = ind_cam1(:,4) == true;
    ind_cam1 = ind_cam1';
    ind_cam2 = ~ind_cam1;
    
    if sum(ind_cam1) == 0
        ind_cam1 = ismember(vertcat(imgs.name),'C1');
        ind_cam1 = ind_cam1(:,17) == true;
        ind_cam1 = ind_cam1';
        ind_cam2 = ~ind_cam1;
    end
    
    ind_cam1 = find(ind_cam1);
    ind_cam2 = find(ind_cam2);

    char_length = length([imgs(end).folder '/' imgs(end).name]);
    imgs_name1 = repmat('a',length(ind_cam1),char_length);
    imgs_name2 = repmat('a',length(ind_cam2),char_length);
    
    count1 = 0;
    for jj = ind_cam1
        count1 = count1 + 1;
        imgs_name1(jj,:) = [imgs(jj).folder '/' imgs(jj).name];
        
        img1 = imrotate(imread(imgs_name1(count1,:)),-90);
        
        if count1 == 1
            sum_img1 = img1;
        else
            sum_img1 = count1 / (count1 + 1) .* sum_img1 + 1 ./ (count1 + 1) .* img1;
        end
    end
    
    count2 = 0;
    for jj = ind_cam2
        count2 = count2 + 1;
        imgs_name2(count2,:) = [imgs(jj).folder '/' imgs(jj).name];
        
        img2 = imrotate(imread(imgs_name2(count2,:)),90);
        
        if count2 == 1
            sum_img2 = img2;
        else
            sum_img2 = count2 / (count2 + 1) .* sum_img2 + 1 ./ (count2 + 1) .* img2;
        end
    end
    
    pattern1{ii} = sum_img1;
    pattern2{ii} = sum_img2;
end

num_img1 = length(pattern1);
num_img2 = length(pattern2);


%% Mask creation

% The calibration trarget is a dot pattern. We intend on using imfindcircle
% to detect the position of the dots on the different calibration images.
% However, the dot size is a bit small for imfindcircles to detect them
% easily. In order to facilitate imfindcircles's work, we create masks of
% the calibration targets, leaving only the dot pattern visible.
% Warning, this step is tedious.

% For each image, the user has to draw a rectangle encompassing the dot
% pattern. It is best if the rectangle closely fits the pattern. The user
% can do this with a pointer dirctly on the image. Each mouse click creates
% a new point. An eroneous point can be deletedby rightè=-clicking it and
% selecting the option from the drop-down menu. Once the contour is closed,
% the user needs to right-click it and select the "create mask" option.

% This stp has to be done for each calibration target, that is
% 2*num_target_plans times, because there are two targets per plane. The
% user is encouraged to save the masks so that they aare not lost should
% the procedure be interrupted.

if ~exist(paths.mask1,'file')
    mask1_u = cell(1,num_target_plans);
    for ii = 1:num_target_plans
        mask1_u{ii} = roipoly(imadjust(pattern1{ii}));
    end
    
    mask1_d = cell(1,num_target_plans);
    for ii = 1:num_target_plans
        mask1_d{ii} = roipoly(imadjust(pattern1{ii}));
    end
    
    mask1.u = mask1_u;
    mask1.d = mask1_d;
    
    save('mask1',paths.mask1)
else
    load(paths.mask1)
end

if ~exist(paths.mask2,'file')
    mask2_u = cell(1,num_target_plans);
    for ii = 1:num_target_plans
        mask2_u{ii} = roipoly(imadjust(pattern2{ii}));
    end
    
    mask2_d = cell(1,num_target_plans);
    for ii = 1:num_target_plans
        mask2_d{ii} = roipoly(imadjust(pattern2{ii}));
    end
    
    mask2.u = mask2_u;
    mask2.d = mask2_d;
    
    save('mask2',paths.mask2)
else
    load(paths.mask2)
end


%% Calibration pattern detection
se = strel('disk',15);

% First detection of the calibration target, using imfindcircles.

calib_pattern1_u = cell(1,num_target_plans);
calib_pattern1_d = cell(1,num_target_plans);

c1_u = cell(1,num_target_plans);
r1_u = cell(1,num_target_plans);

c1_d = cell(1,num_target_plans);
r1_d = cell(1,num_target_plans);
for ii = 1:num_target_plans
    sens_thresh = 0.85;
    
    calib_pattern1_u{ii} = pattern1{ii} .* uint16(mask1.u{ii}); % mask is applied
    calib_pattern1_u{ii} = double(imadjust(calib_pattern1_u{ii})); % contrast enhancement
    
    irregular_bg = imclose(calib_pattern1_u{ii},se); % the dots are morphologically "erased" from the image to create a background
    calib_pattern1_u{ii} = abs(calib_pattern1_u{ii}-irregular_bg); % background subtraction to enhance contrast
    calib_pattern1_u{ii}(calib_pattern1_u{ii}<0.17*max(max(calib_pattern1_u{ii}))) = 0; % dark grey areas are set to black
    calib_pattern1_u{ii} = imbinarize(calib_pattern1_u{ii});

    
    [c1_u{ii},r1_u{ii}] = imfindcircles(calib_pattern1_u{ii},d_px,'ObjectPolarity',...
        'bright','Sensitivity',sens_thresh);
    
    while length(r1_u{ii}) < rows_u * columns && sens_thresh < 0.99
        % the detectin threshold of imfindcircles is increased as long as
        % the correct number of detected dots has not been found
        sens_thresh = sens_thresh + 0.01;
        [c1_u{ii},r1_u{ii}] = imfindcircles(calib_pattern1_u{ii},d_px,'ObjectPolarity',...
            'bright','Sensitivity',sens_thresh);
    end
    
    % Then the trehsold is reset and the same operation is carried out for
    % the bottom dot pattern
    sens_thresh = 0.85;
    
    calib_pattern1_d{ii} = pattern1{ii} .* uint16(mask1.d{ii});
    calib_pattern1_d{ii} = double(imadjust(calib_pattern1_d{ii}));
    
    irregular_bg = imclose(calib_pattern1_d{ii},se);
    calib_pattern1_d{ii} = abs(calib_pattern1_d{ii}-irregular_bg);
    calib_pattern1_d{ii}(calib_pattern1_d{ii}<0.17*max(max(calib_pattern1_d{ii}))) = 0;
    calib_pattern1_d{ii} = imbinarize(calib_pattern1_d{ii});

    
    [c1_d{ii},r1_d{ii}] = imfindcircles(calib_pattern1_d{ii},d_px,'ObjectPolarity',...
        'bright','Sensitivity',sens_thresh);
    
    while length(r1_d{ii}) < rows_d * columns && sens_thresh < 0.99
        sens_thresh = sens_thresh + 0.01;
        [c1_d{ii},r1_d{ii}] = imfindcircles(calib_pattern1_d{ii},d_px,'ObjectPolarity',...
            'bright','Sensitivity',sens_thresh);
    end
end

%% Uncomment and run to see the detected dots
% figure('Position',[100 100 1000 1000])
% for ii = 1:num_target_plans
%    imshow(pattern1{ii},[])
%    viscircles(c1_u{ii},r1_u{ii})
%    pause
% end
% close

%% Preliminary dot detection, second camera
calib_pattern2_u = cell(1,num_target_plans);
calib_pattern2_d = cell(1,num_target_plans);

c2_u = cell(1,num_target_plans);
r2_u = cell(1,num_target_plans);

c2_d = cell(1,num_target_plans);
r2_d = cell(1,num_target_plans);
for ii = 1:num_target_plans
    sens_thresh = 0.85;
    
    calib_pattern2_u{ii} = pattern2{ii} .* uint16(mask2.u{ii});
    calib_pattern2_u{ii} = double(imadjust(calib_pattern2_u{ii}));
    
    irregular_bg = imclose(calib_pattern2_u{ii},se);
    calib_pattern2_u{ii} = abs(calib_pattern2_u{ii}-irregular_bg);
    calib_pattern2_u{ii}(calib_pattern2_u{ii}<0.17*max(max(calib_pattern2_u{ii}))) = 0;
    calib_pattern2_u{ii} = imbinarize(calib_pattern2_u{ii});

    
    [c2_u{ii},r2_u{ii}] = imfindcircles(calib_pattern2_u{ii},d_px,'ObjectPolarity',...
        'bright','Sensitivity',sens_thresh);
    
    while length(r2_u{ii}) < rows_u * columns && sens_thresh < 0.99
        sens_thresh = sens_thresh + 0.01;
        [c2_u{ii},r2_u{ii}] = imfindcircles(calib_pattern2_u{ii},d_px,'ObjectPolarity',...
            'bright','Sensitivity',sens_thresh);
    end
    
    sens_thresh = 0.85;
    
    calib_pattern2_d{ii} = pattern2{ii} .* uint16(mask2.d{ii});
    calib_pattern2_d{ii} = double(imadjust(calib_pattern2_d{ii}));
    
    irregular_bg = imclose(calib_pattern2_d{ii},se);
    calib_pattern2_d{ii} = abs(calib_pattern2_d{ii}-irregular_bg);
    calib_pattern2_d{ii}(calib_pattern2_d{ii}<0.17*max(max(calib_pattern2_d{ii}))) = 0;
    calib_pattern2_d{ii} = imbinarize(calib_pattern2_d{ii});

    
    [c2_d{ii},r2_d{ii}] = imfindcircles(calib_pattern2_d{ii},d_px,'ObjectPolarity',...
        'bright','Sensitivity',sens_thresh);
    
    while length(r2_d{ii}) < rows_d * columns && sens_thresh < 0.99
        sens_thresh = sens_thresh + 0.01;
        [c2_d{ii},r2_d{ii}] = imfindcircles(calib_pattern2_d{ii},d_px,'ObjectPolarity',...
            'bright','Sensitivity',sens_thresh);
    end
end

%% Uncomment and run to see the detected dots
% figure('Position',[100 100 1000 1000])
% for ii = 1:num_target_plans
%    imshow(pattern2{ii},[])
%    viscircles(c2_d{ii},r2_d{ii})
%    pause
% end
% close

%% Synthetic grid : preliminary step

% In order to perform stereo camera calibration, we need to match the
% position of the detected dots between the two cameras. The issue here is
% that imfindcircles does not sort the dots it dectects. Therefore the dot
% position do not match between the two cameras.
% We solve this issue by creating a "synthetic" grid in which the position
% of the dots is known. We match the detected dots to the ones of the grid
% by computing the distance between them, and taking the minimum distance
% between synthetic and detected dots.
% Because the calibration target is not always orthogonal to the camera
% axis, the dot grid does not look straight: it has an apparent
% inclination. Therefore, in order to properly match the detected and
% synthetic dots, we need to correct the effect of that apparent
% inclination.


% In the following, the user is asked to click on the position of three
% dots, in order:
%   - upper left corner
%   - upper right corner
%   - middle dot

if ~exist(paths.target_infos,'file')
    %% Camera 1 - top
    c1_matched_u = cell(1,num_target_plans);
    r1_matched_u = cell(1,num_target_plans);
    
    [d_px1_u,O_1_u,inclination1_u] = get_target_infos(calib_pattern1_u,c1_u,r1_u);
    gamma1_u = d_px1_u ./ dist_mm; %[px.mm^{-1}]
    
    %% Camera 1 - bot
    c1_matched_d = cell(1,num_target_plans);
    r1_matched_d = cell(1,num_target_plans);
    
    [d_px1_d,O_1_d,inclination1_d] = get_target_infos(calib_pattern1_d,c1_d,r1_d);
    gamma1_d = d_px1_d ./ dist_mm; %[px.mm^{-1}]
    
    
    gamma_1 = mean([gamma1_u gamma1_d]);
    
    %% Camera 2 - top
    c2_matched_u = cell(1,num_target_plans);
    r2_matched_u = cell(1,num_target_plans);
    
    [d_px2_u,O_2_u,inclination2_u] = get_target_infos(calib_pattern2_u,c2_u,r2_u);
    gamma2_u = d_px2_u ./ dist_mm; %[px.mm^{-1}]
    
    %% Camera 2 - bot
    c2_matched_d = cell(1,num_target_plans);
    r2_matched_d = cell(1,num_target_plans);
    
    [d_px2_d,O_2_d,inclination2_d] = get_target_infos(calib_pattern2_d,c2_d,r2_d);
    gamma2_d = d_px2_d ./ dist_mm; %[px.mm^{-1}]
    
    
    gamma_2 = mean([gamma2_u gamma2_d]);
    
    
    %% Data
    target_infos.d_px1_u = d_px1_u;
    target_infos.O_1_u = O_1_u;
    target_infos.inclination1_u = inclination1_u;
    target_infos.gamma1_u = gamma1_u;
    
    target_infos.d_px1_d = d_px1_d;
    target_infos.O_1_d = O_1_d;
    target_infos.inclination1_d = inclination1_d;
    target_infos.gamma1_d = gamma1_d;
    
    target_infos.d_px2_u = d_px2_u;
    target_infos.O_2_u = O_2_u;
    target_infos.inclination2_u = inclination2_u;
    target_infos.gamma2_u = gamma2_u;
    
    target_infos.d_px2_d = d_px2_d;
    target_infos.O_2_d = O_2_d;
    target_infos.inclination2_d = inclination2_d;
    target_infos.gamma2_d = gamma2_d;
    
    save(paths.target_infos,'target_infos')
else
    load(paths.target_infos);
    
    d_px1_u = target_infos.d_px1_u;
    O_1_u = target_infos.O_1_u;
    inclination1_u = target_infos.inclination1_u;
    gamma1_u = target_infos.gamma1_u;
    
    d_px1_d = target_infos.d_px1_d;
    O_1_d = target_infos.O_1_d;
    inclination1_d = target_infos.inclination1_d;
    gamma1_d = target_infos.gamma1_d;
    
    d_px2_u = target_infos.d_px2_u;
    O_2_u = target_infos.O_2_u;
    inclination2_u = target_infos.inclination2_u;
    gamma2_u = target_infos.gamma2_u;
    
    d_px2_d = target_infos.d_px2_d;
    O_2_d = target_infos.O_2_d;
    inclination2_d = target_infos.inclination2_d;
    gamma2_d = target_infos.gamma2_d;
end


%% Synthetic grid : creation

% Pattern generation
calib_pattern_u = zeros(columns*rows_u,2);
calib_pattern_d = zeros(columns*rows_d,2);
count = 0;
for ii = 1:columns
    for jj = 1:rows_u
        count = count + 1;
        
        calib_pattern_u(count,1) = (ii - 1) * d_pts;
        calib_pattern_u(count,2) = (jj - 1) * d_pts;
    end
end
count = 0;
for ii = 1:columns
    for jj = 1:rows_d
        count = count + 1;
        
        calib_pattern_d(count,1) = (ii - 1) * d_pts;
        calib_pattern_d(count,2) = (jj - 1) * d_pts;
    end
end

% A "perfect", synthetic grid is generated for matching the points, and
% excluding artefacts from imfindcircles

% Camera 0
c1_synth_u = synthetic_grid_creation(inclination1_u,d_pts,gamma1_u,O_1_u,num_target_plans,columns,rows_u);
c1_synth_d = synthetic_grid_creation(inclination1_d,d_pts,gamma1_d,O_1_d,num_target_plans,columns,rows_d);

% Camera 1
c2_synth_u = synthetic_grid_creation(inclination2_u,d_pts,gamma2_u,O_2_u,num_target_plans,columns,rows_u);
c2_synth_d = synthetic_grid_creation(inclination2_d,d_pts,gamma2_d,O_2_d,num_target_plans,columns,rows_d);


%% Points matching

% Camera 1
if ~exist(paths.c1_matched,'file')
    [c1_matched_u,r1_matched_u] = points_matching(c1_u,r1_u,c1_synth_u,c1_matched_u,r1_matched_u,num_target_plans);
    [c1_matched_d,r1_matched_d] = points_matching(c1_d,r1_d,c1_synth_d,c1_matched_d,r1_matched_d,num_target_plans);
    
    c1_matched.u = c1_matched_u;
    c1_matched.d = c1_matched_d;
    
    r1_matched.u = r1_matched_u;
    r1_matched.d = r1_matched_d;
else
    load(paths.c1_matched)
    load(paths.r1_matched)
    
    c1_matched_u = c1_matched.u;
    c1_matched_d = c1_matched.d;
end


% Camera 2
if ~exist(paths.c2_matched,'file')
    [c2_matched_u,r2_matched_u] = points_matching(c2_u,r2_u,c2_synth_u,c2_matched_u,r2_matched_u,num_target_plans);
    [c2_matched_d,r2_matched_d] = points_matching(c2_d,r2_d,c2_synth_d,c2_matched_d,r2_matched_d,num_target_plans);
    
    c2_matched.u = c2_matched_u;
    c2_matched.d = c2_matched_d;
    
    r2_matched.u = r2_matched_u;
    r2_matched.d = r2_matched_d;
else
    load(paths.c2_matched)
    load(paths.r2_matched)
    
    c2_matched_u = c2_matched.u;
    c2_matched_d = c2_matched.d;
end


%% Uncomment and run to see the detected dots

% figure('Position',[100 100 1000 1000])
% for ii = 1:num_target_plans
%    imshow(pattern2{ii},[])
%    for jj = 1:length(c2_matched_d{ii})
%         viscircles(c2_matched_d{ii}(jj,:),r2_matched_d{ii}(jj))
%         pause(0.1)
%    end
%    pause
% end
% close

%% Creation of the transform matrix with still distorted images
% 
%   Od, unused
%
% c1_matched = [c1_matched_u{1}; c1_matched_d{1}];
% c2_matched = [c2_matched_u{1}; c2_matched_d{1}];
% 
% tform_1to2_d = fitgeotrans(c1_matched,c2_matched,'polynomial',4);
% tform_2to1_d = fitgeotrans(c2_matched,c1_matched,'polynomial',4);
% 
% if ~exist(paths.tform_1to2_d,'file')
%     save(paths.tform_1to2_d,'tform_1to2_d')
%     save(paths.tform_2to1_d,'tform_2to1_d')
% end


%% Calibration
% This is the main claibration step, doine using an adaptation of JY
% Bouguet's MATLAB calibraiton toolbox.
%   http://www.vision.caltech.edu/bouguetj/

% The toolbox was adapted to use a dot pattern -instead of a checkerboard-,
% based on the work of:
%   http://george-vogiatzis.org/calib/

% The routine is ran once for each camera, then the stereo calibration is
% performed. Because it is more accurate, we keep the results from stereo
% calibration.
% The toolbox will display the position of the calibration targets in the
% camera frame, which is a good step to check that the is no aberrant
% result. It also estimates the accuracy of the calibration in pixels. As
% of the writing of this comment, an accuracy of 0.29px for Cam1 and 0.21px
% for Cam2 has been achieved.

imageSize = size(img1);
c1_matched = [c1_matched_u c1_matched_d];
c2_matched = [c2_matched_u c2_matched_d];
calib_pattern = calib_pattern_d;

num_target_plans = length(c1_matched_d);

[params_cam1,params_cam2] = calib_bouguet(2.*num_target_plans,imageSize,calib_pattern,c1_matched,c2_matched);

if ~exist(paths.params_cam1,'file')
    save(paths.params_cam1,'params_cam1')
end

if ~exist(paths.params_cam2,'file')
    save(paths.params_cam2,'params_cam2')
end


%% Undistortion of the images

% Now that we have properly calbrated the cameras, we can remove the
% distottions from the images. We create _ud (un-distoreted) varaibales.

% Here, defining both calib_pattern_ud as cells might seem pointless, but
% it allows us to keep using the same functions as before for determining
% the target pattern orientation.

pattern1_ud{1} = rect(double(pattern1{1}),eye(3),...
    params_cam1.fc,params_cam1.cc,params_cam1.kc,params_cam1.KK_new);

pattern2_ud{1} = rect(double(pattern2{1}),eye(3),...
    params_cam2.fc,params_cam2.cc,params_cam2.kc,params_cam2.KK_new);


calib_pattern1_ud_u{1} = rect(double(calib_pattern1_u{1}),eye(3),...
    params_cam1.fc,params_cam1.cc,params_cam1.kc,params_cam1.KK_new);

calib_pattern2_ud_u{1} = rect(double(calib_pattern2_u{1}),eye(3),...
    params_cam2.fc,params_cam2.cc,params_cam2.kc,params_cam2.KK_new);


calib_pattern1_ud_d{1} = rect(double(calib_pattern1_d{1}),eye(3),...
    params_cam1.fc,params_cam1.cc,params_cam1.kc,params_cam1.KK_new);

calib_pattern2_ud_d{1} = rect(double(calib_pattern2_d{1}),eye(3),...
    params_cam2.fc,params_cam2.cc,params_cam2.kc,params_cam2.KK_new);



%% Detection of the pattern on undistorted images

% We now perform the same steps as before on undistorted images. Now, the
% goal is to compute the transform matyric allowing us to match the fields
% of the two cameras. We therefore need to detect the calibration pattern
% in the plane of the laser sheet, and match the position of the dots.
% The calibration targets 1 are the ones in the plane of the laser sheet.

sens_thresh = 0.85;
[c1_ud_d{1},r1_ud_d{1}] = imfindcircles(calib_pattern1_ud_d{1},d_px,'ObjectPolarity',...
    'bright','Sensitivity',sens_thresh);
while length(r1_ud_d{1}) < rows_d * columns
    sens_thresh = sens_thresh + 0.01;
    [c1_ud_d{1},r1_ud_d{1}] = imfindcircles(calib_pattern1_ud_d{1},d_px,'ObjectPolarity',...
        'bright','Sensitivity',sens_thresh);
end

sens_thresh = 0.85;
[c1_ud_u{1},r1_ud_u{1}] = imfindcircles(calib_pattern1_ud_u{1},d_px,'ObjectPolarity',...
    'bright','Sensitivity',sens_thresh);
while length(r1_ud_u{1}) < rows_u * columns
    sens_thresh = sens_thresh + 0.01;
    [c1_ud_u{1},r1_ud_u{1}] = imfindcircles(calib_pattern1_ud_u{1},d_px,'ObjectPolarity',...
        'bright','Sensitivity',sens_thresh);
end


sens_thresh = 0.85;
[c2_ud_d{1},r2_ud_d{1}] = imfindcircles(calib_pattern2_ud_d{1},d_px,'ObjectPolarity',...
    'bright','Sensitivity',sens_thresh);
while length(r2_ud_d{1}) < rows_d * columns
    sens_thresh = sens_thresh + 0.01;
    [c2_ud_d{1},r2_ud_d{1}] = imfindcircles(calib_pattern2_ud_d{1},d_px,'ObjectPolarity',...
        'bright','Sensitivity',sens_thresh);
end

sens_thresh = 0.85;
[c2_ud_u{1},r2_ud_u{1}] = imfindcircles(calib_pattern2_ud_u{1},d_px,'ObjectPolarity',...
    'bright','Sensitivity',sens_thresh);
while length(r2_ud_u{1}) < rows_u * columns
    sens_thresh = sens_thresh + 0.01;
    [c2_ud_u{1},r2_ud_u{1}] = imfindcircles(calib_pattern2_ud_u{1},d_px,'ObjectPolarity',...
        'bright','Sensitivity',sens_thresh);
end


% Get target infos - Camera 0
[d_px1_ud_u,O_1_ud_u,inclination1_ud_u] = get_target_infos(calib_pattern1_ud_u,c1_ud_u,r1_ud_u);
gamma1_ud_u = d_px1_ud_u ./ dist_mm; %[px.mm^{-1}]
[d_px1_ud_d,O_1_ud_d,inclination1_ud_d] = get_target_infos(calib_pattern1_ud_d,c1_ud_d,r1_ud_d);
gamma1_ud_d = d_px1_ud_d ./ dist_mm; %[px.mm^{-1}]


% Get target infos - Camera 1
[d_px2_ud_u,O_2_ud_u,inclination2_ud_u] = get_target_infos(calib_pattern2_ud_u,c2_ud_u,r2_ud_u);
gamma2_ud_u = d_px2_ud_u ./ dist_mm; %[px.mm^{-1}]
[d_px2_ud_d,O_2_ud_d,inclination2_ud_d] = get_target_infos(calib_pattern2_ud_d,c2_ud_d,r2_ud_d);
gamma2_ud_d = d_px2_ud_d ./ dist_mm; %[px.mm^{-1}]


% Synthetic grid - Camera 0
c1_synth_ud_u = synthetic_grid_creation(inclination1_ud_u,d_pts,gamma1_ud_u,O_1_ud_u,1,columns,rows_u);
c1_synth_ud_d = synthetic_grid_creation(inclination1_ud_d,d_pts,gamma1_ud_d,O_1_ud_d,1,columns,rows_d);

% Synthetic grid - Camera 1
c2_synth_ud_u = synthetic_grid_creation(inclination2_ud_u,d_pts,gamma2_ud_u,O_2_ud_u,1,columns,rows_u);
c2_synth_ud_d = synthetic_grid_creation(inclination2_ud_d,d_pts,gamma2_ud_d,O_2_ud_d,1,columns,rows_d);



%% Pattern detection on the undistorted image
c1_matched_ud_u = cell(size(c1_ud_u));
r1_matched_ud_u = cell(size(r1_ud_u));
c1_matched_ud_d = cell(size(c1_ud_d));
r1_matched_ud_d = cell(size(r1_ud_d));
if ~exist(paths.c1_matched_ud,'file')
    [c1_matched_ud_u,r1_matched_ud_u] = points_matching(c1_ud_u,r1_ud_u,c1_synth_ud_u,c1_matched_ud_u,r1_matched_ud_u,1);
    [c1_matched_ud_d,r1_matched_ud_d] = points_matching(c1_ud_d,r1_ud_d,c1_synth_ud_d,c1_matched_ud_d,r1_matched_ud_d,1);
    
    c1_matched_ud.u = c1_matched_ud_u;
    c1_matched_ud.d = c1_matched_ud_d;
    
    r1_matched_ud.u = r1_matched_ud_u;
    r1_matched_ud.d = r1_matched_ud_d;
else
    load(paths.c1_matched_ud)
    load(paths.r1_matched_ud)
    
    c1_matched_ud_u = c1_matched_ud.u;
    c1_matched_ud_d = c1_matched_ud.d;
    
    r1_matched_ud_u = r1_matched_ud.u;
    r1_matched_ud_d = r1_matched_ud.d;
end

c2_matched_ud_u = cell(size(c2_ud_u));
c2_matched_ud_d = cell(size(c2_ud_d));
r2_matched_ud_u = cell(size(r2_ud_u));
r2_matched_ud_d = cell(size(r2_ud_d));
if ~exist(paths.c2_matched_ud,'file')
    [c2_matched_ud_u,r2_matched_ud_u] = points_matching(c2_ud_u,r2_ud_u,c2_synth_ud_u,c2_matched_ud_u,r2_matched_ud_u,1);
    [c2_matched_ud_d,r2_matched_ud_d] = points_matching(c2_ud_d,r2_ud_d,c2_synth_ud_d,c2_matched_ud_d,r2_matched_ud_d,1);
    
    c2_matched_ud.u = c2_matched_ud_u;
    c2_matched_ud.d = c2_matched_ud_d;
    
    r2_matched_ud.u = r2_matched_ud_u;
    r2_matched_ud.d = r2_matched_ud_d;
else
    load(paths.c2_matched_ud)
    load(paths.r2_matched_ud)
    
    c2_matched_ud_u = c2_matched_ud.u;
    c2_matched_ud_d = c2_matched_ud.d;
    
    r2_matched_ud_u = r2_matched_ud.u;
    r2_matched_ud_d = r2_matched_ud.d;
end


%% Uncomment and run to check the dot detection

% figure('Position',[100 100 1000 1000])
% for ii = 1:num_target_plans
%    imshow(pattern1_ud{ii},[])
%    for jj = 1:length(c1_matched_ud_d{ii})
%         viscircles(c1_matched_ud_d{ii}(jj,:),r1_matched_ud_d{ii}(jj));
%         pause(0.01)
%    end
%    pause
% end
% close

%% Creation of the transform matrix

% We now create aa trasnform matrci to match the points of the two cameras.

c1_matched_ud_full{1} = [c1_matched_ud_u{1}; c1_matched_ud_d{1}];
c2_matched_ud_full{1} = [c2_matched_ud_u{1}; c2_matched_ud_d{1}];

tform_1to2 = fitgeotrans(c1_matched_ud_full{1},c2_matched_ud_full{1},'polynomial',4);
tform_2to1 = fitgeotrans(c2_matched_ud_full{1},c1_matched_ud_full{1},'polynomial',4);

if ~exist(paths.tform_1to2,'file')
    save(paths.tform_1to2,'tform_1to2')
    save(paths.tform_2to1,'tform_2to1')
end


%% Computation of the px / mm ratio

% We compute the mean distance between two grid points in pixels. We
% separate points belonging to the corners, borders and centers, because
% they do not have the same number of neighbours.
% We get the px/mm ratio by averaging the distance between two dot.

num_pts = length(c1_matched_ud_d{1});

ind_corners = false(1,num_pts);
ind_corners([1 rows_d num_pts-rows_d+1 num_pts]) = true;

ind_borders = false(1,num_pts);
ind_borders([1:rows_d rows_d+1:rows_d:num_pts-rows_d rows_d:rows_d:num_pts-rows_d+1 num_pts-rows_d+1:num_pts]) = true;
ind_borders(ind_corners) = false;

ind_centers = ~ind_borders & ~ind_corners;


coord_cam1 = c1_matched_ud_d{1};
coord_cam2 = c2_matched_ud_d{1};


pwd_cam1 = pdist2([coord_cam1(:,1) coord_cam1(:,2)],...
    [coord_cam1(:,1) coord_cam1(:,2)],'euclidean');

dist_corners_cam1 = sort(pwd_cam1(:,ind_corners));
dist_corners_cam1 = dist_corners_cam1(1:3,:);
dist_corners_cam1 = dist_corners_cam1(dist_corners_cam1 ~= 0);
dist_corners_cam1 = unique(dist_corners_cam1);

dist_borders_cam1 = sort(pwd_cam1(:,ind_borders));
dist_borders_cam1 = dist_borders_cam1(1:4,:);
dist_borders_cam1 = dist_borders_cam1(dist_borders_cam1 ~= 0);
dist_borders_cam1 = unique(dist_borders_cam1);

dist_centers_cam1 = sort(pwd_cam1(:,ind_centers));
dist_centers_cam1 = dist_centers_cam1(1:5,:);
dist_centers_cam1 = dist_centers_cam1(dist_centers_cam1 ~= 0);
dist_centers_cam1 = unique(dist_centers_cam1);

dist_cam1 = vertcat(dist_corners_cam1,dist_borders_cam1,dist_centers_cam1);


pwd_cam2 = pdist2([coord_cam2(:,1) coord_cam2(:,2)],...
    [coord_cam2(:,1) coord_cam2(:,2)],'euclidean');

dist_corners_cam2 = sort(pwd_cam2(:,ind_corners));
dist_corners_cam2 = dist_corners_cam2(1:3,:);
dist_corners_cam2 = dist_corners_cam2(dist_corners_cam2 ~= 0);
dist_corners_cam2 = unique(dist_corners_cam2);

dist_borders_cam2 = sort(pwd_cam2(:,ind_borders));
dist_borders_cam2 = dist_borders_cam2(1:4,:);
dist_borders_cam2 = dist_borders_cam2(dist_borders_cam2 ~= 0);
dist_borders_cam2 = unique(dist_borders_cam2);

dist_centers_cam2 = sort(pwd_cam2(:,ind_centers));
dist_centers_cam2 = dist_centers_cam2(1:5,:);
dist_centers_cam2 = dist_centers_cam2(dist_centers_cam2 ~= 0);
dist_centers_cam2 = unique(dist_centers_cam2);

dist_cam2 = vertcat(dist_corners_cam2,dist_borders_cam2,dist_centers_cam2);


if ~isfield(params_cam1,'gamma')
    params_cam1.gamma = mean(dist_cam1) / d_pts;
    save(paths.params_cam1,'params_cam1')
end

if ~isfield(params_cam2,'gamma')
    params_cam2.gamma = mean(dist_cam2) / d_pts;
    save(paths.params_cam2,'params_cam2')
end


%% Test stitching

% Here we test that the points matching can be done. the two images should
% be the same in the rgion of interest (i.e. inside the cell)

img1_ud = rect(double(img1),eye(3),...
    params_cam1.fc,params_cam1.cc,params_cam1.kc,params_cam1.KK_new);
img2_ud = rect(double(img2),eye(3),...
    params_cam2.fc,params_cam2.cc,params_cam2.kc,params_cam2.KK_new);


ref2Dinput = imref2d(size(img2_ud));
img1_t = imwarp(img1_ud,tform_1to2,'OutputView',ref2Dinput);

figure,imshow(img1_t,[])
figure,imshow(img2_ud,[])
figure,imshowpair(img1_t,img2_ud,'montage')
% figure,imshow(imsubtract(img0_t,img1_ud),[])


%%

ref2Dinput = imref2d(size(img1_ud));
img2_t = imwarp(img2_ud,tform_2to1,'OutputView',ref2Dinput);

figure,imshow(img2_t,[])
figure,imshow(img1_ud,[])
figure,imshow(imsubtract(img2_t,img1_ud),[])