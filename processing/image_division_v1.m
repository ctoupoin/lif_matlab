function im = image_division_v1(imgs_name1,imgs_name2,params_cam1,...
    params_cam2,mask_common_roi,tform_2to1,transfer_function_2over1,...
    bkgrd_no_laser1,bkgrd_no_laser2,filter_size,ii)

% Changelog:
%   - v1: Hello world !

% From the raw images, Regions Of Interest are created, which
% correspond to the part of the cell illuminated by the laser
[roi1,roi2] = geometric_transforms(imgs_name1,imgs_name2,params_cam1,...
    params_cam2,mask_common_roi,tform_2to1,ii);

roi1 = imsubtract(roi1,bkgrd_no_laser1);
roi2 = imsubtract(roi2,bkgrd_no_laser2);

% Intensity corrections are applied to the ROIs in roder to create a
% single image, im
roi1_c = roi1 .* transfer_function_2over1;

im = imdivide(roi1_c,roi2);
im = medfilt2(im,[filter_size,filter_size]);

