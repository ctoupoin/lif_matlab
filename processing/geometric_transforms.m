function [roi1,roi2] = geometric_transforms(imgs_name1,imgs_name2,params_cam1,...
    params_cam2,mask_common_roi,tform,ii)

img1 = imrotate(double(imread(imgs_name1(ii,:))),-90);
img2 = imrotate(double(imread(imgs_name2(ii,:))),90);

img1_ud = rect(img1,eye(3),...
    params_cam1.fc,params_cam1.cc,params_cam1.kc,params_cam1.KK_new);
img2_ud = rect(img2,eye(3),...
    params_cam2.fc,params_cam2.cc,params_cam2.kc,params_cam2.KK_new);

ref2Dinput = imref2d(size(img1_ud));
img2_ud_t = imwarp(img2_ud,tform,'OutputView',ref2Dinput);

% ref2Dinput = imref2d(size(img2_ud));
% img1_ud_t = imwarp(img1_ud,tform,'OutputView',ref2Dinput);

% roi1 = img1_ud_t(mask_common_roi(1):mask_common_roi(2),mask_common_roi(3):mask_common_roi(4));
% roi2 = img2_ud(mask_common_roi(1):mask_common_roi(2),mask_common_roi(3):mask_common_roi(4));

roi1 = img1_ud(mask_common_roi(1):mask_common_roi(2),mask_common_roi(3):mask_common_roi(4));
roi2 = img2_ud_t(mask_common_roi(1):mask_common_roi(2),mask_common_roi(3):mask_common_roi(4));

end