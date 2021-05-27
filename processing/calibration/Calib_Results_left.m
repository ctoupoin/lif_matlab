% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly executed under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 3029.560690116482419 ; 3029.560690116482419 ];

%-- Principal point:
cc = [ 1067.030590537458238 ; 993.181249343633567 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.106832062955886 ; 0.057756063006174 ; 0.002480803617577 ; -0.005331947380848 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 93.676718408053333 ; 93.676718408053333 ];

%-- Principal point uncertainty:
cc_error = [ 4.018409990244852 ; 3.378443251690364 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.006949452779252 ; 0.023939300725728 ; 0.000138629213306 ; 0.000281871094872 ; 0.000000000000000 ];

%-- Image size:
nx = 2048;
ny = 2048;


%-- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):
%-- Those variables are used to control which intrinsic parameters should be optimized

n_ima = 20;						% Number of calibration images
est_fc = [ 1 ; 1 ];					% Estimation indicator of the two focal variables
est_aspect_ratio = 0;				% Estimation indicator of the aspect ratio fc(2)/fc(1)
center_optim = 1;					% Estimation indicator of the principal point
est_alpha = 0;						% Estimation indicator of the skew coefficient
est_dist = [ 1 ; 1 ; 1 ; 1 ; 0 ];	% Estimation indicator of the distortion coefficients


%-- Extrinsic parameters:
%-- The rotation (omc_kk) and the translation (Tc_kk) vectors for every calibration image and their uncertainties

%-- Image #1:
omc_1 = [ 1.926157e-02 ; 5.304016e-03 ; -2.453958e-02 ];
Tc_1  = [ -1.162747e+02 ; -1.178097e+02 ; 8.071497e+02 ];
omc_error_1 = [ 1.534549e-03 ; 9.113257e-04 ; 7.497985e-05 ];
Tc_error_1  = [ 1.072125e+00 ; 9.013735e-01 ; 2.496371e+01 ];

%-- Image #2:
omc_2 = [ 2.528314e-02 ; -3.916320e-02 ; -2.388628e-02 ];
Tc_2  = [ -1.181912e+02 ; -1.181399e+02 ; 7.920608e+02 ];
omc_error_2 = [ 1.468330e-03 ; 1.513162e-03 ; 8.266645e-05 ];
Tc_error_2  = [ 1.055549e+00 ; 8.843599e-01 ; 2.448008e+01 ];

%-- Image #3:
omc_3 = [ 2.532114e-02 ; 6.193670e-02 ; -2.522958e-02 ];
Tc_3  = [ -1.093887e+02 ; -1.182181e+02 ; 8.098644e+02 ];
omc_error_3 = [ 1.423237e-03 ; 2.114550e-03 ; 9.315965e-05 ];
Tc_error_3  = [ 1.065480e+00 ; 9.007065e-01 ; 2.497754e+01 ];

%-- Image #4:
omc_4 = [ 1.745932e-02 ; 1.019382e-02 ; -2.445441e-02 ];
Tc_4  = [ -9.720381e+01 ; -1.180874e+02 ; 8.129172e+02 ];
omc_error_4 = [ 1.488800e-03 ; 9.628732e-04 ; 7.702123e-05 ];
Tc_error_4  = [ 1.079237e+00 ; 9.075982e-01 ; 2.513890e+01 ];

%-- Image #5:
omc_5 = [ 1.994269e-02 ; 1.018818e-02 ; -2.388286e-02 ];
Tc_5  = [ -1.103343e+02 ; -1.181788e+02 ; 7.869488e+02 ];
omc_error_5 = [ 1.478357e-03 ; 9.320465e-04 ; 7.387211e-05 ];
Tc_error_5  = [ 1.044813e+00 ; 8.787496e-01 ; 2.433688e+01 ];

%-- Image #6:
omc_6 = [ 2.352555e-02 ; -2.499653e-02 ; -2.390148e-02 ];
Tc_6  = [ -1.427676e+02 ; -1.180440e+02 ; 8.048823e+02 ];
omc_error_6 = [ 1.604443e-03 ; 1.188319e-03 ; 7.970748e-05 ];
Tc_error_6  = [ 1.072625e+00 ; 8.990897e-01 ; 2.488843e+01 ];

%-- Image #7:
omc_7 = [ 2.882892e-02 ; -4.966381e-02 ; -2.338259e-02 ];
Tc_7  = [ -1.367299e+02 ; -1.181465e+02 ; 7.914969e+02 ];
omc_error_7 = [ 1.521179e-03 ; 1.780949e-03 ; 8.795469e-05 ];
Tc_error_7  = [ 1.055631e+00 ; 8.833191e-01 ; 2.445066e+01 ];

%-- Image #8:
omc_8 = [ 2.884551e-02 ; -4.965834e-02 ; -2.338617e-02 ];
Tc_8  = [ -1.367304e+02 ; -1.181463e+02 ; 7.914974e+02 ];
omc_error_8 = [ 1.521469e-03 ; 1.780843e-03 ; 8.795232e-05 ];
Tc_error_8  = [ 1.055632e+00 ; 8.833196e-01 ; 2.445067e+01 ];

%-- Image #9:
omc_9 = [ 2.171419e-02 ; -2.648764e-02 ; -2.401028e-02 ];
Tc_9  = [ -9.327147e+01 ; -1.181159e+02 ; 8.019556e+02 ];
omc_error_9 = [ 1.414302e-03 ; 1.239724e-03 ; 7.979654e-05 ];
Tc_error_9  = [ 1.067638e+00 ; 8.956459e-01 ; 2.479477e+01 ];

%-- Image #10:
omc_10 = [ 1.931304e-02 ; 6.712335e-03 ; -2.442965e-02 ];
Tc_10  = [ -1.199922e+02 ; -1.180363e+02 ; 8.135295e+02 ];
omc_error_10 = [ 1.554843e-03 ; 9.260205e-04 ; 7.546996e-05 ];
Tc_error_10  = [ 1.080380e+00 ; 9.084605e-01 ; 2.516105e+01 ];

%-- Image #11:
omc_11 = [ 3.825060e-03 ; 4.530457e-03 ; -2.403071e-02 ];
Tc_11  = [ -1.112459e+02 ; 9.327677e+01 ; 8.092903e+02 ];
omc_error_11 = [ 9.846925e-04 ; 8.758020e-04 ; 1.140014e-04 ];
Tc_error_11  = [ 1.073414e+00 ; 9.019835e-01 ; 2.502761e+01 ];

%-- Image #12:
omc_12 = [ -1.356788e-03 ; -3.779424e-02 ; -2.294084e-02 ];
Tc_12  = [ -1.134328e+02 ; 9.313727e+01 ; 7.952429e+02 ];
omc_error_12 = [ 1.000395e-03 ; 1.494306e-03 ; 1.152579e-04 ];
Tc_error_12  = [ 1.058542e+00 ; 8.908425e-01 ; 2.455956e+01 ];

%-- Image #13:
omc_13 = [ -2.615609e-03 ; 5.959793e-02 ; -2.487204e-02 ];
Tc_13  = [ -1.041128e+02 ; 9.326132e+01 ; 8.123112e+02 ];
omc_error_13 = [ 1.080307e-03 ; 2.038876e-03 ; 1.208499e-04 ];
Tc_error_13  = [ 1.068273e+00 ; 9.097689e-01 ; 2.500132e+01 ];

%-- Image #14:
omc_14 = [ 4.345516e-03 ; 1.016251e-02 ; -2.416078e-02 ];
Tc_14  = [ -9.221646e+01 ; 9.299953e+01 ; 8.147568e+02 ];
omc_error_14 = [ 9.864567e-04 ; 9.294598e-04 ; 1.165067e-04 ];
Tc_error_14  = [ 1.080090e+00 ; 9.082495e-01 ; 2.519099e+01 ];

%-- Image #15:
omc_15 = [ 6.338778e-03 ; 1.007665e-02 ; -2.345851e-02 ];
Tc_15  = [ -1.054392e+02 ; 9.290517e+01 ; 7.893172e+02 ];
omc_error_15 = [ 9.561899e-04 ; 9.000178e-04 ; 1.135137e-04 ];
Tc_error_15  = [ 1.046327e+00 ; 8.798352e-01 ; 2.440505e+01 ];

%-- Image #16:
omc_16 = [ -7.737717e-04 ; -2.437506e-02 ; -2.317190e-02 ];
Tc_16  = [ -1.379050e+02 ; 9.318205e+01 ; 8.076807e+02 ];
omc_error_16 = [ 1.023287e-03 ; 1.173973e-03 ; 1.138244e-04 ];
Tc_error_16  = [ 1.075014e+00 ; 9.020598e-01 ; 2.496934e+01 ];

%-- Image #17:
omc_17 = [ -4.717975e-03 ; -4.825266e-02 ; -2.253739e-02 ];
Tc_17  = [ -1.320832e+02 ; 9.325536e+01 ; 7.952062e+02 ];
omc_error_17 = [ 1.037597e-03 ; 1.773090e-03 ; 1.167874e-04 ];
Tc_error_17  = [ 1.059757e+00 ; 8.927627e-01 ; 2.453418e+01 ];

%-- Image #18:
omc_18 = [ -4.685737e-03 ; -4.825069e-02 ; -2.253744e-02 ];
Tc_18  = [ -1.320832e+02 ; 9.325550e+01 ; 7.952055e+02 ];
omc_error_18 = [ 1.037547e-03 ; 1.772904e-03 ; 1.167890e-04 ];
Tc_error_18  = [ 1.059756e+00 ; 8.927600e-01 ; 2.453419e+01 ];

%-- Image #19:
omc_19 = [ 6.403328e-04 ; -2.541745e-02 ; -2.328808e-02 ];
Tc_19  = [ -8.847501e+01 ; 9.308522e+01 ; 8.046321e+02 ];
omc_error_19 = [ 9.755088e-04 ; 1.224609e-03 ; 1.163610e-04 ];
Tc_error_19  = [ 1.069557e+00 ; 8.997366e-01 ; 2.486734e+01 ];

%-- Image #20:
omc_20 = [ 2.774011e-03 ; 6.633745e-03 ; -2.394062e-02 ];
Tc_20  = [ -1.149830e+02 ; 9.307241e+01 ; 8.156019e+02 ];
omc_error_20 = [ 9.977270e-04 ; 8.906595e-04 ; 1.139872e-04 ];
Tc_error_20  = [ 1.081463e+00 ; 9.088648e-01 ; 2.522234e+01 ];

