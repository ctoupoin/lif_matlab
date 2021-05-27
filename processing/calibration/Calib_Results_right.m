% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly executed under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 3314.746471540921902 ; 3314.746471540921902 ];

%-- Principal point:
cc = [ 1003.896259605088517 ; 992.747621506495534 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.137653944106492 ; 0.104830913443659 ; 0.001896613138047 ; -0.007219352507646 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 75.123624699174798 ; 75.123624699174798 ];

%-- Principal point uncertainty:
cc_error = [ 4.789708602743419 ; 3.520442335786830 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.006724414600754 ; 0.026195283836280 ; 0.000133051782412 ; 0.000269612139678 ; 0.000000000000000 ];

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
omc_1 = [ 2.730246e-02 ; 1.303390e-02 ; -2.154608e-02 ];
Tc_1  = [ -1.806832e+02 ; -1.186403e+02 ; 8.824568e+02 ];
omc_error_1 = [ 1.467355e-03 ; 1.110380e-03 ; 7.093082e-05 ];
Tc_error_1  = [ 1.276784e+00 ; 9.389609e-01 ; 2.002383e+01 ];

%-- Image #2:
omc_2 = [ 3.069760e-02 ; -3.347207e-02 ; -2.067452e-02 ];
Tc_2  = [ -1.829529e+02 ; -1.189775e+02 ; 8.661986e+02 ];
omc_error_2 = [ 1.668661e-03 ; 1.207579e-03 ; 8.282997e-05 ];
Tc_error_2  = [ 1.260472e+00 ; 9.224615e-01 ; 1.963535e+01 ];

%-- Image #3:
omc_3 = [ 2.609561e-02 ; 7.676887e-02 ; -2.212475e-02 ];
Tc_3  = [ -1.738837e+02 ; -1.190578e+02 ; 8.856609e+02 ];
omc_error_3 = [ 1.134977e-03 ; 1.995072e-03 ; 9.070514e-05 ];
Tc_error_3  = [ 1.272734e+00 ; 9.391961e-01 ; 2.000362e+01 ];

%-- Image #4:
omc_4 = [ 2.414289e-02 ; 1.834965e-02 ; -2.151565e-02 ];
Tc_4  = [ -1.615300e+02 ; -1.188696e+02 ; 8.885829e+02 ];
omc_error_4 = [ 1.504319e-03 ; 1.171500e-03 ; 7.187398e-05 ];
Tc_error_4  = [ 1.285080e+00 ; 9.449857e-01 ; 2.015676e+01 ];

%-- Image #5:
omc_5 = [ 2.652521e-02 ; 1.954677e-02 ; -2.091167e-02 ];
Tc_5  = [ -1.752150e+02 ; -1.189839e+02 ; 8.605690e+02 ];
omc_error_5 = [ 1.385683e-03 ; 1.160672e-03 ; 7.003986e-05 ];
Tc_error_5  = [ 1.244322e+00 ; 9.154829e-01 ; 1.952445e+01 ];

%-- Image #6:
omc_6 = [ 2.915511e-02 ; -2.000001e-02 ; -2.077051e-02 ];
Tc_6  = [ -2.072312e+02 ; -1.189670e+02 ; 8.801775e+02 ];
omc_error_6 = [ 1.567501e-03 ; 1.074670e-03 ; 7.902249e-05 ];
Tc_error_6  = [ 1.278886e+00 ; 9.375074e-01 ; 1.996899e+01 ];

%-- Image #7:
omc_7 = [ 3.462681e-02 ; -4.466811e-02 ; -2.008421e-02 ];
Tc_7  = [ -2.015197e+02 ; -1.190350e+02 ; 8.656553e+02 ];
omc_error_7 = [ 1.715940e-03 ; 1.374356e-03 ; 9.224916e-05 ];
Tc_error_7  = [ 1.262171e+00 ; 9.221677e-01 ; 1.960980e+01 ];

%-- Image #8:
omc_8 = [ 3.466121e-02 ; -4.473740e-02 ; -2.008907e-02 ];
Tc_8  = [ -2.015190e+02 ; -1.190344e+02 ; 8.656472e+02 ];
omc_error_8 = [ 1.716589e-03 ; 1.375495e-03 ; 9.230191e-05 ];
Tc_error_8  = [ 1.262172e+00 ; 9.221592e-01 ; 1.960950e+01 ];

%-- Image #9:
omc_9 = [ 2.745620e-02 ; -1.974545e-02 ; -2.089923e-02 ];
Tc_9  = [ -1.578433e+02 ; -1.188755e+02 ; 8.767832e+02 ];
omc_error_9 = [ 1.673462e-03 ; 1.083169e-03 ; 7.560392e-05 ];
Tc_error_9  = [ 1.273135e+00 ; 9.332631e-01 ; 1.988357e+01 ];

%-- Image #10:
omc_10 = [ 2.631847e-02 ; 1.439670e-02 ; -2.143738e-02 ];
Tc_10  = [ -1.842971e+02 ; -1.188950e+02 ; 8.894706e+02 ];
omc_error_10 = [ 1.445720e-03 ; 1.127959e-03 ; 7.184465e-05 ];
Tc_error_10  = [ 1.286604e+00 ; 9.463587e-01 ; 2.018299e+01 ];

%-- Image #11:
omc_11 = [ -6.748410e-03 ; 1.222392e-02 ; -2.093990e-02 ];
Tc_11  = [ -1.761863e+02 ; 9.258781e+01 ; 8.849139e+02 ];
omc_error_11 = [ 1.048718e-03 ; 1.109314e-03 ; 1.102098e-04 ];
Tc_error_11  = [ 1.278837e+00 ; 9.394082e-01 ; 2.007411e+01 ];

%-- Image #12:
omc_12 = [ -1.151936e-02 ; -3.273752e-02 ; -2.025546e-02 ];
Tc_12  = [ -1.787612e+02 ; 9.241929e+01 ; 8.694410e+02 ];
omc_error_12 = [ 1.141263e-03 ; 1.227799e-03 ; 1.132746e-04 ];
Tc_error_12  = [ 1.263373e+00 ; 9.256836e-01 ; 1.970408e+01 ];

%-- Image #13:
omc_13 = [ -1.129069e-02 ; 7.519158e-02 ; -2.184734e-02 ];
Tc_13  = [ -1.691730e+02 ; 9.259023e+01 ; 8.881649e+02 ];
omc_error_13 = [ 1.005282e-03 ; 2.033010e-03 ; 1.184722e-04 ];
Tc_error_13  = [ 1.275197e+00 ; 9.461293e-01 ; 2.000201e+01 ];

%-- Image #14:
omc_14 = [ -6.193196e-03 ; 1.766591e-02 ; -2.089248e-02 ];
Tc_14  = [ -1.570559e+02 ; 9.232094e+01 ; 8.906491e+02 ];
omc_error_14 = [ 1.071612e-03 ; 1.157879e-03 ; 1.105161e-04 ];
Tc_error_14  = [ 1.286595e+00 ; 9.454316e-01 ; 2.019629e+01 ];

%-- Image #15:
omc_15 = [ -3.415018e-03 ; 1.910126e-02 ; -2.032045e-02 ];
Tc_15  = [ -1.708445e+02 ; 9.223464e+01 ; 8.631283e+02 ];
omc_error_15 = [ 9.889132e-04 ; 1.161793e-03 ; 1.101245e-04 ];
Tc_error_15  = [ 1.246520e+00 ; 9.164573e-01 ; 1.957311e+01 ];

%-- Image #16:
omc_16 = [ -1.096136e-02 ; -1.923443e-02 ; -2.035807e-02 ];
Tc_16  = [ -2.029397e+02 ; 9.237984e+01 ; 8.831488e+02 ];
omc_error_16 = [ 1.120307e-03 ; 1.110756e-03 ; 1.133375e-04 ];
Tc_error_16  = [ 1.281353e+00 ; 9.386978e-01 ; 2.003824e+01 ];

%-- Image #17:
omc_17 = [ -1.437012e-02 ; -4.420394e-02 ; -1.997096e-02 ];
Tc_17  = [ -1.974563e+02 ; 9.248699e+01 ; 8.694400e+02 ];
omc_error_17 = [ 1.165380e-03 ; 1.420815e-03 ; 1.171557e-04 ];
Tc_error_17  = [ 1.265792e+00 ; 9.267800e-01 ; 1.968727e+01 ];

%-- Image #18:
omc_18 = [ -1.447568e-02 ; -4.421644e-02 ; -1.998361e-02 ];
Tc_18  = [ -1.974574e+02 ; 9.248864e+01 ; 8.694429e+02 ];
omc_error_18 = [ 1.166461e-03 ; 1.421209e-03 ; 1.171519e-04 ];
Tc_error_18  = [ 1.265797e+00 ; 9.267807e-01 ; 1.968730e+01 ];

%-- Image #19:
omc_19 = [ -9.119573e-03 ; -2.015723e-02 ; -2.036058e-02 ];
Tc_19  = [ -1.535846e+02 ; 9.241230e+01 ; 8.793747e+02 ];
omc_error_19 = [ 1.140745e-03 ; 1.088766e-03 ; 1.116003e-04 ];
Tc_error_19  = [ 1.275309e+00 ; 9.352922e-01 ; 1.993970e+01 ];

%-- Image #20:
omc_20 = [ -8.206594e-03 ; 1.380453e-02 ; -2.075179e-02 ];
Tc_20  = [ -1.798230e+02 ; 9.234078e+01 ; 8.917829e+02 ];
omc_error_20 = [ 1.060594e-03 ; 1.129201e-03 ; 1.100600e-04 ];
Tc_error_20  = [ 1.288423e+00 ; 9.465844e-01 ; 2.022918e+01 ];

