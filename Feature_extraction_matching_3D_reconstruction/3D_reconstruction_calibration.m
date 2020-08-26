%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIDEO FOR MULTIPLE AND MOVING CAMERAS (VMMC)
% IPCV @ UAM
%JINGWEN YANG FINAL PROJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

% include ACT_lite path
ACT_path = './ACT_lite/';
addpath(genpath(ACT_path));
% include extra funs
extra_funs_path = './extra_funs/';
addpath(genpath(extra_funs_path));

%load('example_real_scene.mat')

%% PARAMETER DEFINITION %%
%%define dataset path
%params.Directory    = fullfile('data_final_project/bed');
params.Directory    = fullfile('data_final_project/bed');

%%Detector parameters
params.detector     =  'SURF'; %'LoG_ss', 'SURF', 'SIFT','DoG_ss', 'KAZE' 'DoH' ...
params.nscales      =       5;
params.noctaves     =        3;
params.sigma0       =      1.6; % as we are using Matlab functions this is the minimum value allowed
params.npoints      =      500;

%%Descriptor parameters (see doc extractFeatures for explanation and additional parameters)
params.descriptor   =  'SURF'; % 'SIFT', 'SURF', 'DSP-SIFT', 'KAZE'...
%params.descriptor   =  'SURF'; % 'SIFT', 'SURF', 'DSP-SIFT', 'KAZE'...
params.desOnDecom   =   false; % describe on scale-space (linear or non-linear) decomposition (if available)
params.Upright      =   false; % set to true to avoid orientation estimation.
% for DSP-SIFT
params.dsp.ns       =      10;% number of sampled scales
params.dsp.sc_min   =     1/6;% smallest scale (relative to detection)
params.dsp.sc_max   =       3;% largest scale (relative to detection);    
%%Matching parameters (see doc matchFeatures for explanation and additional parameters)
params.MaxRatio     =   0.45;  %change here to get more confident points 
params.Metric       =  'SSD';
%% END OF PARAMETER DEFINITION %%

%% addpaths
addpath(genpath('./detectors/'));
addpath(genpath('./descriptors/'));
addpath(genpath('./toolbox/'));

%% preload dataset
params.Scene = imageDatastore(params.Directory);
numImages    = numel(params.Scene.Files);

%% initialize (sort of)
ima{numImages}           = [];
points{numImages}        = [];
decomposition{numImages} = [];
features{numImages}      = [];

%% get sigmas
k = 1;
params.sigmas = zeros(1,params.noctaves*params.nscales);
for o = 0:params.noctaves-1
    params.sigmas(k:(k+params.nscales-1)) = params.sigma0.*pow2([0:(params.nscales-1)]/params.nscales + o);
    k = k+params.nscales;
end

%% detect & describe
for j = 1:numImages
%% Load and convert images %%
ima{j}      =      readimage(params.Scene, j);
%ima{j}     =     imresize(ima{j},0.25); %where we resize the image to reduce computation cost
gima        =      im2double(rgb2gray(ima{j}));

%% PoI Detection %%
sprintf('Detecting for image: %d',j)
[points{j},decomposition{j}] =  myDetector(gima,params);

%% PoI Description %%
sprintf('Describing for image: %d',j)
[features{j},points{j}]      =  myDescriptor(points{j},decomposition{j},params);

%% show detections
% figure(j)
% imshow(ima{j}); hold on;
% plot(points{j},'showOrientation',true);

end

% ------------------------------------------------------------------------
% 1. compute consistent point matches from N views
% ------------------------------------------------------------------------
point_matrix = n_view_matching(points, features, ima, params.MaxRatio, params.Metric);
npoints = size(point_matrix,2);
point_matrix_homo = ones(3,npoints,numImages);
point_matrix_homo(1:2,:,:)= point_matrix;

% ------------------------------------------------------------------------
% 2. Compute the fundamental matrix using the first and last cameras
% of the camera set (N cameras)
% ------------------------------------------------------------------------
q_2cams(:,:,1)=point_matrix_homo(:,:,1); 
q_2cams(:,:,2)=point_matrix_homo(:,:,end);

[F, P_2cam_est,Q_2cam_est,q_2cam_est] = MatFunProjectiveCalib(q_2cams);
%vgg_gui_F(ima{1},ima{numImages},F');
disp(['Residual reprojection error. 8 point algorithm  = ' num2str( ErrorRetroproy(q_2cams,P_2cam_est,Q_2cam_est)/2 )]);
draw_reproj_error(q_2cams,P_2cam_est,Q_2cam_est);

% ------------------------------------------------------------------------
% 3.a.i Resection. Obtain the projection matrices of the rest of cameras using the PDLT_NA function 
% ------------------------------------------------------------------------
P_Ncam_est = zeros(3,4,numImages);
P_Ncam_est(:,:,1) = P_2cam_est(:,:,1);
P_Ncam_est(:,:,numImages) = P_2cam_est(:,:,2);
for i=1:numImages
    P_Ncam_est(:,:,i) = PDLT_NA(point_matrix_homo(:,:,i),Q_2cam_est,0,0);
end
% Compute the statistics of the reprojection error for the initial projective reconstruction
disp(['Residual reprojection err resectioning  = ' num2str( ErrorRetroproy(point_matrix_homo,P_Ncam_est,Q_2cam_est)/2 )]);
draw_reproj_error(point_matrix_homo,P_Ncam_est,Q_2cam_est);
% ------------------------------------------------------------------------
% 3.a.ii Projective Bundle Adjustment. Use BAProjectiveCalib function
% Coordinates of 3D and 2D points are given in homogeneus coordinates
% ------------------------------------------------------------------------
% auxiliary matrix that indicates that all points are visible in all the cameras
vp = ones(npoints,numImages);
[P_BA,Q_BA,x_reprojected_BA] = BAProjectiveCalib(point_matrix_homo,P_Ncam_est,Q_2cam_est,vp);
%Compute the statistics of the reprojection error for the improved projective reconstruction
disp(['Residual reprojection err Bundle Adjustment = ' num2str( ErrorRetroproy(point_matrix_homo,P_BA,Q_BA)/2)]);
draw_reproj_error(point_matrix_homo,P_BA,Q_BA);
%---------------------------------------------------------------------
% 4. Recompute F matrix between two of the cameras, using the projection
% matrices obtained after the projective Bundle Adjustment. use
% vgg_F_from_P 
% ------------------------------------------------------------------------
F_recomputed = vgg_F_from_P(P_BA(:,:,1),P_BA(:,:,end));
%vgg_gui_F(ima{1},ima{numImages},F_recomputed');

%---------------------------------------------------------------------
% 5. use the properties of essential matrix(between two cameras) to obtain
% an euclidean reconstruction of teh scene (use reprojected points after
% Bundle Adjustment).
% ------------------------------------------------------------------------
x_reprojected_BA_2cam(:,:,1) = x_reprojected_BA(:,:,1);
x_reprojected_BA_2cam(:,:,2) = x_reprojected_BA(:,:,end);

%K from the small checkerboard
%K(:,:,1) = [1266.02028234893, -8.07083221315575, 0.810107434933; 0, 1267.55898122649, 0.787326298565; 0,0, 1];
%K(:,:,2) = [1266.02028234893, -8.07083221315575, 0.810107434933; 0, 1267.55898122649, 0.787326298565; 0,0, 1];
K_s = [1266.02028234893, -8.07083221315575, 812.810107434933; 0, 1267.55898122649, 588.787326298565; 0,0,1];
K_s(1,2) = acos(-K_s(1,2)/K_s(1,1));
K_s(2,2) = K_s(2,2)/sin(K_s(1,2));

K(:,:,1) = K_s;
K(:,:,2) = K_s;
 
%if we scale down the image, we should also scale down the intrinsic params
% K(:,:,1) = K(:,:,1)/4;
% K(:,:,2) = K(:,:,2)/4;
% 
% K(3,3,1)= 1;
% K(3,3,2)= 1;
  

[E,P_euc,Q_euc] = EuclideanReconstructionUsingEssentialMatrix(q_2cams,x_reprojected_BA_2cam,F_recomputed,K);

disp('************************************* END')


