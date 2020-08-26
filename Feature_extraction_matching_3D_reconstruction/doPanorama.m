%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIDEO FOR MULTIPLE AND MOVING CAMERAS (VMMC)
% IPCV @ UAM
% Marcos Escudero-Vi�olo (marcos.escudero@uam.es)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%% PARAMETER DEFINITION %%
%%define dataset path
params.Directory    = fullfile('Datasets/malaga/');
params.resize       = 0.8;

%%Detector parameters
params.detector     =  'LoG_ss'; %'LoG_ss', 'SURF', 'SIFT', ...
params.nscales      =        5;
params.noctaves     =        3;
params.sigma0       =      1.6; % as we are using Matlab functions this is the minimum value allowed
params.npoints      =      300;

%%Descriptor parameters (see doc extractFeatures for explanation and additional parameters)
params.descriptor   =  'SIFT'; % 'SIFT', 'SURF', 'DSP-SIFT'...
params.desOnDecom   =    true; % describe on scale-space (linear or non-linear) decomposition (if available)
params.Upright      =   false; % set to true to avoid orientation estimation.
% for DSP-SIFT
params.dsp.ns       =      10;% number of sampled scales
params.dsp.sc_min   =     1/6;% smallest scale (relative to detection)
params.dsp.sc_max   =       3;% largest scale (relative to detection);    

%%Matching parameters (see doc matchFeatures for explanation and additional parameters)
params.MaxRatio     =   0.6;
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
ima{numImages}                     = [];
points{numImages}                  = [];
decomposition{numImages}           = [];
features{numImages}                = [];
tforms{numImages}                  = [];
matchedPoints{numImages,numImages} = [];
imageSize                          = rand(numImages,2);

%% get sigmas
k = 1;
params.sigmas = zeros(1,params.noctaves*params.nscales);
for o = 0:params.noctaves-1
    params.sigmas(k:(k+params.nscales-1)) = params.sigma0.*pow2([0:(params.nscales-1)]/params.nscales + o);
    k = k+params.nscales;
end
%% detect, describe, matching & transform
for j = 1:numImages
%% Load and convert images %%
ima{j}         =    readimage(params.Scene, j);
ima{j}         =     imresize(ima{j},params.resize);
gima           =    im2double(rgb2gray(ima{j}));
imageSize(j,:) =    size(gima);

%% initialize tform
tforms{j}      = projective2d(eye(3));

%% PoI Detection %%
sprintf('Detecting for image: %d',j)
[points{j},decomposition{j}] =  myDetector(gima,params);

%% PoI Description %%
sprintf('Describing for image: %d',j)
[features{j},points{j}]      =  myDescriptor(points{j},decomposition{j},params);

%% show detections
figure(j)
imshow(ima{j}); hold on;
plot(points{j});

% if j is not the first image, extract matchings and transformation with
% respect to the previous image
if j>1
indexPairs       = matchFeatures(features{j},features{j-1},'MaxRatio',params.MaxRatio,'Metric',params.Metric) ;
matchedPoints{j,j-1} =   points{j}(indexPairs(:,1));
matchedPoints{j-1,j} = points{j-1}(indexPairs(:,2));

% Estimate the transformation between I(j) and I(j-1).
tforms{j} = estimateGeometricTransform(matchedPoints{j,j-1} , matchedPoints{j-1,j},...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
    
% Compute T(n) * T(n-1) * ... * T(1)
tforms{j}.T = tforms{j}.T * tforms{j-1}.T; 
end

end


%% From Matlab's Panorama image stiching example:
% At this point, all the transformations in |tforms| are relative to the
% first image. This was a convenient way to code the image registration
% procedure because it allowed sequential processing of all the images.
% However, using the first image as the start of the panorama does not
% produce the most aesthetically pleasing panorama because it tends to
% distort most of the images that form the panorama. A nicer panorama can
% be created by modifying the transformations such that the center of the
% scene is the least distorted. This is accomplished by inverting the
% transform for the center image and applying that transform to all the
% others.
%
% Start by using the |projective2d| |outputLimits| method to find the
% output limits for each transform. The output limits are then used to
% automatically find the image that is roughly in the center of the scene.

%% Compute the output limits for each transform
xlim = imageSize;
ylim = imageSize;
for j = 1:numImages        
    [xlim(j,:), ylim(j,:)] = outputLimits(tforms{j}, [1 imageSize(j,2)], [1 imageSize(j,1)]);    
end

%% Extract center image
avgXLim = mean(xlim, 2);
avgYLim = mean(ylim, 2);

d = (avgXLim - mean(avgXLim)).^2 + (avgYLim - mean(avgYLim)).^2;
[~, idx] = sort(d);

centerIdx = floor((numel(tforms)+1)/2);
centerImageIdx = idx(centerIdx);

%% Apply the center image's inverse transform to all the others.
% i.e. align to middle image
Tinv = invert(tforms{centerImageIdx});
for j = 1:numel(tforms)    
    tforms{j}.T = tforms{j}.T * Tinv.T;
end

%% Create panorama
%% Initialize the Panorama
% Now, create an initial, empty, panorama into which all the images are
% mapped. 
% 
% Use the |outputLimits| method to compute the minimum and maximum output
% limits over all transformations. These values are used to automatically
% compute the size of the panorama.

for j = 1:numel(tforms)           
    [xlim(j,:), ylim(j,:)] = outputLimits(tforms{j}, [1 imageSize(j,2)], [1 imageSize(j,1)]);
end

maxImageSize = max(imageSize);

% Find the minimum and maximum output limits 
xMin = min([1; xlim(:)]);
xMax = max([maxImageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([maxImageSize(1); ylim(:)]);

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.
panorama = zeros([height width 3], 'like', ima{1});

%% Create the Panorama
% Use |imwarp| to map images into the panorama and use
% |vision.AlphaBlender| to overlay the images together.

blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');  

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

% Create the panorama.
for j = 1:numImages
   
    % Transform current image into the panorama.
    warpedImage = imwarp(ima{j}, tforms{j}, 'OutputView', panoramaView);
                  
    % Generate a binary mask.    
    mask = imwarp(true(imageSize(j,1),imageSize(j,2)), tforms{j}, 'OutputView', panoramaView);
    
    % Overlay the warpedImage onto the panorama.
    panorama = step(blender, panorama, warpedImage, mask);
end

figure
imshow(panorama)
title(sprintf('Panorama with %d views',numImages))






