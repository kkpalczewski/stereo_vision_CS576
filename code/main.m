% For this project, you will need to report performance of camera
% calibration and depth estimation.
% The starter code is initialized to placeholder just so that the starter
% code does not crash when run unmodified and you can get a preview of how
% results are presented.

% You should not modify 'main.m' file as we will grade your homework by
% replacing functions by your implementations. 

%% Generate checkerboard points
numView = 6;  % Number of checkerboard images

leftImageFiles = cell(1, numView);
rightImageFiles = cell(1, numView);
for i=1:numView
    leftImageFiles{1,i} = sprintf('checkerboard/left%02d.png', i);
    rightImageFiles{1,i} = sprintf('checkerboard/right%02d.png', i);
end

[imagePoints, boardSize] = detectCheckerboardPoints(leftImageFiles, rightImageFiles);

%% Estimate stereo camera parameters
patchSize = 30;  % patch size of checkerboard in mm.
imageSize = size(imread(leftImageFiles{1}));  % image size in pixels.
imageSize = imageSize(1:2);

cameraParameters1 = estimateSingleCameraParameters(imagePoints(:, :, :, 1), boardSize, patchSize, imageSize);
cameraParameters2 = estimateSingleCameraParameters(imagePoints(:, :, :, 2), boardSize, patchSize, imageSize);

% Calculate relative R, t between two cameras.
[R, t] = vision.internal.calibration.estimateInitialTranslationAndRotation(cameraParameters1, cameraParameters2);

% Generate stereo parameters structure
stereoParams = stereoParameters(cameraParameters1, cameraParameters2, R, t);

%% Evaluate rectification
figure; showReprojectionErrors(stereoParams);

for i=1:2
    fprintf('Evaluation [%02d]\n', i);
    leftEval = imread(sprintf('checkerboard/left_eval%02d.png', i));
    rightEval = imread(sprintf('checkerboard/right_eval%02d.png', i));

    [unrectifiedImagePoints, ~] = detectCheckerboardPoints(leftEval, rightEval);
    unrectifiedDiff = abs(unrectifiedImagePoints(:, 2, :, 1) - unrectifiedImagePoints(:, 2, :, 2));
    fprintf('Original y coordinate mean difference: %.4f\n', mean(unrectifiedDiff(:)));

    [leftRrectifiedEval, rightRectifiedEval] = rectifyStereoImages(leftEval, rightEval, stereoParams);
    [rectifiedImagePoints, ~] = detectCheckerboardPoints(leftRrectifiedEval, rightRectifiedEval);
    rectifiedDiff = abs(rectifiedImagePoints(:, 2, :, 1) - rectifiedImagePoints(:, 2, :, 2));
    fprintf('Rectified y coordinate mean difference: %.4f\n', mean(rectifiedDiff(:)));
end
fprintf('\n');

%% Estimate depth of stereo images
%scenes = {'scene1', 'scene2'};
scenes = {'scene2'}
for scene=scenes
    scene = scene{1};
    fprintf('Start depth estimation [%s]\n', scene);
    %% Retify scene images
    leftImage = imread(sprintf('%s/left.png', scene));
    rightImage = imread(sprintf('%s/right.png', scene));
    [leftRectified, rightRectified] = rectifyStereoImages(leftImage, rightImage, stereoParams);
    figure; imshow([leftRectified, rightRectified]);
    
    imwrite(leftRectified, sprintf('output/left_%s_rectified.png', scene));
    imwrite(rightRectified, sprintf('output/right_%s_rectified.png', scene));
    
    %% Prepare rectified depth
    % NOTE: what it is about??
    load(sprintf('%s/gt_depthmap', scene));
    [leftGtDepthRectified, rightGtDepthRectified] = rectifyStereoImages(leftGtDepth, rightGtDepth, stereoParams);
    
    %% Stereo matching
    [depthMap, disparityMap] = estimateDepth(leftRectified, rightRectified, stereoParams);
    figure; imagesc(disparityMap); colorbar;
    figure; imagesc(depthMap); colorbar;
    % NOTE: cross validation for finding parameters
    %W = [3,5];
    %minDisparity = [10];
    %maxDisparity = [100];
    %CrossVal = estimateDepth2(leftRectified, rightRectified, stereoParams,...
    %    W, minDisparity, maxDisparity, leftGtDepthRectified);
    
    % TODO: check if images which are closer has higher or lower disparity
    depthDiff = leftGtDepthRectified(:, 140:end-140) - depthMap(:, 140:end-140);
    fprintf('Depth mean difference: %.2f\n', mean(abs(depthDiff(:))));
    fprintf('\n');
end