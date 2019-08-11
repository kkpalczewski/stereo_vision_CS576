function [depthMap, disparityMap] = estimateDepth(leftImage, rightImage, stereoParameters)
% This function estimate disparity and depth values from left and right
% images. You should calculate disparty map first and then convert the
% disparity map to depth map using left camera parameters.

% Function inputs:
% - 'leftImage': rectified left image.
% - 'rightImage': rectified right image.
% - 'stereoParameters': stereo camera parameters.

% Function outputs:
% - 'depth': depth map of left camera.
% - 'disparity': disparity map of left camera.

MinDisparity = 6;
MaxDisparity = 68;

W1 = 5;
W2 = 12;
W2matrix = ones(W2, W2);

leftImageGray = rgb2gray(im2double(leftImage));
rightImageGray = rgb2gray(im2double(rightImage));

% TODO: comment when mature algorithm
%leftImageGray = leftImageGray(380:420,:);
%rightImageGray = rightImageGray(380:420,:);

NCCMatrix = zeros([size(leftImageGray),MaxDisparity-MinDisparity+1]);

% For debugging!!!


translation = stereoParameters.TranslationOfCamera2;
baseline = norm(translation); %in Pixels
focalLength = stereoParameters.CameraParameters1.FocalLength(1); %in Pixels

disparityMap = zeros(size(leftImageGray));
depthMap = zeros(size(leftImageGray));
% ----- Your code here (10) -----
% Left image is strict in a place
% For each depth plane (so for each translation on x)
%  For each pixel in the comoposite image, compute variance
d_index = 1;
w = (W1-1)/2;
leftImagePadded = padarray(leftImageGray(:,1:size(leftImageGray, 2)), [w w], 'replicate');

for  depth = MinDisparity:MaxDisparity
    fprintf('Start depth NCC estimation [%f.0]\n', depth);
    translatedRight = imtranslate(rightImageGray, [depth 0]);
    rightImagePadded = padarray(translatedRight, [w w], 'replicate');

    for x = 1+w : size(leftImagePadded,2)-w
        for y = 1+w : size(leftImagePadded, 1)-w
            leftPatch = leftImagePadded(y-w:y+w, x-w:x+w);
            rightPatch = rightImagePadded(y-w:y+w, x-w:x+w);
            NCCMatrix(y - w, x - w, d_index) = corr2(leftPatch, rightPatch); 
        end
    end
    NCCMatrix(:, :, d_index) = imfilter(NCCMatrix(:, :, d_index), W2matrix, 'replicate');
    d_index = d_index + 1;
end

% take maximum
disparityMapTmp = ones(size(disparityMap,1), size(disparityMap,2))*MinDisparity;

NCCMax = ones(size(NCCMatrix, 1), size(NCCMatrix, 2))*-1;

for x = MinDisparity+1:size(disparityMap,2)
   for y = 1:size(disparityMap,1)
       for z = 1:size(NCCMatrix,3)
           if NCCMatrix(y, x, z) >= NCCMax(y, x)
               NCCMax(y, x) = NCCMatrix(y, x, z);
               disparityMap(y, x) = z;
           end
       end
   end
end

disparityMap = disparityMap + disparityMapTmp;
%print(disparityMap)
%imshow(disparityMap)
depthMap = (focalLength*baseline)./disparityMap;


