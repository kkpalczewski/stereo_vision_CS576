function [meanError] = estimateDepth2(leftImage, rightImage, ...
    stereoParameters, wVec, mindispVec, maxdispVec, trueDepth)
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

leftImageGray = rgb2gray(im2double(leftImage));
rightImageGray = rgb2gray(im2double(rightImage));
% For debugging!!!
% TODO: comment when mature algorithm
%leftImageGray = leftImageGray(380:420,:);
%rightImageGray = rightImageGray(380:420,:);

translation = stereoParameters.TranslationOfCamera2;
baseline = norm(translation); %in Pixels
focalLength = stereoParameters.CameraParameters1.FocalLength(1); %in Pixels

disparityMap = zeros(size(leftImageGray));
depthMap = zeros(size(leftImageGray));
% ----- Your code here (10) -----
meanError = zeros(numel(wVec), numel(maxdispVec), numel(mindispVec));

%padded fixed matrix

% Left image is strict in a place
% For each depth plane (so for each translation on x)
%  For each pixel in the comoposite image, compute variance
W2 = 12;
W2matrix = ones(W2, W2);

MinDisp_index = 1;
for MinDisparity = mindispVec
    MaxDisp_index = 1;
    for MaxDisparity = maxdispVec
        W_index = 1;
        for W1 = wVec
            NCCMatrix = zeros([size(leftImageGray),MaxDisparity-MinDisparity+1]);
            d_index = 1;
            for  depth = MinDisparity:MaxDisparity
                % padded sliding matrix
                rightImageSlide = [zeros(size(rightImageGray,1),depth), ...
                           rightImageGray(:,1:size(rightImageGray,2)-depth)];
                W1matrix = ones(W1,W1);
                A = imfilter(leftImageGray, W1matrix, 'replicate');
                AA = imfilter(leftImageGray.*leftImageGray, W1matrix, 'replicate');
                AE = A/(W1*W1);
                B = imfilter(rightImageSlide, W1matrix, 'replicate');
                BB = imfilter(rightImageSlide.*rightImageSlide, W1matrix, 'replicate');
                BE = B/(W1*W1);

                NOM = A.*B - B.*AE - A.*BE + AE.*BE;
                DOM1 = sqrt(AA);
                DOM2 = sqrt(BB);
                DOM3 = DOM2.*DOM1;
                NEW_NCC = NOM./DOM3;
                NEW_NCC = imfilter(NEW_NCC, W2matrix, 'replicate');
                %NCCMatrix(:,:,d_index = NEW_NCC
                NEW_NCC = NEW_NCC(:,d_index:size(NEW_NCC, 2));
                NCCMatrix(:,d_index:size(NCCMatrix, 2),d_index) = NEW_NCC;
                d_index = d_index + 1;
            end

            % take maximum
            disparityMapTmp = ones(size(NCCMatrix,1), size(NCCMatrix,2))*MinDisparity;
            [~,disparityMap] = max(NCCMatrix,[],3);
            disparityMap = disparityMap + disparityMapTmp-1;
            depthMap = (focalLength*baseline)./disparityMap;
            depthMap = fillmissing(depthMap, 'movmean', 5);
            %depthMap = wiener2(depthMap,[5 5]);
            depthDiff = trueDepth(:, 140:end-140) - depthMap(:, 140:end-140);
            meanError(W_index, MaxDisp_index, MinDisp_index) = mean(abs(depthDiff(:)));
                
            fprintf('W: [%.f], MinDisparity: [%.f], MaxDisparity: [%.f]\n', W1, MinDisparity, MaxDisparity);
            W_index = W_index + 1;
        end
        MaxDisp_index = MaxDisp_index + 1;
    end
    MinDisp_index = MinDisp_index + 1;
end