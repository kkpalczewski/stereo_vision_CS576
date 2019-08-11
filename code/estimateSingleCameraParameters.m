function [cameraParams] = estimateSingleCameraParameters(imagePoints, boardSize, patchSize, imageSize)
% This function will estimate camera parameters (intrinsic, extrinsic) from
% checkerboard image points.

% Zhang's method consists of 5 parts
% 1. Estimate homography from checkerboard plane to screen space.
% 2. Calculate B matrix by solving Vb = 0.
% 3. Extract intrinsic parameters from B matrix.
% 4. Calculate extrinsic parameters from intrinsic parameters and homography.
% 5. Refine parameters using the maximum likelihood estimation.

% Function inputs:
% - 'imagePoints': positions of checkerboard points in a screen space.
% - 'boardSize': the number of horizontal, vertical patchs in the checkerboard.
% - 'patchSize': the size of the checkerboard patch in mm.
% - 'imageSize': the size of the checkerboard image in pixels.

% Function outputs:
% - 'cameraParams': a camera parameter includes intrinsic and extrinsic.

numView = size(imagePoints, 3);
numVerticalPatch = boardSize(1) - 1;
numHorizontalPatch = boardSize(2) - 1;
numCorner = size(imagePoints, 1);

%% Estimate a homography (appendix A)

% Generate checkerboard world points
worldPoints = zeros(size(imagePoints,1), size(imagePoints,2));
% Fill worldPoints (positions of checkerboard corners)
% ----- Your code here (1) ----- (slide 6)
for i = 1:length(worldPoints)
    worldPoints(i,:) = [floor((i - 1)/numVerticalPatch) * patchSize, mod((i-1), numVerticalPatch) * patchSize];
end

% Build L matrix
L = zeros(2 * numCorner, 9, numView);

% Fill L matrix
% ----- Your code here (2) ----- (slide 13)
for nv = 1:numView
    for cor = 1:numCorner
        % coordinates of points in model space (camera space) (u v)
        u = imagePoints(cor, 1, nv);
        v = imagePoints(cor, 2, nv);
        
        % NOTE: To normalization (2)
        % m_h = [u; v; 1];
        % m_h_norm = N*m_h;
        %Mm = num2cell(m_h_norm');
        %[u, v, ~] = deal(Mm{:});
        
        % coordinates of point in image space (world space)(X Y)
        X = worldPoints(cor, 1);
        Y = worldPoints(cor, 2);
        p = [X, Y, 1];
        % L matrix
        L(2*cor - 1, :, nv) = [p*-1, zeros(1,3), p*u];
        L(2*cor, :, nv) = [zeros(1, 3), p*-1, p*v];
    end
end

% Calculate a homography using SVD
% Fill homography matrix
% ----- Your code here (3) ----- (slide 15)
homography = zeros(3,3,numView);
for nv = 1:numView
    [~, Sh, Vh] = svd(L(:,:,nv));
    % search for index of min. singular value
    [~, index] = min(diag(Sh));
    Vht = Vh';
    homography(:, :, nv) = reshape(Vht(index, :), [3,3])';
end
% NOTE: To normalization (3)
for nv = 1:numView
   homography(:, :, nv) = homography(:, :, nv)/homography(3, 3, nv); 
end


%% Solve closed-form (section 3.1)
V = zeros(2 * numView, 6);
b = zeros(6, 1);

% Fill V matrix and calculate b vector
% ----- Your code here (4) ----- (slide 19, 23)

%Lambda functions for better read
% TODO: check if you are right
Vrow = @(k, l, H) [
    H(1,k)*H(1,l);
    H(1,k)*H(2,l) + H(2,k)*H(1,l);
    H(1,k)*H(3,l)+H(3,k)*H(1,l);
    H(2,k)*H(2,l);
    H(2,k)*H(3,l)+H(3,k)*H(2,l);
    H(3,k)*H(3,l)];
VfromH = @(H) [Vrow(1, 2, H)'; Vrow(1, 1, H)' - Vrow(2, 2, H)'];

for nv = 1:numView
    V(2*nv - 1:2*nv,:) = VfromH(homography(:,:,nv));
end
% compute b from SVD
[~, Sv, Vv] = svd(V);
% search for index of min. singular value
[~, index] = min(diag(Sv));
Vvt = Vv';
b = Vvt(index,:);

%% Extraction of the intrinsic parameters from matrix B (appendix B)
% ----- Your code here (5) ----- (slide 24)
B = num2cell(b);
% TODO: consult if this is a typo here
[B11, B12, B13, B22, B23, B33] = deal(B{:}); %this is in our homework
%[B11, B12, B22, B13, B23, B33] = deal(B{:});
v0 = (B12*B13-B11*B23)/(B11*B22-B12*B12);  % modify this line
lambda = B33 - (B13*B13+v0*(B12*B13-B11*B23))/B11; % modify this line
alpha = sqrt(lambda/B11);  % modify this line
beta = sqrt(lambda*B11/(B11*B22-B12*B12));  % modify this line
gamma = -B12*alpha*alpha*beta/lambda;  % modify this line
u0 = gamma*v0/beta-B13*alpha*alpha/lambda;  % modify this line

%% Estimate initial RT (section 3.1)
Rt = zeros(3, 4, numView);

% Fill Rt matrix
% ----- Your code here (6) ----- (slide 25, 26)
%camera intrinistic parameters
K = [[alpha, gamma, u0]; [0, beta, v0]; [0, 0, 1]];

KTest = K;

FunScaleFactor = @(K,H) (1/norm(K\H(:,1)) + 1/norm(K\H(:,2)))/2;
% TODO: if we check it
% make homography once again
%for nv = 1:numView
%   homogrpahy(:,:,nv) = inv(N)*homography(:,:,nv); 
%end
rvecsTest = zeros(numView, 3);
tvecsTest = zeros(numView, 3);
for nv = 1:numView
    lambdaPrim = FunScaleFactor(K, homography(:,:,nv)); 
    r1 = lambdaPrim*inv(K)*homography(:,1,nv);
    r2 = lambdaPrim*inv(K)*homography(:,2,nv);
    r3 = cross(r1, r2);
    t = lambdaPrim*inv(K)*homography(:,3,nv);
    [U, ~, V] = svd([r1, r2, r3]);
    %TODO: Explain why we can do this (slide 25)
    Rt(:, :, nv) = [U*V', t];
    rvecsTest(nv,:) = rotationMatrixToVector(Rt(:,1:3,nv))';
    tvecsTest(nv,:) = t';
end
%% Maximum likelihood estimation (section 3.2)
options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', ...
    'TolX', 1e-32, 'TolFun', 1e-32, 'MaxFunEvals', 1e64, ...
    'MaxIter', 1e64, 'UseParallel', true, 'Display', 'iter');

% Build initial x value as x0
% ----- Your code here (7) ----- (slide 29)

% 5 for intrinsic
% 3 for translation, 3 for rotation, total 6 for each checkerboard image
x0 = zeros(5 + 6 * size(imagePoints, 3), 1);  % modify this line
x0(1:5,1) = [alpha; beta; gamma; u0; v0];
for nv = 1:numView
    x0(6+(nv-1)*size(imagePoints, 3) : 6+nv*size(imagePoints, 3)-1, 1) = ...
    [rotationMatrixToVector(Rt(:,1:3,nv))'; Rt(:,4,nv)]; 
end

% Non-least square optimization
% Read [https://mathworks.com/help/optim/ug/lsqnonlin.html] for more information
[objective] = @(x) func_calibration(imagePoints, worldPoints, x);

[x_hat, ~, ~, ~, ~] = lsqnonlin(objective,x0,[],[],options);

Rttest = zeros(size(Rt));
for nv = 1:numView
   Rttest(:,:,nv) = [rotationVectorToMatrix(x_hat(6+(nv-1)*size(imagePoints, 3):...
       6+(nv-1)*size(imagePoints, 3)+2, 1)'),...
       x_hat(6+(nv-1)*size(imagePoints, 3)+3:6+(nv-1)*size(imagePoints, 3)+5, 1)];
end
%% Build camera parameters
rvecs = zeros(numView, 3);
tvecs = zeros(numView, 3);
K = [1, 0, 0
     0, 1, 0
     0, 0, 1];

% Extract intrinsic matrix K, rotation vectors and translation vectors from x_hat
% ----- Your code here (8) -----
K = [[x_hat(1,1), x_hat(3,1), x_hat(4,1)];...
    [0, x_hat(2,1), x_hat(5,1)];...
    [0, 0, 1]];

for nv=1:numView
    rvecs(nv, :) = x_hat(6+(nv-1)*6 : 6+(nv-1)*6+2, 1)';
    tvecs(nv, :) = x_hat(6+(nv-1)*6+3 : 6+(nv-1)*6+5, 1)';
    Rrr = rotationVectorToMatrix(x_hat(6+(nv-1)*6 : 6+(nv-1)*6+2, 1));
    rvecs(nv, :) = rotationMatrixToVector(Rrr');
end

% Generate cameraParameters structure
cameraParams = cameraParameters('IntrinsicMatrix', K', ...
    'RotationVectors', rvecs, 'TranslationVectors', tvecs, ...
    'WorldPoints', worldPoints, 'WorldUnits', 'mm', ...
    'imageSize', imageSize) ; 

cameraParamsTest = cameraParameters('IntrinsicMatrix', KTest', ...
    'RotationVectors', rvecsTest, 'TranslationVectors', tvecsTest, ...
    'WorldPoints', worldPoints, 'WorldUnits', 'mm', ...
    'imageSize', imageSize) ; 
% Uncomment this line after you implement this function to calculate
% reprojection errors of your camera parameters.
%cameraParams = cameraParamsTest

reprojected_errors = imagePoints - cameraParams.ReprojectedPoints;

cameraParams = cameraParameters('IntrinsicMatrix', K', ...
    'RotationVectors', rvecs, 'TranslationVectors', tvecs, ...
    'WorldPoints', worldPoints, 'WorldUnits', 'mm', ...
    'imageSize', imageSize, 'ReprojectionErrors', reprojected_errors) ; 

testImg = zeros(54,2);
R = rotationVectorToMatrix(x_hat(12:14));
t = x_hat(15:17);
for im = 1:54
   p = [worldPoints(im,1);worldPoints(im,2); 1];
   q = K*[R(:,1:2),t]*p;
   Q = num2cell(q);
   [u_hat, v_hat, w_hat] = deal(Q{:});
   testImg(im,:) = [u_hat/w_hat, v_hat/w_hat]; 
end
%cameraParams = cameraParamsTest
%mOrig = imread('..\code\checkerboard\left02.png');
%imshow(mOrig,'InitialMagnification',30);
%hold on
%plot(cameraParams.ReprojectedPoints(:,1,2),cameraParams.ReprojectedPoints(:,2,2),'g*-');
%hold on
%w =worldToImage(cameraParams,R',t,[worldPoints, ones(54,1)]);
%plot(w(:,1),w(:,2), 'm*-');
%hold on
%plot(cameraParamsTest.ReprojectedPoints(:,1,2),cameraParamsTest.ReprojectedPoints(:,2,2),'c*-');
%hold on
%plot(imagePoints(:,1,2),imagePoints(:,2,2),'r*-')
%hold on
%plot(testImg(:,1),testImg(:,2),'b*-')
%plot(cameraParams.ReprojectedPoints(:,1,3),cameraParams.ReprojectedPoints(:,2,3),'g*-');
%hold off

