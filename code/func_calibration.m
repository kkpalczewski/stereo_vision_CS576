function [objective] = func_calibration(imagePoints, worldPoints, x)
% Objective function to minimize eq.10 in Zhang's paper. 
% Size of input variable x is 5+6*n where n is number of checkerboard 
% images. An intrinsic matrix can be reconstructed from first five
% parameters, and the extrinsic matrix can be reconstructed from remain
% parameters.

% You should fill the variable hat_m which contains reprojected positions 
% of checkerboard points in screen coordinate.

% Function inputs:
% - 'imagePoints': positions of checkerboard points in a screen space.
% - 'worldPoints': positions of checkerboard points in a model space.
% - 'x': parameters to be optimized.

% Function outputs:
% - 'objective': difference of estimated values and real values.
    
numView = size(imagePoints,3);
hat_m = zeros(size(imagePoints));
% ----- Your code here (9) -----
% x0(1:5,1) = [alpha; beta; gamma; u0; v0]
% K = [[alpha, gamma, u0]; [0, beta, v0]; [0, 0, 1]];
K = [[x(1,1), x(3,1), x(4,1)];...
    [0, x(2,1), x(5,1)];...
    [0, 0, 1]];

for nv=1:numView
    R = rotationVectorToMatrix(x(6+(nv-1)*6 : 6+(nv-1)*6+2, 1)');
    t = x(6+(nv-1)*6+3 : 6+(nv-1)*6+5, 1);
    for im=1:size(imagePoints,1)
        p = [worldPoints(im,1);worldPoints(im,2); 1];
        q = K*[R(:,1:2),t]*p;
        Q = num2cell(q);
        [u_hat, v_hat, w_hat] = deal(Q{:});
        hat_m(im,:,nv) = [u_hat/w_hat, v_hat/w_hat];
    end
end
objective = imagePoints - hat_m;
