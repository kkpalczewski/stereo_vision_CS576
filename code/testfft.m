clear all
I = imread('checkerboard/right_eval01.png');
I = rgb2gray(I);
[r, c] = size(I);
for i = 1:r
    X(i, :) = fft(I(i,:));
end

for j = 1:c
    Y(:,j) = fft(X(:,j));
end
figure(1)
imshow(I)
figure(2)
M = Y;
M = fftshift(M);
Ab = abs(M);
Ab = ((Ab - min(min(Ab)))/(max(max(Ab)))).*255;
imshow(Ab)