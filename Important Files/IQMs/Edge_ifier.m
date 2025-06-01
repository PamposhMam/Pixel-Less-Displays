%Take the edge-image of some of these photos
% Load and preprocess the image
cd Images\O
I = imread('Dot full.bmp');
I_gray = rgb2gray(I);

% Apply Gaussian blur to smooth out small details
I_blur = imgaussfilt(I_gray, 4); % Increase the sigma for stronger blurring

% Detect edges
edges = edge(I_blur, 'Canny');

% Display the result
imshow(edges);
title('Outline with Small Details Suppressed');
imwrite(edges, 'O.bmp')