% Load and preprocess the image
I = imread('lenna.jpg');
I_gray = rgb2gray(I);

% Apply Gaussian blur to smooth out small details
I_blur = imgaussfilt(I_gray, 4); % Increase the sigma for stronger blurring, was ogn sigma=4,biccy2=10

% Detect edges
edges = edge(I_blur, 'Canny');

% Display the result
imshow(edges);
title('Outline with Small Details Suppressed');
imwrite(edges, 'Lenna_Edge.bmp');