% PixelHoloFromGreyscale.m
% Generate hologram from every edge pixel, modulated by grayscale intensity

clear; clc;

% --- Load grayscale image ---
filename = 'Biccy.jpg';
[~, name, ~] = fileparts(filename);  % name = 'Biccy'
img = im2double(imread('Biccy.jpg'));
if size(img,3) > 1
    img = rgb2gray(img);
end
[M, N] = size(img);

% --- Edge detection ---
edges = edge(img, 'Canny');  % Logical matrix
[rowIdx, colIdx] = find(edges);  % Get edge pixel coordinates

% --- Frequency grid ---
[u, v] = meshgrid(linspace(-0.5, 0.5, N), linspace(-0.5, 0.5, M));
U = u;
V = v;

% --- Initialize hologram ---
H_total = zeros(M, N);
count=0
N=[1,2,3,4,5,6,7,8,9,10];

% --- Loop through each edge pixel ---
for k = 1:length(rowIdx)
    y = rowIdx(k);
    x = colIdx(k);
    
    amp = img(y, x);  % Get grayscale intensity

    % Define a small rectangle (e.g. 2×2 region centered at the pixel)
    coords = [x-1, y-1; x+1, y+1];
    
    % Generate tiny rectangle hologram
    H = RECT_colour(amp, coords, U, V);
    
    H_total = H_total + H;
    count=count+1;
    if mod(count,400)==1
        disp(count)
    end
end

% --- Display final hologram ---
figure;
imagesc(abs(H_total));
colormap gray; axis image off;
title('Pixel-by-pixel Greyscale Hologram');
%%
replay=ifft2(H_total);
replay_abs=mat2gray(abs(replay));
imshow(replay_abs);
imwrite(replay_abs, [name '_replay.bmp']);
Holo=mat2gray(abs(H_total));
imwrite(Holo, [name '_hologram.bmp']);
