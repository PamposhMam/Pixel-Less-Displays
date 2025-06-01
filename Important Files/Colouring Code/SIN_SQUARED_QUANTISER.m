% QuantizerComparisonHolo_SinSquared.m
% Generate a single hologram using sin^2 quantization, and save the hologram matrix

clear; clc;

% --- Load grayscale image ---
filename = 'Biccy.jpg';
[~, name, ~] = fileparts(filename);
img = im2double(imread(filename));
if size(img, 3) > 1
    img = rgb2gray(img);
end

% --- Embed image onto top half of canvas with double the number of rows ---
%T_total = zeros(size(img,1)*2, size(img,2));  % Double rows only
%T_total(1:size(img,1), :) = img;              % Place original image at the top
%img = T_total;

% --- Resize to 1024x1280 using bilinear interpolation ---
%target_height = 1024;
%target_width = 1280;
%img = imresize(img, [target_height, target_width], 'bilinear');

% --- Image size AFTER resizing ---
[M, N] = size(img);
total_pixels = M * N;

% --- Frequency grid ---
[u, v] = meshgrid(linspace(-0.5, 0.5, N), linspace(-0.5, 0.5, M));
U = u;
V = v;

% --- Create output folder ---
outdir = 'SinSquared_Only_Output';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% --- Apply sin-squared quantizer ---
nLevels = 128;
mapped = (sin(pi * img / 2)).^2;
quantized = floor(mapped * nLevels);
quantized = min(quantized, nLevels - 1);

% --- Build hologram from quantized regions ---
H_partial = zeros(M, N);
for i = 0:(nLevels - 1)
    mask = quantized == i;
    if ~any(mask(:)), continue; end

    CC = bwconncomp(mask, 8);
    for regionIdx = 1:CC.NumObjects
        pixelList = CC.PixelIdxList{regionIdx};
        [rows, cols] = ind2sub([M, N], pixelList);

        if numel(pixelList) < 3
            continue;
        end

        coords = [min(cols), min(rows); max(cols), max(rows)];
        mean_intensity = mean(img(pixelList));
        scale = 1 + numel(pixelList) / total_pixels;
        amp = mean_intensity * scale;

        H = RECT_colour(amp, coords, U, V);
        H_partial = H_partial + H;
    end
end

% --- Save hologram matrix to .mat file ---
save(fullfile(outdir, [name '_Holo_part.mat']), 'H_partial');

% --- Optional: visualize replay field
Replay = ifft2(H_partial);
imwrite(mat2gray(abs(Replay)), fullfile(outdir, [name '_replay.bmp']));

disp('✅ Hologram (H_partial) generated and saved.');
