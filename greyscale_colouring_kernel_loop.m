% MultiKernelPixelHolo.m
% Generate holograms from edge pixels using different kernel sizes

clear; clc;

% --- Load grayscale image ---
filename = 'Biccy.jpg';
[~, name, ~] = fileparts(filename);  % name = 'Biccy'
img = im2double(imread(filename));
if size(img, 3) > 1
    img = rgb2gray(img);
end
[M, N] = size(img);


% --- Create output folder ---
outdir = 'HologramOutputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% --- Edge detection (original size) ---
edges = edge(img, 'Canny');
[rowIdx, colIdx] = find(edges);  % Edge coordinates

% --- Kernel sizes to loop over ---
kernel_sizes = 1:10;

% --- Loop over each kernel size ---
for kSize = kernel_sizes

    fprintf('Processing kernel size %dx%d...\n', kSize, kSize);
    pad = ceil(kSize / 2);

    % --- Pad image and edge map ---
    padded_img = padarray(img, [pad pad], 0, 'both');
    padded_edges = padarray(edges, [pad pad], 0, 'both');

    % --- Frequency grid (same size as original image) ---
    [u, v] = meshgrid(linspace(-0.5, 0.5, N), linspace(-0.5, 0.5, M));
    U = u;
    V = v;

    % --- Initialize hologram ---
    H_total = zeros(M, N);
    count = 0;

    % --- Loop through edge pixels ---
    for idx = 1:length(rowIdx)
        y = rowIdx(idx);
        x = colIdx(idx);

        % Adjust for padding
        y_pad = y + pad;
        x_pad = x + pad;

        % Grayscale intensity
        amp = padded_img(y_pad, x_pad);

        % Rectangle around (x,y) with kernel size
        half = floor(kSize / 2);
        coords = [colIdx(idx) - half, rowIdx(idx) - half;
          colIdx(idx) + half, rowIdx(idx) + half];  % CORRECT


        % Clip coords to valid image bounds
        coords(:,1) = max(min(coords(:,1), N), 1);
        coords(:,2) = max(min(coords(:,2), M), 1);

        % Generate hologram
        H = RECT_colour(amp, coords, U, V);
        H_total = H_total + H;

        count = count + 1;
        if mod(count, 400) == 1
            fprintf('  Processed %d pixels\n', count);
        end
    end

    % --- Compute replay and normalize ---
    replay = ifft2(H_total);
    replay_abs = mat2gray(abs(replay));
    imshow(replay_abs);
    
    % --- Save replay image ---
    replay_filename = fullfile(outdir, sprintf('%s_replay_kernel%dx%d.bmp', name, kSize, kSize));
    imwrite(replay_abs, replay_filename);
    
    % --- Save normalized hologram magnitude ---
    Holo = mat2gray(abs(H_total));
    holo_filename = fullfile(outdir, sprintf('%s_hologram_kernel%dx%d.bmp', name, kSize, kSize));
    imwrite(Holo, holo_filename);


end

disp('All kernel sizes processed and saved.');
