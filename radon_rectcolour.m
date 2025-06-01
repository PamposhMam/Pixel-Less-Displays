clear all; close all; clc;

% Grid size
Ny = 1024;
Nx = 1280;

% Frequency grid (normalized)
[u_grid, v_grid] = meshgrid(...
    linspace(-0.5, 0.5, Nx), ...
    linspace(-0.5, 0.5, Ny));

% Define square aperture coordinates
coords = [400, 500;
          400, 600;
          500, 600;
          500, 500];  % [y, x]

% Convert to bounding box format
bbox = [min(coords(:,2)), min(coords(:,1));  % x, y
        max(coords(:,2)), max(coords(:,1))];

% Amplitude
amp = 1;

% Generate hologram (Fourier domain)
Holo_FT = RECT_colour(amp, bbox, u_grid, v_grid);

% Inverse Fourier transform to reconstruct field
replay = ifft2(Holo_FT);
replay_abs = abs(replay);

% Normalize for display
replay_img = mat2gray(replay_abs);

% Display reconstructed field
figure;
imshow(replay_img, []);
title('Reconstructed Replay Intensity');

% Apply Radon transform
theta = 0:1:179;  % angles in degrees
[R, xp] = radon(replay_abs, theta);

% Display Radon projection
figure;
imagesc(theta, xp, R); colormap hot;
xlabel('Angle (degrees)'); ylabel('Projection Position');
title('Radon Transform of Replay Field');

% Optional: inverse Radon reconstruction
recon = iradon(R, theta, 'linear', 'Ram-Lak', 1, Ny);

% Display reconstruction
figure;
imagesc(recon); axis image; colormap hot;
title('Reconstruction via Inverse Radon Transform');
