clear all; close all;

% --- Define two rectangular shapes in x,y domain ---
coords1 = [120, 650;
           120, 700;
           200, 700;
           200, 120] / 2;

coords2 = [100, 650;
           100, 800;
           300, 800;
           300, 650] / 2;

coords_list = {coords1, coords2};
amp_list = [10, 8];

% --- Grid parameters ---
Ny = 1024;
Nx = 1280;

% Spatial frequency domain (u, v), unshifted
[u_grid, v_grid] = meshgrid(linspace(0, 1, Nx), linspace(0, 1, Ny));

% Hologram accumulator
H_total = zeros(Ny, Nx);

% Delta position (can vary per shape if you want)
decu = -0.5;
decv = -0.5;

% --- Loop over each rectangle and call RectXPONential ---
for i = 1:length(coords_list)
    coords = coords_list{i};
    amp = amp_list(i);
    
    % Use your original function as-is
    H_component = RectXPONential(amp, coords, u_grid, v_grid, decu, decv);

    % Accumulate
    H_total = H_total + H_component;
end

% --- Replay field ---
Replay = ifft2(H_total);  % NO fftshift

% --- Visualization ---
figure;
imagesc(abs(H_total)); axis image; colormap hot;
title('Superposed Hologram |H(u,v)|');

figure;
imshow(mat2gray(abs(Replay))); title('Replay Field |ifft2(H)|');

figure;
imshow(angle(Replay), []); title('Replay Field Phase');
