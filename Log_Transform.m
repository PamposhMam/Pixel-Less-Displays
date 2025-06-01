clc; clear; close all;

% Define the size of the hologram plane in pixels
num_rows = 1280; % Vertical resolution
num_cols = 1024; % Horizontal resolution

% Create the spatial coordinate arrays
u = linspace(-num_cols / 2, num_cols / 2, num_cols); % Horizontal pixel indices
v = linspace(-num_rows / 2, num_rows / 2, num_rows); % Vertical pixel indices

% Create the 2D grid
[U, V] = meshgrid(u, v);

% Convert Cartesian coordinates to polar coordinates
Theta = atan2(V, U); % Compute theta for each (U, V)
Theta(Theta < 0) = Theta(Theta < 0) + 2*pi; % Ensure theta is in [0, 2π]

% Define fixed parameters
rho = sqrt(U.^2 + V.^2); % Fixed radius
theta1_fixed = pi/2; % Max theta

% Limit Theta to [0, theta1_fixed]
Theta(Theta > theta1_fixed) = 0; % Set out-of-range values to 0

% Compute J1(theta) dynamically
J1_value = besselj(1, Theta); % J1 as a function of Theta

% Compute the aperture function
aperture = J1_value ./ rho; % Aperture varies with theta

% Normalize aperture for visualization
aperture_normalized = mat2gray(aperture); 

% Display the aperture function
figure;
imshow(aperture_normalized, 'XData', u, 'YData', v, 'InitialMagnification', 'fit');
title('Aperture with J1(Theta)');
xlabel('X');
ylabel('Y');
colormap hot;
colorbar;
imwrite(aperture_normalized, 'aperture_image.png'); % Save image

% Compute the 2D Fourier Transform (Hologram)
aperture_ft = fftshift(fft2(aperture)); % Compute FFT and shift zero frequency to center
aperture_ft_mag = abs(aperture_ft); % Get magnitude

% Normalize hologram for visualization
hologram_normalized = mat2gray(aperture_ft_mag);

% Display the hologram
figure;
imshow(hologram_normalized, 'XData', u, 'YData', v, 'InitialMagnification', 'fit');
title('Hologram with J1(Theta)');
xlabel('X');
ylabel('Y');
colormap hot;
colorbar;
imwrite(hologram_normalized, 'hologram_image.png'); % Save image
