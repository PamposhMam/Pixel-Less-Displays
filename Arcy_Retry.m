% Define the size of the hologram plane in pixels
num_rows = 780; % Number of pixels along the vertical dimension
num_cols = 1024;  % Number of pixels along the horizontal dimension

% Create the spatial coordinate arrays
u = linspace(-num_cols / 2, num_cols / 2, num_cols); % Horizontal pixel indices
v = linspace(-num_rows / 2, num_rows / 2, num_rows); % Vertical pixel indices

% Create the 2D grid
[U, V] = meshgrid(u, v);

% Parameters
a = num_cols / 2; % Offset, x-axis
b = num_rows / 2; % Offset, y-axis
theta1 = pi;
theta0 = 0.5;

r1 = 15; % Length of the arc
r0 = 0;  % Starting point of the arc
alpha = [r1, r0]; % Define alpha as a vector

% Compute ftr2
ftr2 = 2 * pi * 1j * alpha;

% Compute L1, L2, L3, L4 using element-wise operations
L1 = (alpha(1) * exp(ftr2(1) .* U)) ./ (ftr2(1) .* U);
L2 = (-exp(ftr2(1) .* U)) ./ ((ftr2(1) .* U).^2);
L3 = (exp(ftr2(1) .* theta1 .* V)) ./ (((ftr2(1) .* U).^2) * theta1);
L4 = (-exp(ftr2(1) .* theta0 .* V)) ./ (((ftr2(1) .* U).^2) * theta0);

% Compute L as the grid
L_hold = L1 + L2 + L3 + L4;
L = L_hold; % Simplify the computation for this context

% Normalize L to [0, 1] for image display (if necessary)
L_normalized = mat2gray(abs(L)); % Normalize using absolute value

% Display the hologram image
imshow(L_normalized);
title('Hologram');

% Save the hologram image to a file (e.g., PNG)
imwrite(L_normalized, 'L_image.png');

% Example: Compute the replay field (ifft2 is just an example here)
replay_field = ifft2(L); % Apply inverse FFT

% Normalize the replay field to [0, 1] for display
replay_field_normalized = mat2gray(abs(replay_field));

% Display the replay field
figure;
imshow(replay_field_normalized);
title('Replay Field');

% Save the replay field as an image
imwrite(replay_field_normalized, 'replay_field_image.png');
