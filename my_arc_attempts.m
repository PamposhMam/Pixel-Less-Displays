% Define the size of the hologram plane in pixels
num_rows = 780; % Number of pixels along the vertical dimension
num_cols = 1024; % Number of pixels along the horizontal dimension

% Create the spatial coordinate arrays
u = linspace(-num_cols / 2, num_cols / 2, num_cols); % Horizontal pixel indices
v = linspace(-num_rows / 2, num_rows / 2, num_rows); % Vertical pixel indices

% Create the 2D grid
[U, V] = meshgrid(u, v);

% Set the length of the arc, r, and the angle of rotation, theta1
r = 10;
theta1 = pi;

% Compute polar coordinates (rho, phi) from (U, V)
rho = U; %sqrt(U.^2 + V.^2);
phi = V; %atan2(V, U);

% Define alpha as a complex exponential
alpha = 1j * rho .* cos(phi - theta1);

% Evaluate expressions at r and 0
var_r_r = r; % Upper limit of the integral
var_r_0 = 0; % Lower limit of the integral

%% The Maths:

L_1_r = (theta1^2) .* exp(-var_r_r .* alpha) / 2;
L_2_r = var_r_r.^2 ./ alpha;
L_3_r = 2 * var_r_r ./ (alpha.^2);
L_4_r = 2 ./ (alpha.^3);
L_r = L_1_r .* (L_2_r - L_3_r + L_4_r); % L evaluated at r

L_1_0 = (theta1^2) .* exp(-var_r_0 .* alpha) / 2;
L_2_0 = var_r_0.^2 ./ alpha;
L_3_0 = 2 * var_r_0 ./ (alpha.^2);
L_4_0 = 2 ./ (alpha.^3);
L_0 = L_1_0 .* (L_2_0 - L_3_0 + L_4_0); % L evaluated at 0

% Compute L_final as the difference L(r) - L(0)
L_final = L_r - L_0;

% Normalize L to [0, 1] for image display
L_normalized = mat2gray(abs(L_final));  

% Display the image
imshow(L_normalized);

% Save the image to a file (e.g., PNG)
imwrite(L_normalized, 'L_image.png');

% Compute the replay field using inverse FFT
replay_field = ifft2(L_final);  % Apply inverse FFT

% Normalize the replay field to [0, 1] for display
replay_field_normalized = mat2gray(abs(replay_field));

% Display the replay field
imshow(replay_field_normalized);

% Save the replay field as an image
imwrite(replay_field_normalized, 'replay_field_image.png');
