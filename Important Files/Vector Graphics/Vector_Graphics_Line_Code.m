% % Define the size of the hologram plane in pixels
% num_rows = 780; % Number of pixels along the vertical dimension
% num_cols = 1024;  % Number of pixels along the horizontal dimension
% 
% % Create the spatial coordinate arrays
% u = linspace(-num_cols / 2, num_cols / 2, num_cols); % Horizontal pixel indices
% v = linspace(-num_rows / 2, num_rows / 2, num_rows); % Vertical pixel indices
% 
% % Create the 2D grid
% [U, V] = meshgrid(u, v);
% 
% % Define L_num as a complex exponential for a range (-2, 0) -> (2, 0)
% x1 = num_cols/2 +100; x0 = num_cols/2 -100; % Define horizontal points for the range
% y1 = num_rows/2; y0 = num_rows/2;  % Define vertical points for the range
% 
% % Compute L_num1 using complex exponential
% L_num1 = 1j * exp(1j * 2 * pi * (U * x0 + V * y0));  % Using U and V coordinates
% 
% % Define the denominator term L_denom
% L_denom = 2 * pi * ((x1 - x0) * U + (y1 - y0) * V);
% 
% % Compute L_num2
% L_num2 = 1 - exp(1j * 2 * pi * ((x1 - x0) * U + (y1 - y0) * V));
% 
% % Compute L as the grid, using element-wise operations
% L = (L_num1 .* L_num2) ./ L_denom;
% 
% % Normalize L to [0, 1] for image display (if necessary)
% L_normalized = mat2gray(abs(L));  % Normalize using absolute value
% 
% % Display the image
% imshow(L_normalized);
% 
% % Save the image to a file (e.g., PNG)
% imwrite(L_normalized, 'L_image.png');
% 
% % Example: Compute the replay field (ifft2 is just an example here)
% replay_field = ifft2(L);  % Apply inverse FFT (or any other computation)
% 
% % Normalize the replay field to [0, 1] for display
% replay_field_normalized = mat2gray(abs(replay_field));
% 
% % Display the replay field
% imshow(replay_field_normalized);
% 
% % Save the replay field as an image
% imwrite(replay_field_normalized, 'replay_field_image.png');

%% For multiple Points
clear all

shape= "Star";
outputFolder= shape;

% Step 2: Create a Folder with the Same Name
outputFolder = fullfile(pwd, shape); % Create folder path
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % Make directory if it doesn't exist
end

% Define the size of the hologram plane in pixels
num_rows = 1280; % Number of pixels along the vertical dimension
num_cols = 1024;  % Number of pixels along the horizontal dimension

% Create the spatial coordinate arrays
u = linspace(-num_cols / 2, num_cols / 2, num_cols); % Horizontal pixel indices
v = linspace(-num_rows / 2, num_rows / 2, num_rows); % Vertical pixel indices

% Create the 2D grid
[U, V] = meshgrid(u, v);


%Might need a scale factor
alpha=30;
alpha2=2*alpha;

% Define coordinates for multiple points
% Fixed starting points (bottom-left corner of the square)
a0= num_cols/(2*alpha);
b0= num_rows/(2*alpha);

%coordinates for the shape
x1= alpha*[a0 + 6*cos(pi/10), a0 + 2*cos(3*pi/10), a0 + 6*cos(pi/2), a0 + 2*cos(7*pi/10), a0 + 6*cos(9*pi/10), a0 + 2*cos(11*pi/10), a0 + 6*cos(13*pi/10),a0 + 2*cos(15*pi/10), a0 + 6*cos(17*pi/10), a0 + 2*cos(19*pi/10), a0+ 6*cos(21*pi/10)];
y1= alpha*[b0 + 6*sin(pi/10), b0 + 2*sin(3*pi/10), b0 + 6*sin(pi/2), b0 + 2*sin(7*pi/10), b0 + 6*sin(9*pi/10), b0 + 2*sin(11*pi/10), b0 + 6*sin(13*pi/10),b0 + 2*sin(15*pi/10), b0 + 6*sin(17*pi/10), b0 + 2*sin(19*pi/10), b0+ 6*sin(21*pi/10)];

% % Ensure the starting points for each segment
 x0 = [a0 + 2*cos(19*pi/10), x1(1:end-1)]; % Start of each segment
 y0 = [b0 + 2*sin(19*pi/10), y1(1:end-1)]; % Start of each segment


% Ensure the size of x1 and y1 match
%assert(numel(x1) == numel(y1), 'x1 and y1 must have the same number of elements.');

% Initialize the output grid L
L = zeros(size(U));
% Initialize the grid to store the accumulated contributions
L_accumulated = zeros(size(U));

% Loop through each pair of (x1, y1)
for k = 1:numel(x1)
    % Compute L_num1 using complex exponential
    L_num1 = 1j * exp(1j * 2 * pi * (U * x0(k) + V * y0(k)));

    % Define the denominator term L_denom
    L_denom = 2 * pi * ((x1(k) - x0(k)) * U + (y1(k) - y0(k)) * V);
    epsilon = 1e-9; % Small value to prevent division by zero
    L_denom = L_denom + epsilon;

    % Compute L_num2
    L_num2 = 1 - exp(1j * 2 * pi * ((x1(k) - x0(k)) * U + (y1(k) - y0(k)) * V));

    % Compute the hologram for this point
    L_temp = (L_num1 .* L_num2) ./ L_denom;

    % Accumulate the holograms
    L_accumulated = L_accumulated + L_temp;

    % Debugging output
    disp(['Iteration ', num2str(k)]);
    disp(['Mean of L_temp (iteration ', num2str(k), '): ', num2str(mean(abs(L_temp(:))))]);
end

% Compute the average of the accumulated holograms
L = L_accumulated / numel(x1);

% Normalize L to [0, 1] for image display
L_normalized = mat2gray(abs(L));  % Normalize using absolute value

% Display the image
imshow(L_normalized);

% Save the image to a file (e.g., PNG)
imwrite(L_normalized, 'L_star_holo.png');
imwrite(L_normalized, 'L_star_holo.bmp');

% Example: Compute the replay field (ifft2 is just an example here)
replay_field = ifft2(L);  % Apply inverse FFT (or any other computation)

% Normalize the replay field to [0, 1] for display
replay_field_normalized = mat2gray(abs(replay_field));

% Display the replay field
imshow(replay_field_normalized);

% Save the replay field as an image
imwrite(replay_field_normalized, 'replay_field_image_star.png');


%%
% Create an empty image (black) of the same size as the hologram plane
target_image = zeros(num_rows, num_cols);

% Number of points to interpolate for each line
num_points = 1000; % Adjusted to a reasonable value for efficiency

% Loop through each pair of (x0, y0) and (x1, y1)
for k = 1:numel(x0)
    % Generate evenly spaced points between (x0(k), y0(k)) and (x1(k), y1(k))
    x_points = linspace(x0(k), x1(k), num_points);
    y_points = linspace(y0(k), y1(k), num_points);

    % Convert continuous coordinates to pixel indices
    x_pixels = round(x_points); % Adjust for image center
    y_pixels = round(y_points); % Adjust for image center

    % Filter valid indices (stay within image bounds)
    valid_indices = (x_pixels >= 1 & x_pixels <= num_cols) & ...
                    (y_pixels >= 1 & y_pixels <= num_rows);

    % Apply valid indices to filter out-of-bound points
    x_pixels = x_pixels(valid_indices);
    y_pixels = y_pixels(valid_indices);

    % Mark the points on the target image
    linear_indices = sub2ind(size(target_image), y_pixels, x_pixels);
    target_image(linear_indices) = 1; % Set pixels to white (1)
end

% Display the target image
imshow(target_image);

% Save the target image as a file (e.g., BMP)
imwrite(target_image, 'target_image_star.bmp');


%%
% Step 1: Generate a Phase Mask for "JONES APERTURE"
textSize = 50; % Adjust for visibility
position = [round(num_cols * 0.3), round(num_rows * 0.5)]; % Center it

% Create a blank mask to overlay text
textMask = ones(num_rows, num_cols); 
textMask = insertText(textMask, position, 'JONES APERTURE', 'FontSize', textSize, 'BoxColor', 'black', 'TextColor', 'white');
textMask = im2bw(rgb2gray(textMask)); % Convert to binary mask

% Convert binary mask into a phase modulation (0 or pi)
phaseModulation = pi * textMask; 

% Step 2: Apply the Phase Modulation to the Hologram
Holo = L .* exp(1j * phaseModulation); 

% Step 3: Define the Jones Matrix
tau = pi;
thetajones = pi/2;
J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
    -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];

% Step 1: Apply the Jones Matrix Directly to Holo
Ex_mod = J(1,1) * Holo + J(1,2) * Holo;
Ey_mod = J(2,1) * Holo + J(2,2) * Holo;

% Step 2: Compute the Final Hologram (Before FFT)
Holo_mod = Ex_mod + Ey_mod;

% Step 3: Convert to Binary Phase (0 or pi)
Holo_mod = mod(angle(Holo_mod), pi) > (pi/2);  % Threshold to 0 or 1
Holo_mod = Holo_mod * pi;  % Convert binary 0,1 → 0,π

% Step 4: Compute the Replay Field
Replay = ifft2(exp(1j * Holo_mod)); % FFT of binary phase modulated hologram

% Step 5: Display Results
figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or π)');
figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT) from Binary Phase');

% Step 6: Save the Binary Hologram
imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "Binary_Phase_Hologram.bmp"));
imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "Binary_Replay.bmp"));

%Step 7: Also save a correct-sized hologram to the folder
Holo_mod_sized= imresize(Holo_mod, [1024, 1280], 'bilinear');
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.bmp"));