% Define the size of the hologram plane in pixels
num_rows = 780; % Number of pixels along the vertical dimension
num_cols = 1024;  % Number of pixels along the horizontal dimension

% Create the spatial coordinate arrays
u = linspace(-num_cols / 2, num_cols / 2, num_cols); % Horizontal pixel indices
v = linspace(-num_rows / 2, num_rows / 2, num_rows); % Vertical pixel indices

% Create the 2D grid
[U, V] = meshgrid(u, v);
a= num_cols/2;         %offset, x-axis
b= num_rows/2;         %offset, y-axis
theta0=0;        %starting angle
theta1=pi/10;        %ending angle
r=15  ;           %length of the arc

% Define L_num as a complex exponential for a range (-2, 0) -> (2, 0)
x1 =  a+ r*cos(theta1); x0 = a+r*cos(theta0); % Define horizontal points for the range
y1 =  b+r*sin(theta1); y0 = b+r*sin(theta0);  % Define vertical points for the range

% Compute L_num1 using complex exponential
L_num1 = 1j * exp(1j * 2 * pi * (U * x0 + V * y0));  % Using U and V coordinates

% Define the denominator term L_denom
L_denom = 2 * pi * ((x1 - x0) * U + (y1 - y0) * V);

% Compute L_num2
L_num2 = 1 - exp(1j * 2 * pi * ((x1 - x0) * U + (y1 - y0) * V));

angles=linspace(theta1, theta0, 1000);

% Compute L as the grid, using element-wise operations
L = (L_num1 .* L_num2) ./ L_denom;

% Normalize L to [0, 1] for image display (if necessary)
L_normalized = mat2gray(abs(L));  % Normalize using absolute value

% Display the image
imshow(L_normalized);

% Save the image to a file (e.g., PNG)
imwrite(L_normalized, 'L_image.png');

% Example: Compute the replay field (ifft2 is just an example here)
replay_field = ifft2(L);  % Apply inverse FFT (or any other computation)

% Normalize the replay field to [0, 1] for display
replay_field_normalized = mat2gray(abs(replay_field));

% Display the replay field
imshow(replay_field_normalized);

% Save the replay field as an image
imwrite(replay_field_normalized, 'replay_field_image.png');
