%Using the sin and cos transforms from the databook

% rho= 5;
% angles= [0, pi/3];  %w0, w1
% const= 2*pi*rho/1j; %pre-multiplier
% 
% eqn=zeros(2,1);
% for i=1:2
%     eqn= const*dirac(angles(i))/angles(i);
% end
% 
% Hologram= eqn(2)-eqn(1);
% Shift it to the centre
% Hologram=Holohram*exp(-2*1j*pi*(num_rows*v/2+num_cols*u/2))
%%
clc; clear; close all;

% Define the size of the hologram plane in pixels
num_rows = 1080;  % Height
num_cols = 1280;  % Width

% Create spatial coordinate arrays
u = linspace(-num_cols / 2, num_cols / 2, num_cols); % Horizontal pixel indices
v = linspace(-num_rows / 2, num_rows / 2, num_rows); % Vertical pixel indices

% Create the 2D grid
[U, V] = meshgrid(u, v);

% Hologram parameters
rho = 5;
angles = [0, pi/3];  % w0, w1
const = (2 * pi * rho) / 1j;  % Pre-multiplier

% Initialize the hologram
H = zeros(num_rows, num_cols);

for i = 1:2
    if angles(i) == 0
        eqn(i) = 0; % Handling dirac(0), which ideally should be infinite, but we approximate it as 0
    else if angles(i)
        eqn(i) = const / angles(i); % Approximating dirac as a division by the angle
    end
end

% Compute the hologram
Hologram = eqn(2) - eqn(1);

% Normalize hologram for display
H_normalized = mat2gray(abs(H));

% Display the hologram on a black grid
figure;
imshow(H_normalized);
title('Hologram on Black Grid');

% Save the hologram image
imwrite(H_normalized, 'hologram.png');

% Compute the Fourier Transform (Replay Image)
replay_field = fftshift(fft2(H));
replay_intensity = abs(replay_field);

% Normalize the replay field for display
replay_normalized = mat2gray(replay_intensity);

% Display the replay image
figure;
imshow(replay_normalized);
title('Fourier Transform (Replay Image)');

% Save the replay image
imwrite(replay_normalized, 'replay_image.png');
