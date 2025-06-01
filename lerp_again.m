% Define the hologram resolution
alias = 1;
Holo = zeros(alias * 500, alias * 500);
[numrows, numcols] = size(Holo);

% Define parameters
omega = 1;
r1 = 150;
r2 =190;

% x1 = r1 * cos(0 + omega)* alias; % ax
% y1 = r1 * sin(0 + omega) * alias; % ay
% x0 = r2 * cos(0) * alias;          % mx
% y0 = r2 * sin(0) * alias;          % my
% x2 = r1 * cos(0 + 2*omega) * alias; % bx
% y2 = r1 * sin(0 + 2*omega) * alias; % by


x0 = 10;          % ax
y0 = 10;          % ay
x1 = 180; % mx
y1 = 130; % my
x2 = 210; % bx
y2 = 150; % by

n = 100;

% Output folder setup
name = sprintf("lerp_x0_%d_y0_%d_x1_%d_y1_%d_x2_%d_y2_%d", x0, y0, x1, y1, x2, y2);
outputFolder = fullfile(pwd, name);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Define the U and V axes
U = linspace(-0.5, 0.5, numcols);
V = linspace(-0.5, 0.5, numrows);
[U, V] = meshgrid(U, V);

% Simpson's Rule setup
t_values = linspace(0, 1, 2*n + 1); % 2n + 1 points

% Define the CLERP rational expression (corrected)
f = @(t) exp(2 * pi * 1j * ( ...
    U .* ( ((1 - t) .* (x2 - x1) .* x0 + t .* (x1 - x0) .* x2) ./ ...
           ((1 - t) .* (x2 - x1) + t .* (x1 - x0)) ) + ...
    V .* ( ((1 - t) .* (y2 - y1) .* y0 + t .* (y1 - y0) .* y2) ./ ...
           ((1 - t) .* (y2 - y1) + t .* (y1 - y0)) ) ));

% Use the external Simpson's Rule integrator
Holo_Func = Simpsons_Rule(f, 0, 1, n);

% Combine result with the original hologram
Holo = Holo + Holo_Func;

% Compute the Replay Field
Replay = fft2(Holo);

% Define the Jones Matrix
tau = pi;
thetajones = pi/2;
J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, ...
     -1j*sin(tau/2)*sin(2*thetajones); ...
     -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];

% Apply the Jones Matrix to the hologram
Ex_mod = J(1,1) * Holo + J(1,2) * Holo;
Ey_mod = J(2,1) * Holo + J(2,2) * Holo;

% Final hologram modulation
Holo_mod = Ex_mod + Ey_mod;
Holo_mod = mod(angle(Holo_mod), pi) > (pi/2); % Binary phase
Holo_mod = Holo_mod * pi;

% Display results
figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or \pi)');
figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT) from Binary Phase');

% Save outputs
imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "Binary_Phase_Hologram.bmp"));
imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "Binary_Replay.bmp"));

% Resize and save hologram
Holo_mod_sized = imresize(Holo_mod, [1024, 1280], 'bilinear');
Holo_mod_sized = fliplr(Holo_mod_sized);
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.bmp"));
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.png"));