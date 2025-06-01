% Define the hologram resolution
alias = 1;
Holo = zeros(alias * 300, alias * 300);
[numrows, numcols] = size(Holo);

% Define parameters
omega = pi/3;
r1 = 50;
r2 = 250;

x1 = r1 * cos(0 + omega)* alias; %ax
y1 = r1 * sin(0 + omega) * alias; %ay
x0 = r2 * cos(0) * alias; %mx
y0 = r2 * sin(0) * alias; %my
x2 = r1* cos(0+2*omega)*alias; %bx
y2= r1*sin(0+2*omega)*alias; %by

% x1=100;
% y1=100;
% x0= 40;
% y0= 40;
% x2= 110;
% y2= 110;

n=100;
%%
name= "lerp" + num2str(r1) + "_" + num2str(r2) + "_" + num2str(omega);

% Step 2: Create a Folder with the Same Name
outputFolder = fullfile(pwd, name); % Create folder path
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % Make directory if it doesn't exist
end

%%

% Define the U and V axes
U = linspace(-0.5, 0.5, numcols);
V = linspace(-0.5, 0.5, numrows);
[U, V] = meshgrid(U, V);

% Number of steps for Simpson's Rule
n = 100; % Adjust for accuracy
num_intervals = 2 * n;
h = 1 / num_intervals; % Step size

% Precompute alpha and beta for all required t values
t_values = linspace(0, 1, num_intervals + 1); % 2n + 1 points
alpha_values = ((1 - (x2 - x1) * x0) .* (x1 - x0) .* x2 + (t_values) .* (x2 - x1) .* x0) ./ ...
               ((1 - x2 + x1) .* (x1 - x0) + (t_values) .* (x1 - x0));

beta_values = ((1 - (y2 - y1) * y0) .* (y1 - y0) .* y2 + (t_values) .* (y2 - y1) .* y0) ./ ...
              ((1 - y2 + y1) .* (y1 - y0) + (t_values) .* (y1 - y0));


% Define the integrand as a function handle
f = @(t) exp(2 * pi * 1j * ( ...
    U .* ( ((1 - (x2 - x1) * x0) * (x1 - x0) * x2 + t .* (x2 - x1) * x0) ./ ...
           ((1 - x2 + x1) * (x1 - x0) + t .* (x1 - x0)) ) + ...
    V .* ( ((1 - (y2 - y1) * y0) * (y1 - y0) * y2 + t .* (y2 - y1) * y0) ./ ...
           ((1 - y2 + y1) * (y1 - y0) + t .* (y1 - y0)) ) ));

% Use the external Simpson's rule function
Holo_Func = Simpsons_Rule(f, 0, 1, n);

% Combine result with the original hologram
Holo = Holo + Holo_Func;

% Step 4: Compute the Replay Field
Replay = (fft2(Holo)); % FFT of binary phase modulated hologram

% % Save the result
% imwrite(mat2gray(abs(Holo)), 'Arc_Hologram.bmp'); % Save the computed hologram

%% encoding for the hologram:
% Step 3: Define the Jones Matrix
tau = pi;
thetajones = pi/2;
J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
    -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];

% % Step 1: Apply the Jones Matrix Directly to Holo
Ex_mod = J(1,1) * Holo + J(1,2) * Holo;
Ey_mod = J(2,1) * Holo + J(2,2) * Holo;

% Step 2: Compute the Final Hologram (Before FFT)
Holo_mod = Ex_mod + Ey_mod;

% % Step 3: Convert to Binary Phase (0 or pi)
Holo_mod = mod(angle(Holo_mod), pi) > (pi/2);  % Threshold to 0 or 1
Holo_mod = Holo_mod * pi;  % Convert binary 0,1 → 0,π



% Step 5: Display Results
figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or π)');
figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT) from Binary Phase');

% Step 6: Save the Binary Hologram
imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "Binary_Phase_Hologram.bmp"));
imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "Binary_Replay.bmp"));

%Step 7: Also save a correct-sized hologram to the folder
 Holo_mod_sized= imresize(Holo_mod, [1024, 1280], 'bilinear');
% Holo_mod_

Holo_mod_sized=fliplr(Holo_mod_sized);
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.bmp"));
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.png"));

