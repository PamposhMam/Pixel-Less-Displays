% Define the hologram resolution
alias = 2;
Holo = zeros(alias * 1080, alias * 1280);
[numrows, numcols] = size(Holo);

% Define parameters
omega = pi;
r1 = 60;
r2 = 60;

x0 = r1 * cos(0 + omega)* alias;
y0 = r1 * sin(0 + omega) * alias;
x1 = r1 * cos(0) * alias;
y1 = r1 * sin(0) * alias;

%%
name= num2str(r1) + "_" + num2str(r2) + "_" + num2str(omega);

% Step 2: Create a Folder with the Same Name
outputFolder = fullfile(pwd, name); % Create folder path
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % Make directory if it doesn't exist
end

%%

% Define the U and V axes
U = linspace(-0.5, 0.5, numcols); %-1,1 originally
V = linspace(-0.5, 0.5, numrows); %-1,1 originally
[U, V] = meshgrid(U, V);

% Number of steps for Simpson's Rule
n = 350; % Adjust for accuracy
num_intervals = 2 * n;
h = 1 / num_intervals; % Step size

% Precompute alpha and beta for all required t values
t_values = linspace(0, 1, num_intervals + 1); % 2n + 1 points
alpha_values = (sin((1 - t_values) * omega) * x0 + sin(t_values * omega) * x1) / sin(omega);
beta_values = (sin((1 - t_values) * omega) * y0 + sin(t_values * omega) * y1) / sin(omega);

% Precompute Simpson's weights
weights = 2 * ones(1, num_intervals + 1);
weights(1) = 1; % First weight
weights(end) = 1; % Last weight
weights(2:2:end-1) = 4; % Odd indices get weight 4

% Compute the hologram integral using precomputed values
Holo_Func = zeros(size(Holo)); % Initialize sum

for i = 1:length(t_values)
    % Compute the phase term for all pixels at once
    Holo_Deriv = exp(2 * pi * 1j * (U .* alpha_values(i) + V .* beta_values(i)));
    
    % Apply Simpson's weight and sum
    Holo_Func = Holo_Func + weights(i) * Holo_Deriv;
end

% Final integral calculation
Holo_Func = (h / 3) * Holo_Func;

% Combine result with the original hologram
Holo = Holo + Holo_Func;

% Save the result
imwrite(mat2gray(abs(Holo)), 'Arc_Hologram.bmp'); % Save the computed hologram

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

% Step 4: Compute the Replay Field
Replay = fftshift(fft2(Holo)); % FFT of binary phase modulated hologram

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

% %% --- Encoding for the hologram ---
% 
% % Step 1: Convert the hologram to binary phase (0 or π) based on phase angle
% Holo_phase = angle(Holo);
% Holo_binary = uint8(mod(Holo_phase, pi) > pi/2);  % 0 or 1
% 
% % Step 2: Resize to projector resolution (e.g., 1024×1280)
% Holo_binary_resized = imresize(double(Holo_binary), [1024, 1280], 'bilinear');
% 
% % Step 3: Convert to 8-bit (0 = 0, 255 = π)
% Holo_binary_uint8 = uint8(round(Holo_binary_resized * 255));
% 
% % Step 4: Optional flip (e.g., horizontal for SLM alignment)
% Holo_binary_uint8 = fliplr(Holo_binary_uint8);
% 
% % Step 5: Save final binary phase hologram
% imwrite(Holo_binary_uint8, fullfile(outputFolder, "Holo_Sized.bmp"));
% imwrite(Holo_binary_uint8, fullfile(outputFolder, "Holo_Sized.png"));
% 
% % Optional: simulate replay (if needed)
% Replay = fftshift(fft2(exp(1j * pi * double(Holo_binary))));  % Use original binary map
% imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "Binary_Replay.bmp"));







%% Alt Save Everything Code:
% Display the image
figure; imshow(Holo); title('Hologram');

% Save the image
imwrite(Holo, fullfile(outputFolder, "_holo.png"));
imwrite(Holo, fullfile(outputFolder, "_holo.bmp"));

% Compute and save the replay field
replay_field = fft2(Holo);
replay_field_normalized = mat2gray(abs(replay_field));

figure; imshow(replay_field_normalized); title('Reconstructed Image');
imwrite(replay_field_normalized, fullfile(outputFolder, "_replay_field.png"));

%% Step 2: Overlay "JONES APERTURE"
textSize = 50;
position = [round(grid_cols * 0.3), round(grid_rows * 0.5)];

% Create a blank mask to overlay text
textMask = ones(grid_rows, grid_cols); 
textMask = insertText(textMask, position, 'JONES APERTURE', 'FontSize', textSize, 'BoxColor', 'black', 'TextColor', 'white');
textMask = im2bw(rgb2gray(textMask)); % Convert to binary mask

% Convert to phase modulation (0 or pi)
phaseModulation = pi * textMask; 

% Apply the phase modulation to the hologram
Holo = Holo .* exp(1j * phaseModulation);

%% Step 3: Apply the Jones Matrix
tau = pi;
thetajones = pi/2;
J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
    -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];

% Apply the Jones Matrix
Ex_mod = J(1,1) * Holo + J(1,2) * Holo;
Ey_mod = J(2,1) * Holo + J(2,2) * Holo;

% Compute the final hologram before FFT
Holo_mod = Ex_mod + Ey_mod;

% Convert to binary phase (0 or π)
Holo_mod = mod(angle(Holo_mod), pi) > (pi/2);
Holo_mod = Holo_mod * pi;

% Compute the Replay Field
Replay = ifft2(exp(1j * Holo_mod));

% Display Results
figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or π)');
figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image from Binary Phase');

% Save the Binary Hologram
imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "_Binary_Phase_Hologram.bmp"));
imwrite(mat2gray(abs(Replay)), fullfile(outputFolder,  "_Binary_Replay.bmp"));

% Step 7: Also save a correctly sized hologram
Holo_mod_sized = imresize(Holo_mod, [1024, 1280], 'bilinear');
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "_Holo_Sized.bmp"));