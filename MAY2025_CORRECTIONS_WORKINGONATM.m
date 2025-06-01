% Fully corrected version: image coordinates adjusted to (x, y) everywhere
clear all;
close all;

% Step 1: Read and Resize the Image
filename = "Biccy.bmp";
[~, name, ~] = fileparts(filename);
target = imread(filename);
Alpha= 250;
scalefac = 1;
name = name + 'FINAL' + num2str(Alpha) + "_scalefac_" + num2str(scalefac);
outputFolder = fullfile(pwd, name);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Step 2: Prepare Image and Declare Variables
target = im2bw(target);  % Ensure binary
%target = imresize(target, [1024, 1280], 'bilinear');


% % Step 0: Embed original target into a 4x4 grid and mirror it
% [orig_rows, orig_cols] = size(target);
% grid_rows = orig_rows * 2;
% grid_cols = orig_cols * 2;
% 
% % Create blank 4x4 canvas
% grid_target = zeros(grid_rows, grid_cols);
% 
% % Place original in top-left
% grid_target(1:orig_rows, 1:orig_cols) = target;
% % 
% % %Create 180° rotated version
% % rotated_target = rot90(target, 2);
% % 
% % %Place rotated version in bottom-right
% % grid_target(end-orig_rows+1:end, end-orig_cols+1:end) = rotated_target;
% 
% %Use this as the new target
% target = grid_target;





[num_y, num_x] = size(target);  % rows = y, cols = x

PixelEquations = cell(num_y, num_x);
curve_or_line = repmat(' ', num_y, num_x);
StoredPoints = nan(num_y, num_x, 6);  % [x0 y0 x1 y1 x2 y2]

% Step 3: Detect Connected Components
CC = bwconncomp(target, 4);
numObjects = CC.NumObjects;

% Step 4: Process Connected Components %% this code gives you the working version given to Wilkinson
for objIdx = 1:numObjects
    pixelIdxList = CC.PixelIdxList{objIdx};
    [y_coords, x_coords] = ind2sub(size(target), pixelIdxList);

    if length(x_coords) > 1
        [orderedCoords, ~] = ordered([x_coords, y_coords]);
    else
        orderedCoords = [x_coords, y_coords];
    end

    if size(orderedCoords, 1) >= 3
        A = orderedCoords(1, :);
        B = orderedCoords(round(end/2), :);
        C = orderedCoords(end, :);

        dAB = norm(A - B);
        dBC = norm(B - C);
        dAC = norm(A - C);

        if (dAB + dBC) < Alpha * dAC
            label = 'l';
            equation = sprintf('y = %.3fx + %.3f', (C(2)-A(2))/(C(1)-A(1)), A(2) - (C(2)-A(2))/(C(1)-A(1))*A(1));
            storedData = [A(1), A(2), C(1), C(2), NaN, NaN];
        else
            label = 'c';
            coeffs = polyfit([A(1); B(1); C(1)], [A(2); B(2); C(2)], 2);
            equation = sprintf('y = %.3fx^2 + %.3fx + %.3f', coeffs);
            storedData = [A(1), A(2), B(1), B(2), C(1), C(2)];
        end

        x = A(1); y = A(2);
        if x > 0 && y > 0 && x <= num_x && y <= num_y
            curve_or_line(y, x) = label;
            PixelEquations{y, x} = equation;
            StoredPoints(y, x, :) = storedData;
        end
    end
end

% %Step 4: merging of consecutive lines if they're too short:
% % Adjustable merge threshold
% merge_thresh = 10;  % Change as needed
% 
% % Store merged segments for processing later
% mergedSegments = {};  
% mergedLabels = {};
% mergedEquations = {};
% 
% for objIdx = 1:numObjects
%     pixelIdxList = CC.PixelIdxList{objIdx};
%     [y_coords, x_coords] = ind2sub(size(target), pixelIdxList);
% 
%     if length(x_coords) > 1
%         [orderedCoords, ~] = ordered([x_coords, y_coords]);
%     else
%         orderedCoords = [x_coords, y_coords];
%     end
% 
%     % Combine segments that are very close together
%     while objIdx < numObjects
%         nextIdx = objIdx + 1;
%         nextPixelIdxList = CC.PixelIdxList{nextIdx};
%         [ny, nx] = ind2sub(size(target), nextPixelIdxList);
%         if length(nx) > 1
%             [nextOrdered, ~] = ordered([nx, ny]);
%         else
%             nextOrdered = [nx, ny];
%         end
% 
%         if isempty(nextOrdered) || isempty(orderedCoords)
%             break;
%         end
% 
%         % Check if start of next segment is close to end of current
%         dist = norm(orderedCoords(end,:) - nextOrdered(1,:));
%         if dist < merge_thresh
%             orderedCoords = [orderedCoords; nextOrdered];  % Merge
%             objIdx = objIdx + 1;  % Move to next segment
%         else
%             break;
%         end
%     end
% 
%     % Classification
%     if size(orderedCoords, 1) >= 3
%         A = orderedCoords(1, :);
%         B = orderedCoords(round(end/2), :);
%         C = orderedCoords(end, :);
% 
%         dAB = norm(A - B);
%         dBC = norm(B - C);
%         dAC = norm(A - C);
% 
%         if (dAB + dBC) < Alpha * dAC
%             label = 'l';
%             equation = sprintf('y = %.3fx + %.3f', (C(2)-A(2))/(C(1)-A(1)), A(2) - (C(2)-A(2))/(C(1)-A(1))*A(1));
%             storedData = [A(1), A(2), C(1), C(2), NaN, NaN];
%         else
%             label = 'c';
%             coeffs = polyfit([A(1); B(1); C(1)], [A(2); B(2); C(2)], 2);
%             equation = sprintf('y = %.3fx^2 + %.3fx + %.3f', coeffs);
%             storedData = [A(1), A(2), B(1), B(2), C(1), C(2)];
%         end
% 
%         mergedSegments{end+1} = storedData;
%         mergedLabels{end+1} = label;
%         mergedEquations{end+1} = equation;
%     end
% end
% 
% % Store results into full-resolution arrays
% for k = 1:length(mergedSegments)
%     data = mergedSegments{k};
%     x = round(data(1)); y = round(data(2));
%     if x > 0 && y > 0 && x <= num_x && y <= num_y
%         curve_or_line(y, x) = mergedLabels{k};
%         PixelEquations{y, x} = mergedEquations{k};
%         StoredPoints(y, x, :) = data;
%     end
% end
% 






% Step 5: Visualize classifications
figure; imshow(target); hold on;
for x = 1:num_x
    for y = 1:num_y
        if curve_or_line(y,x) == 'l'
            plot(x, y, 'g.', 'MarkerSize', 5);
        elseif curve_or_line(y,x) == 'c'
            plot(x, y, 'r.', 'MarkerSize', 5);
        end
    end
end
hold off; title('Classified Lines (Green) and Curves (Red)');

% Step 6: Hologram Grid
scaled_x = round(num_x * scalefac);
scaled_y = round(num_y * scalefac);
[U, V] = meshgrid(linspace(-0.5, 0.5, scaled_x), linspace(-0.5, 0.5, scaled_y));% was -1,1 originally
%[U, V] = meshgrid(linspace(0, 1, scaled_x), linspace(0, 1, scaled_y));
%[U, V] = meshgrid(1:scaled_x, 1:scaled_y);


Holo = zeros(scaled_y, scaled_x);
ContributionCount = zeros(scaled_y, scaled_x);
processed_segments = zeros(scaled_y, scaled_x);

% Step 7: Apply Hologram Equations
for x = 1:num_x
    for y = 1:num_y
        if y <= scaled_y && x <= scaled_x && ~isnan(StoredPoints(y,x,1)) && processed_segments(y,x) == 0
        if ~isnan(StoredPoints(y,x,1)) && processed_segments(y,x) == 0
            pts = squeeze(StoredPoints(y,x,:));

            if curve_or_line(y,x) == 'l' && ~isnan(pts(3))
                H = LineHolo(pts(1), pts(2), pts(3), pts(4), U, V);
                H(isnan(H)) = 0;
                Holo = Holo + H;
                ContributionCount(H ~= 0) = ContributionCount(H ~= 0) + 1;
                processed_segments(round(pts(2)), round(pts(1))) = 1;

            elseif curve_or_line(y,x) == 'c' && ~isnan(pts(5))
                H = Arc2(pts(1), pts(2), pts(3), pts(4), pts(5), pts(6), U, V, 10);
                H(isnan(H)) = 0;
                Holo = Holo + H;
                ContributionCount(H ~= 0) = ContributionCount(H ~= 0) + 1;
                processed_segments(round(pts(2)), round(pts(1))) = 1;
            end
        end
        end
    end
end
%%
Holo(ContributionCount > 0) = Holo(ContributionCount > 0) ./ ContributionCount(ContributionCount > 0);
save(fullfile(outputFolder, ["Var_Holo_.mat"]), "Holo")

%%
% %% Step 8: FFT and Display
% Replay = ifft2(Holo);
% figure; imshow(mat2gray(abs(Holo))); title('Final Hologram');
% figure; imshow(mat2gray(abs(Replay))); title('Replay Field');
% 
% % Optional padding before FFT for cleaner edges
% pad_factor = 1;
% pad_y = round(size(Holo, 1) * pad_factor);
% pad_x = round(size(Holo, 2) * pad_factor);
% Holo_padded = padarray(Holo, [pad_y, pad_x], 0, 'both');
% Replay2 = fft2(Holo_padded);
% Replay2 = Replay2(pad_y+1:end-pad_y, pad_x+1:end-pad_x);

% %% --- Binary Phase Modulation ---
% % Create binary phase mask (0 or π)
% binary_mask = mod(angle(Holo), pi) > pi/2;  % logical
% binary_mask_uint8 = uint8(binary_mask);    % 0 or 1
% 
% % Simulate replay from binary phase mask
% Replay_binary_phase = ifft2(exp(1j * pi * double(binary_mask_uint8)));
% 
% % Display
% figure; imshow(binary_mask_uint8, []); title('Binary Phase Hologram (0 or 1)');
% figure; imshow(mat2gray(abs(Replay_binary_phase))); title('Replay from Binary Phase');
% 
% %% --- Save Everything ---
% 
% %% Convert continuous Holo phase into strict binary phase hologram
% 
% % Step 1: Get binary 0 or 1 mask from Holo's phase
% % Note: this is a simple binarization — not optimal, just consistent
% H_phase = angle(Holo);
% H_binary = uint8(mod(H_phase, pi) > pi/2);  % 0 or 1
% 
% % Step 2: Resize for projector (you said this was essential)
% H_binary_scaled = imresize(double(H_binary), [1024, 1280], 'bilinear');
% 
% % Step 3: Convert to 8-bit grayscale where 0 = 0 phase, 255 = π phase
% H_binary_uint8 = uint8(round(H_binary_scaled * 255));
% 
% % Step 4: Flip vertically if needed
% H_binary_uint8 = flipud(H_binary_uint8);
% 
% % Step 5: Save this final hologram
% imwrite(H_binary_uint8, fullfile(outputFolder, "Binary_Phase_Hologram_Sized.bmp"), 'bmp');
% 
% % (Optional) Also simulate replay
% Replay_binary = ifft2(exp(1j * pi * double(H_binary)));
% figure; imshow(mat2gray(abs(Replay_binary))); title('Replay from Naive Binary Phase');
% imwrite(mat2gray(abs(Replay_binary)), fullfile(outputFolder, "Binary_Phase_Replay_Naive.bmp"));


% %% Alt Save Everything Code:
% % Display the image
% figure; imshow(Holo); title('Hologram');
% 
% % Save the image
% imwrite(Holo, fullfile(outputFolder, shape + "_holo.png"));
% imwrite(Holo, fullfile(outputFolder, shape + "_holo.bmp"));
% 
% % Compute and save the replay field
% replay_field = ifft2(Holo);
% replay_field_normalized = mat2gray(abs(replay_field));
% 
% figure; imshow(replay_field_normalized); title('Reconstructed Image');
% imwrite(replay_field_normalized, fullfile(outputFolder, shape + "_replay_field.png"));

%% Step 2: Overlay "JONES APERTURE"
textSize = 50;
[grid_rows, grid_cols]=size(H);
position = [round(grid_cols * 0.3), round(grid_rows * 0.5)];

% Create a blank mask to overlay text
textMask = ones(grid_rows, grid_cols); 
textMask = insertText(textMask, position, 'JONES APERTURE', 'FontSize', textSize, 'BoxColor', 'black', 'TextColor', 'white');
textMask = im2bw(rgb2gray(textMask)); % Convert to binary mask

% Convert to phase modulation (0 or pi)
phaseModulation = pi * textMask; 

% % Apply the phase modulation to the hologram
% Holo = Holo .* exp(1j * phaseModulation);

% %% Step 3: Apply the Jones Matrix
% tau = pi;
% thetajones = pi/2;
% J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
%     -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];
% 
% % Apply the Jones Matrix
% Ex_mod = J(1,1) * Holo + J(1,2) * Holo;
% Ey_mod = J(2,1) * Holo + J(2,2) * Holo;
% 
% % Compute the final hologram before FFT
% Holo_mod = Ex_mod + Ey_mod;
% 
% % Convert to binary phase (0 or π)
% Holo_mod = mod(angle(Holo_mod), pi) > (pi/2);
% Holo_mod = Holo_mod * pi;
% 
% % Display Results
% figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or π)');
% figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image from Binary Phase');
% 
% % Save the Binary Hologram
% imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "_Binary_Phase_Hologram.bmp"));
% imwrite(mat2gray(abs(Replay)), fullfile(outputFolder,  "_Binary_Replay.bmp"));
% 
% % Step 7: Also save a correctly sized hologram
% Holo_mod_sized = imresize(Holo_mod, [1024, 1280], 'bilinear');
% imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "_Holo_Sized.bmp"));

%% Just the Basic Saving Code: 
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

% % Step 3: Convert to Binary Phase (-pi or pi),,, used to be (0 or pi)
Holo_mod = mod(angle(Holo_mod), pi) > (pi/2);  % Binary mask: 0 or 1
%Holo_mod = Holo_mod * 2 - 1;                  % Map: 0 → -1, 1 → +1
Holo_mod = Holo_mod * pi;                    % Final: -π or +π

% Step 4: Compute the Replay Field
Replay = (fft2(Holo)); % FFT of binary phase modulated hologram

Holo_image = uint8(Holo_mod / pi * 255);  % Converts 0 or π → 0 or 255
imwrite(Holo_image, fullfile(outputFolder, "Binary_Phase_Hologram_8bit.bmp"));


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

% === Step 4b: Create Amplitude-Modulated Hologram (for comparison) ===

% Use same mask as binary phase (or define your own threshold logic if needed)
amplitudeMask = mod(angle(Holo_mod), pi) > (pi/2);  % Same thresholding

% % Convert to -1 / +1 amplitude pattern
% ampHolo = amplitudeMask * 2 - 1;  % 0 → -1, 1 → +1

% % Rescale to 0–255 for display (−1 → 0, +1 → 255)
% ampImage = uint8((ampHolo + 1) / 2 * 255);

% % Save amplitude modulated hologram
% imwrite(ampImage, fullfile(outputFolder, "AmplitudeModulatedHologram.bmp"));
% imwrite(ampImage, fullfile(outputFolder, "AmplitudeModulatedHologram.png"));
