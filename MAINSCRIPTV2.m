clear all
close all

%% Step 1: Read and Resize the Image
filename = "Easy.bmp";
[~, name, ~] = fileparts(filename);
target = imread(filename);
Alpha = 2;  % Classification threshold
name = name + num2str(Alpha);

% Step 2: Create and Descend into Folder
outputFolder = fullfile(pwd, name);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Step 3: Pre-Processing
[num_rows, num_cols] = size(target);
PixelEquations = cell(num_rows, num_cols);
curve_or_line = repmat(' ', num_rows, num_cols);
StoredPoints = nan(num_rows, num_cols, 6);

CC = bwconncomp(target, 4);
numObjects = CC.NumObjects;

%% Step 4: Segment Detection with Per-Point Geometry Tracking
tolerance = Alpha; % how far a point can deviate to still be "on the line/curve"

for objIdx = 1:numObjects
    pixelIdxList = CC.PixelIdxList{objIdx};
    [rows, cols] = ind2sub(size(target), pixelIdxList);
    coords = ordered([rows, cols]);
    N = size(coords, 1);

    k = 1;
    while k <= N - 2
        A = coords(k, :);
        B = coords(k + 1, :);
        C = coords(k + 2, :);
        current_segment = [A; B; C];
        idx = k + 3;

        % Try to fit a line first
        while idx <= N
            D = coords(idx, :);
            temp = [current_segment; D];
            x = temp(:,1);
            y = temp(:,2);
            p_line = polyfit(x([1 end]), y([1 end]), 1);
            y_pred_line = polyval(p_line, x);
            error_line = mean(abs(y - y_pred_line));

            p_curve = polyfit(x, y, 2);
            y_pred_curve = polyval(p_curve, x);
            error_curve = mean(abs(y - y_pred_curve));

            if error_line < tolerance
                current_segment = [current_segment; D];
                label = 'l';
                equation = sprintf('y = %.3fx + %.3f', p_line(1), p_line(2));
            elseif error_curve < tolerance
                current_segment = [current_segment; D];
                label = 'c';
                equation = sprintf('y = %.3fx^2 + %.3fx + %.3f', p_curve(1), p_curve(2), p_curve(3));
            else
                break
            end
            idx = idx + 1;
        end

        % Save segment
        A = current_segment(1, :);
        C = current_segment(end, :);
        B = current_segment(round(end/2), :);
        if label == 'l'
            storedData = [A(1), A(2), C(1), C(2), NaN, NaN];
        else
            storedData = [A(1), A(2), B(1), B(2), C(1), C(2)];
        end

        for pt = 1:size(current_segment,1)
            x = current_segment(pt, 1);
            y = current_segment(pt, 2);
            curve_or_line(x, y) = label;
            PixelEquations{x, y} = equation;
            StoredPoints(x, y, :) = storedData;
        end

        k = idx; % Move to the next new segment
    end
end

%% Step 5: Visualization
figure; imshow(target); hold on;
for i = 1:num_rows
    for j = 1:num_cols
        if curve_or_line(i,j) == 'l'
            plot(j, i, 'g.', 'MarkerSize', 5);  
        elseif curve_or_line(i,j) == 'c'
            plot(j, i, 'r.', 'MarkerSize', 5);  
        end
    end
end
title('Segmented: Lines (Green), Curves (Red)');
hold off;


%% 🔥 **Step 5.5: Scale Up Before Applying the Hologram**
 scale_factor = 1;  % Upscaling factor
% 
% % Increase Grid Resolution
scaled_rows = num_rows * scale_factor;
scaled_cols = num_cols * scale_factor;
% 
% % Scale stored points (x0, y0, etc.)
% for i = 1:num_rows
%     for j = 1:num_cols
%         if ~isnan(StoredPoints(i, j, 1))
%             StoredPoints(i, j, :) = StoredPoints(i, j, :) * scale_factor;
%         end
%     end
% end

% Define High-Resolution Grid
U = linspace(-1, 1, scaled_cols);
V = linspace(-1, 1, scaled_rows);
[U, V] = meshgrid(U, V);

% Initialize High-Resolution Hologram
Holo = zeros(scaled_rows, scaled_cols);
ContributionCount = zeros(scaled_rows, scaled_cols);
processed_segments = zeros(scaled_rows, scaled_cols);

% **🔥 Step 6: Apply Hologram Equations at Higher Resolution**
equation_counter = 0;  % Initialize equation counter
display_flag = false;  % Flag to track when to display
nty= 400; %how many equations to skip before displaying an intermediate hologram

for i = 1:num_rows
    for j = 1:num_cols
        if ~isnan(StoredPoints(i, j, 1)) && processed_segments(i, j) == 0
            points = squeeze(StoredPoints(i, j, :)); %lmao, learned a new matlab function i guess

            if curve_or_line(i, j) == 'l' && ~isnan(points(3)) %if it's a line
                x0 = points(1); y0 = points(2); %x0,y0 assigned
                x1 = points(3); y1 = points(4); %x1,y1 assigned

                % if processed_segments(x0, y0) == 0
                %     disp(['Processing LineHolo (scaled) with x0 = ', num2str(x0), ', y0 = ', num2str(y0), ...
                %           ', x1 = ', num2str(x1), ', y1 = ', num2str(y1)]);
                    Holo_Func = LineHolo(x0, y0, x1, y1, U, V);
                %     if any(isinf(Holo_Func(:)))
                %         warning('🚨 Infinity detected in Holo_Func at equation %d!', equation_counter);
                %         disp('Max absolute value before correction:');
                %         disp(max(abs(Holo_Func(:))));
                % 
                %         % Cap the values to the max finite number in the array
                %         max_finite_value = max(Holo_Func(~isinf(Holo_Func(:))));
                %         Holo_Func(isinf(Holo_Func)) = max_finite_value;  % Set to largest valid number
                %     end
                    Holo = Holo + Holo_Func;
                    ContributionCount(Holo_Func ~= 0) = ContributionCount(Holo_Func ~= 0) + 1;
                    processed_segments(x0, y0) = 1;
                    processed_segments(x1, y1) = 1;
                % end
                equation_counter = equation_counter + 1;  
                
            elseif curve_or_line(i, j) == 'c' && ~isnan(points(5))
                x0 = points(1); y0 = points(2);
                x1 = points(3); y1 = points(4);
                x2 = points(5); y2 = points(6);

            
                    Holo_Func = Arc2(x0, y0, x1, y1, x2, y2, U, V, 10);
            

                    
                    Holo = Holo + Holo_Func;
                    ContributionCount(Holo_Func ~= 0) = ContributionCount(Holo_Func ~= 0) + 1;
                    processed_segments(x0, y0) = 1;
                    processed_segments(x1, y1) = 1;
                    processed_segments(x2, y2) = 1;
                % end
                equation_counter = equation_counter + 1;  
            end
        end

        % **🔥 Set Flag for Display Every nty Equations**
        if mod(equation_counter, nty) == 0 && equation_counter > 0  
            display_flag = true;  
        end
    end

    % **🔥 Display at the End of Each Row (Prevents Multiple Figures Per Iteration)**
    if display_flag  
        disp(['Displaying intermediate hologram at equation count: ', num2str(equation_counter)]);

        % Normalize Hologram before displaying
        nonzero_idx = ContributionCount > 0;
        Holo(nonzero_idx) = Holo(nonzero_idx) ./ ContributionCount(nonzero_idx);

        % Display Hologram
        figure;
        imshow(mat2gray(abs(Holo)));
        title(['Hologram at ', num2str(equation_counter), ' Equations']);

        % Compute and Display Replay Field
        Replay = fftshift(fft2(Holo));
        figure;
        imshow(mat2gray(abs(Replay)));
        title(['Replay Field at ', num2str(equation_counter), ' Equations']);

        display_flag = false;  % Reset flag so it only triggers every 200 equations
    end
end

% **🔥 Final Normalization**
nonzero_idx = ContributionCount > 0;
Holo(nonzero_idx) = Holo(nonzero_idx) ./ ContributionCount(nonzero_idx);
%% capture the replay field before all the bullshit that happens next:
Replay = fft2(Holo); % FFT of binary phase modulated hologram

%%

% Step 1: Generate a Phase Mask for "JONES APERTURE"
textSize = 50; % Adjust for visibility
position = [round(scaled_cols * 0.3), round(scaled_rows * 0.5)]; % Center it

% Create a blank mask to overlay text
textMask = ones(scaled_rows, scaled_cols); 
textMask = insertText(textMask, position, 'JONES APERTURE', 'FontSize', textSize, 'BoxColor', 'black', 'TextColor', 'white');
textMask = im2bw(rgb2gray(textMask)); % Convert to binary mask

% Convert binary mask into a phase modulation (0 or pi)
phaseModulation = pi * textMask; 

% Step 2: Apply the Phase Modulation to the Hologram
Holo = Holo .* exp(1j * phaseModulation); 

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


%%
% Step 5: Display Results
figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or π)');
figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT) from Binary Phase');

% Step 6: Save the Binary Hologram
imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "Binary_Phase_Hologram.bmp"));
imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "Binary_Replay.bmp"));

%Step 7: Also save a correct-sized hologram to the folder
Holo_mod_sized= imresize(Holo_mod, [1024, 1280], 'bilinear');
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.bmp"));