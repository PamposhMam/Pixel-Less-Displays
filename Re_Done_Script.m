clear all
close all

% Step 1: Read and Resize the Image
filename = "Easy.bmp"; % Change this to the actual input file
[~, name, ~] = fileparts(filename); % Extract just the name without extension
target = imread(filename);
Alpha = 1.5;    
scalefac=1;
name= name + num2str(Alpha) + "_scalefac_" + num2str(scalefac);
% Step 2: Create a Folder with the Same Name
outputFolder = fullfile(pwd, name); % Create folder path
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % Make directory if it doesn't exist
end


figure;
imshow(target);
title('Edge-Defined Image');

% Step 2: Declare Variables
 
[num_rows, num_cols] = size(target);

% Use a cell array to store equations for each pixel
PixelEquations = cell(num_rows, num_cols);  
curve_or_line = repmat(' ', num_rows, num_cols);  

% New Storage: Stores Start, Middle (if needed), and End Points
StoredPoints = nan(num_rows, num_cols, 6);  % Each pixel stores [x0, y0, x1, y1, x2, y2]

% Step 3: Detect Connected Components
CC = bwconncomp(target, 4);
numObjects = CC.NumObjects;

% Step 4: Process Each Detected Edge Component
for objIdx = 1:numObjects
    pixelIdxList = CC.PixelIdxList{objIdx};  
    [rows, cols] = ind2sub(size(target), pixelIdxList);  

    % Only sort if at least two points exist
    if length(rows) > 1  
        orderedCoords = ordered([rows, cols]);  
    else
        orderedCoords = [rows, cols];  % Keep as is for single-point objects
    end

    % Classify as Line or Curve (Only if 3+ points exist)
    if size(orderedCoords, 1) >= 3  
        A = orderedCoords(1, :);   % First point in segment
        B = orderedCoords(round(end/2), :); % Midpoint (only used for curves)
        C = orderedCoords(end, :); % Last point in segment

        dist_AB = norm(A - B);
        dist_BC = norm(B - C);
        dist_AC = norm(A - C);

        if (dist_AB + dist_BC) < Alpha * dist_AC
            label = 'l';  % Line

            % Compute linear equation: y = mx + b
            m = (C(2) - A(2)) / (C(1) - A(1));  
            b = A(2) - m * A(1);  
            equation = sprintf('y = %.3fx + %.3f', m, b);

            % Store only start and end points for lines (one-time storage)
            storedData = [A(1), A(2), C(1), C(2), NaN, NaN];  

        else
            label = 'c';  % Curve

            % Compute quadratic equation: y = ax² + bx + c
            X = [A(1); B(1); C(1)];
            Y = [A(2); B(2); C(2)];
            coeffs = polyfit(X, Y, 2);  
            equation = sprintf('y = %.3fx^2 + %.3fx + %.3f', coeffs(1), coeffs(2), coeffs(3));

            % Store start, middle, and end points for curves
            storedData = [A(1), A(2), B(1), B(2), C(1), C(2)];
        end

        % Store classification, equation, and points **only once per segment**
        x = A(1);
        y = A(2);
        curve_or_line(x, y) = label;
        PixelEquations{x, y} = equation;
        StoredPoints(x, y, :) = storedData;
    end
end

%% The display case
% Step 5: Display Results
figure;
imshow(target);
hold on;
for i = 1:num_rows
    for j = 1:num_cols
        if curve_or_line(i,j) == 'l'
            plot(j, i, 'g.', 'MarkerSize', 5);  
        elseif curve_or_line(i,j) == 'c'
            plot(j, i, 'r.', 'MarkerSize', 5);  
        end
    end
end
hold off;
title('Classified Lines (Green) and Curves (Red)');

% Step 6: Access Example Equations (For Debugging)
fprintf('Equation at (100,200): %s\n', PixelEquations{100,200});

%% 🔥 **Step 5.5: Scale Up Before Applying the Hologram**
 scale_factor = scalefac;  % Upscaling factor
% 
% % Increase Grid Resolution
scaled_rows = round(num_rows * scale_factor);
scaled_cols = round(num_cols * scale_factor);
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

                    Holo_Func = LineHolo(x0, y0, x1, y1, U, V);
                    Holo_Func(isnan(Holo_Func)) = 0;
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
                    Holo_Func = ArcHolo(x0, y0, x1, y1, x2, y2, U, V, 10);
                  Holo_Func(isnan(Holo_Func)) = 0; %lazy fix but idgaf

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
Replay = fft2(Holo'); % FFT of binary phase modulated hologram



%%
% % Debugging: Display final Holo before FFT
% figure; imshow(mat2gray(abs(Holo))); title('Final Hologram Before FFT');
% 
% % Normalize and apply FFT for reconstruction
% Replay = (ifft2(Holo));
% 
% % % Debugging: Display replay field before saving
% % figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image Before Saving');
% 
% % Save and display results
% imwrite(mat2gray(abs(Holo)), 'Final_Hologram.bmp');
% imwrite(mat2gray(abs(Replay)), 'Final_Replay.bmp');
% 
% %Hologram Binarization Modulation
% %Jones Matrix
% tau= pi;
% thetajones= pi/2;
% 
% J= [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
%     -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];
% 
% 
% 
% % Display reconstructed hologram
% figure; imshow(mat2gray(abs(Holo))); title('Final Hologram');
% figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT)');


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