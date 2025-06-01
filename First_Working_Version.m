clear all
close all

% Step 1: Read and Resize the Image
target = imread("Easy.bmp");
target = imresize(target, [1280, 1080]);  % Resize correctly
target = im2bw(target);  % Convert to binary (Assuming edge image)

figure;
imshow(target);
title('Edge-Defined Image');

% Step 2: Declare Variables
Alpha = 15;     
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

% % Step 4: Process Each Detected Edge Component
% for objIdx = 1:numObjects
%     pixelIdxList = CC.PixelIdxList{objIdx};  
%     [rows, cols] = ind2sub(size(target), pixelIdxList);  
% 
%     % Only sort if at least two points exist
%     if length(rows) > 1  
%         orderedCoords = ordered([rows, cols]);  
%     else
%         orderedCoords = [rows, cols];  % Keep as is for single-point objects
%     end
% 
%     % Classify as Line or Curve (Only if 3+ points exist)
%     if size(orderedCoords, 1) >= 3  
%         A = orderedCoords(1, :);   % First point in segment
%         B = orderedCoords(round(end/2), :); % Midpoint (only used for curves)
%         C = orderedCoords(end, :); % Last point in segment
% 
%         dist_AB = norm(A - B);
%         dist_BC = norm(B - C);
%         dist_AC = norm(A - C);
% 
%         if (dist_AB + dist_BC) < Alpha * dist_AC
%             label = 'l';  % Line
% 
%             % Compute linear equation: y = mx + b
%             m = (C(2) - A(2)) / (C(1) - A(1));  
%             b = A(2) - m * A(1);  
%             equation = sprintf('y = %.3fx + %.3f', m, b);
% 
%             % Store only start and end points for lines
%             storedData = [A(1), A(2), C(1), C(2), NaN, NaN];  
% 
%         else
%             label = 'c';  % Curve
% 
%             % Compute quadratic equation: y = ax² + bx + c
%             X = [A(1); B(1); C(1)];
%             Y = [A(2); B(2); C(2)];
%             coeffs = polyfit(X, Y, 2);  
%             equation = sprintf('y = %.3fx^2 + %.3fx + %.3f', coeffs(1), coeffs(2), coeffs(3));
% 
%             % Store start, middle, and end points for curves
%             storedData = [A(1), A(2), B(1), B(2), C(1), C(2)];
%         end
% 
%         % 🔥 **Store classification, equation, AND points for EVERY pixel**
%         for k = 1:length(orderedCoords)
%             x = orderedCoords(k, 1);
%             y = orderedCoords(k, 2);
%             curve_or_line(x, y) = label;  
%             PixelEquations{x, y} = equation;
%             StoredPoints(x, y, :) = storedData;  % ✅ Assign full storedData
%         end
%     end
% end

% Step 5: **Debug: Count Classified Pixels**
disp(['Total Pixels Classified as Lines: ', num2str(sum(curve_or_line(:) == 'l'))]);
disp(['Total Pixels Classified as Curves: ', num2str(sum(curve_or_line(:) == 'c'))]);

% Step 5: Display Results
figure;
imshow(target);
hold on;
for i = 1:num_rows
    for j = 1:num_cols
        if curve_or_line(i,j) == 'l'
            plot(j, i, 'g.', 'MarkerSize', 5);  % Green for lines
        elseif curve_or_line(i,j) == 'c'
            plot(j, i, 'r.', 'MarkerSize', 5);  % Red for curves
        end
    end
end
hold off;
title('Classified Lines (Green) and Curves (Red)');

% Step 6: Access Example Equations (For Debugging)
fprintf('Equation at (100,200): %s\n', PixelEquations{100,200});
fprintf('Stored Points at (100,200): [%s]\n', num2str(squeeze(StoredPoints(100,200,:))'));


% Initialize hologram plane
Holo = zeros(num_rows, num_cols);
ContributionCount = zeros(num_rows, num_cols); % Track how many times a pixel is updated

% Define U, V axes
U = linspace(-1, 1, size(Holo, 2));
V = linspace(-1, 1, size(Holo, 1));
[U, V] = meshgrid(U, V);

% Keep track of processed segments to avoid duplicate processing
processed_segments = zeros(num_rows, num_cols);

% Loop over all pixels in curve_or_line
for i = 1:num_rows
    for j = 1:num_cols
        % If pixel has stored points (not NaN) and has not been processed
        if ~isnan(StoredPoints(i, j, 1)) && processed_segments(i, j) == 0
            % Extract stored coordinate points
            points = squeeze(StoredPoints(i, j, :));

            % Mark segment as processed
            processed_segments(i, j) = 1;

            % Line Processing (Only process start pixel of a segment)
            if curve_or_line(i, j) == 'l' && ~isnan(points(3))  
                x0 = points(1); y0 = points(2);
                x1 = points(3); y1 = points(4);

                if processed_segments(x0, y0) == 0
                    disp(['Processing LineHolo with x0 = ', num2str(x0), ', y0 = ', num2str(y0), ...
                          ', x1 = ', num2str(x1), ', y1 = ', num2str(y1)]);
                    Holo_Func = LineHolo(x0, y0, x1, y1, U, V);
                    Holo = Holo + Holo_Func;

                    % Increase contribution count
                    ContributionCount(Holo_Func ~= 0) = ContributionCount(Holo_Func ~= 0) + 1;

                    % Mark the entire segment as processed
                    processed_segments(x0, y0) = 1;
                    processed_segments(x1, y1) = 1;
                end

            % Curve Processing
            elseif curve_or_line(i, j) == 'c' && ~isnan(points(5))  
                x0 = points(1); y0 = points(2);
                x1 = points(3); y1 = points(4);
                x2 = points(5); y2 = points(6);

                if processed_segments(x0, y0) == 0
                    disp(['Processing ArcHolo with x0 = ', num2str(x0), ', y0 = ', num2str(y0), ...
                          ', x1 = ', num2str(x1), ', y1 = ', num2str(y1), ...
                          ', x2 = ', num2str(x2), ', y2 = ', num2str(y2)]);
                    Holo_Func = ArcHolo(x0, y0, x1, y1, x2, y2, U, V, 10);
                    Holo = Holo + Holo_Func;

                    % Increase contribution count
                    ContributionCount(Holo_Func ~= 0) = ContributionCount(Holo_Func ~= 0) + 1;

                    % Mark the entire segment as processed
                    processed_segments(x0, y0) = 1;
                    processed_segments(x1, y1) = 1;
                    processed_segments(x2, y2) = 1;
                end
            end
        end
    end
end

% **🔥 Normalize Hologram by Contribution Count**
nonzero_idx = ContributionCount > 0;
Holo(nonzero_idx) = Holo(nonzero_idx) ./ ContributionCount(nonzero_idx);


%%
% Debugging: Display final Holo before FFT
figure; imshow(mat2gray(abs(Holo))); title('Final Hologram Before FFT');

% Normalize and apply FFT for reconstruction
Replay = fftshift(fft2(Holo));

% Debugging: Display replay field before saving
figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image Before Saving');

% Save and display results
imwrite(mat2gray(abs(Holo)), 'Final_Hologram.bmp');
imwrite(mat2gray(abs(Replay)), 'Final_Replay.bmp');

% Display reconstructed hologram
figure; imshow(mat2gray(abs(Holo))); title('Final Hologram');
figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT)');