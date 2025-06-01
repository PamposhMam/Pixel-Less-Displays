clc; clear; close all;
% Define the root folder containing the images
rootFolder = "C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\classic-cgh-algorithms\IQM_Tests\Test3_comparisonqs";

% Get a list of all image files recursively
% Get a list of all image files recursively (.jpg, .bmp, .png)
imageFiles = [
    dir(fullfile(rootFolder, '**', '*.jpg'));
    dir(fullfile(rootFolder, '**', '*.bmp'));
    dir(fullfile(rootFolder, '**', '*.png'))
];
  % Change extension if needed

% Loop over each image file
for k = 1:numel(imageFiles)
    %% Load Grayscale Image
    imagePath = fullfile(imageFiles(k).folder, imageFiles(k).name);
    [~, imageName, ~] = fileparts(imagePath); % Extract filename without extension
    img = imread(imagePath);
    [~, imageName, ~] = fileparts(imagePath); % Extract filename
    if size(img,3) == 3
        img = rgb2gray(img); % Convert to grayscale if needed
    end
    img = double(img); % Convert to double for calculations
    
    %% set-up the output
    % Extract folder path from imagePath to save output in the same location
    [hostFolder, imageName, ~] = fileparts(imagePath); % Extract folder path and filename
    sanitizedImageName = regexprep(imageName, '[\/:*?"<>|]', '_'); % Sanitize filename
    
    % Skip processing if 'BinaryHolo' is found in the filename- don't
    % require the evaluation of apertures lmao
    if contains(imageName, 'BinaryHolo', 'IgnoreCase', true)
        fprintf('Skipping analysis for: %s (contains "BinaryHolo")\n', imageName);
        continue;
    end

    %% Image Size and Quadrant Extraction
    [rows, cols] = size(img);
    halfRows = floor(rows / 2);
    halfCols = floor(cols / 2);
    
    % Extract four quadrants
    quadrants = {
        img(1:halfRows, 1:halfCols),       % Top-Left
        img(1:halfRows, halfCols+1:end),  % Top-Right
        img(halfRows+1:end, 1:halfCols),  % Bottom-Left
        img(halfRows+1:end, halfCols+1:end) % Bottom-Right
    };
    
    % Compute total intensity for each quadrant
    totalIntensities = cellfun(@(q) sum(q(:)), quadrants);
    
    % Select quadrant with the highest intensity
    [~, maxIdx] = max(totalIntensities);
     quad_img = quadrants{maxIdx};
    
    % Identify quadrant location (for debugging/logging purposes)
    locations = {'Top-Left', 'Top-Right', 'Bottom-Left', 'Bottom-Right'};
    selectedLocation = locations{maxIdx}; %prev locations{maxIdx}
    fprintf('Selected Quadrant: %s (Highest Intensity)\n', selectedLocation);
    
    %Angle for the CoM stuff later on:
    angle=atan2(rows,cols);
        
    
    %% Metrics Calculation
    % Brightness (Mean Intensity)
    brightness = mean(quad_img(:));
    fullbrightness=mean(img(:));
    
    % Light 'noise' Spread (Standard Deviation)
    noiseSpread = std(quad_img(:));
    noiseSpreadFull=std(img(:));
    
    % Coverage (Percentage of Significant Pixels)
    threshold = max(quad_img(:)) * 0.05; % 5% of max intensity
    coverage = sum(quad_img(:) > threshold) / numel(quad_img) * 100;
    fullthresh= max(img(:))*0.05;
    fullcoverage= sum(img(:) > threshold) / numel(img) * 100;
    
    % Center of Mass (COM)
    [quadRows, quadCols] = size(quad_img);
    [X, Y] = meshgrid(1:quadCols, 1:quadRows);
    totalIntensity = sum(quad_img(:));
    COM_x = sum(X(:) .* quad_img(:)) / totalIntensity;
    COM_y = sum(Y(:) .* quad_img(:)) / totalIntensity;

    %Proposed Matlab fix:
    
        % Flip COM coordinates based on quadrant
        % switch maxIdx
        %     case 1
        %     COM_x = COM_x; % No change needed
        %     COM_y = COM_y;
        %     case 2
            COM_x = cols - quadCols + COM_x;
            COM_y = COM_y;
        %     case 3
        %     COM_x = COM_x;
        %     COM_y = rows - quadRows + COM_y;
        %     case 4
        %     COM_x = cols - quadCols + COM_x;
        %     COM_y = rows - quadRows + COM_y;
        % end

    COM = [COM_x, COM_y];
    
    %Centre of Image
    imCentreX= cols/2;
    imCentreY= rows/2;
    
    %Let's move everything globally now!

    % switch maxIdx
    %     case 1 % Top-Left
    %         Y_dot= imCentreY+COM_y;
    %         X_dot= imCentreX-COM_x;
    %     case 2 % Top-Right
    %         Y_dot= imCentreY+COM_y;
    %         X_dot= imCentreX+COM_x;
    %     case 3 % Bottom-Left
    %          Y_dot= imCentreY+COM_y;
    %         X_dot= imCentreX+COM_x;
    %     case 4 % Bottom-Right
    %          Y_dot= imCentreY+COM_y;
    %         X_dot= imCentreX+COM_x;
    %     otherwise
    %         error('Invalid quadrant index.');
    % end



    % Distance from COM to Quadrant Center as a ratio along the diagonal of the
    % quadrant
    quad_center = [quadCols / 2, quadRows / 2];
    COM_distance = sqrt(sum((COM - quad_center).^2));
    quad_diag_len= sqrt(quadRows^2+quadCols^2);
    quad_diag_angle = atan2(quadRows, quadCols);
    dist_ratio= (COM_distance*cos(quad_diag_angle))/quad_diag_len; %so basically what i'm doing is aligning it to the diag then normalizing it
    dist_ang= atan2(COM_y,COM_x);
    disp_from_maindiag= quad_diag_angle-dist_ang;
    
    %Figure features
    GCF=gcd(cols,rows);     %lol they use 'divisor' it's FACTOR, wenches!
    colsar= cols/GCF;
    rowsar=rows/GCF; %the rows for the aspect ratio.
    aspectrat=sprintf([num2str(colsar) ':' num2str(rowsar)]);

    %% Debugging COM with Red Dot and Diagonal Line
%% Debugging COM with Red Dot and Diagonal Line- just the quadrant
if contains(imageName, 'grid', 'IgnoreCase', true)
        % Add a red dot for COM and a line for the quadrant diagonal
        quadRGB = uint8(repmat(mat2gray(quad_img) * 255, 1, 1, 3)); % Convert quadrant to grayscale and make RGB
        
        % Define the COM coordinates relative to the quadrant
        dotX = round(COM_x);
        dotY = round(COM_y);
        
        % Ensure the red dot stays within the bounds of the quadrant
        dotX = max(min(dotX, quadCols), 1);
        dotY = max(min(dotY, quadRows), 1);
        
        % Overlay the red dot (draw a circle at the COM)
        dotRadius = 10; % Size of the red dot
        [Xgrid, Ygrid] = meshgrid(1:quadCols, 1:quadRows);
        circleMask = (Xgrid - dotX).^2 + (Ygrid - dotY).^2 <= dotRadius^2;
        
        % Apply the red dot to each color channel
        for channel = 1:3
            quadChannel = quadRGB(:, :, channel); % Extract the specific channel
            quadChannel(circleMask) = 255*(channel); % Apply red color to the mask area
            quadRGB(:, :, channel) = quadChannel; % Update the channel in the image
        end
        
        % Add a diagonal line to the quadrant
        lineColor = [0, 255, 0]; % Green
        diagonalMask = abs(Ygrid - Xgrid) <= 1; % Line from top-left to bottom-right
        
        % Apply the diagonal line to each color channel
        for channel = 1:3
            quadChannel = quadRGB(:, :, channel); % Extract the specific channel
            quadChannel(diagonalMask) = lineColor(channel); % Apply line color to the mask area
            quadRGB(:, :, channel) = quadChannel; % Update the channel in the image
        end
        
        % Save the annotated quadrant image
        outputAnnotatedPath = fullfile(hostFolder, sprintf('%s_annotated.png', sanitizedImageName));
        imwrite(quadRGB, outputAnnotatedPath);
        fprintf('Annotated quadrant image saved as: %s\n', outputAnnotatedPath);
end

% % Add a red dot and diagonal line to the entire image if 'grid' is in the filename
% if contains(imageName, 'grid', 'IgnoreCase', true)
%     % Convert the grayscale image to RGB for visualization
%     imgRGB = uint8(cat(3, img, img, img)); % Convert to RGB
% 
%     % Define the coordinates for the red dot (center of mass in the quadrant)
%     dotRadius = 10; % Size of the red dot
%     dotXGlobal = round(COM_x + (maxIdx == 2 || maxIdx == 4) * halfCols); % Adjust X based on quadrant
%     dotYGlobal = round(COM_y + (maxIdx == 3 || maxIdx == 4) * halfRows); % Adjust Y based on quadrant
%     [Xgrid, Ygrid] = meshgrid(1:cols, 1:rows);
%     circleMask = (Xgrid - dotXGlobal).^2 + (Ygrid - dotYGlobal).^2 <= dotRadius^2;
% 
%     % Apply the red dot to each color channel
%     for channel = 1:3
%         imgChannel = imgRGB(:, :, channel); % Extract the specific channel
%         imgChannel(circleMask) = 255 * (channel == 1); % Red channel gets 255, others get 0
%         imgRGB(:, :, channel) = imgChannel; % Update the channel in the image
%     end

%     % Add a diagonal line for the selected quadrant
%     diagStart = [1, 1]; % Top-left of the quadrant
%     diagEnd = [halfCols, halfRows]; % Bottom-right of the quadrant
%     if maxIdx == 2
%         diagStart = [halfCols + 1, 1]; diagEnd = [cols, halfRows]; % Top-right quadrant
%     elseif maxIdx == 3
%         diagStart = [1, halfRows + 1]; diagEnd = [halfCols, rows]; % Bottom-left quadrant
%     elseif maxIdx == 4
%         diagStart = [halfCols + 1, halfRows + 1]; diagEnd = [cols, rows]; % Bottom-right quadrant
%     end
% 
%     % Draw the diagonal line on the image
%     lineMask = drawLineMask(rows, cols, diagStart, diagEnd);
%     for channel = 1:3
%         imgChannel = imgRGB(:, :, channel); % Extract the specific channel
%         imgChannel(lineMask) = 255 * (channel == 2); % Green line
%         imgRGB(:, :, channel) = imgChannel; % Update the channel in the image
%     end
% 
%     % Save the modified image
%     outputFilePath = fullfile(hostFolder, sprintf('%s_fullwithRedDotAndDiagonal.jpg', sanitizedImageName));
%     imwrite(imgRGB, outputFilePath);
%     fprintf('Red dot and diagonal line added and saved as: %s\n', outputFilePath);
% end



    %% Display Results
    
    
metricsTable = table( ...
    {'Selected Quadrant'; 'Quadrant Intensity Spread (std)'; 'Quadrant Brightness (mean)'; 'Quadrant Coverage (%)'; 'COM Distance Ratio'; 'COM Angle (for main quadrant, wrt centre (radians))'; ...
    'COM Angle from the Centreline'; 'Total Pixels'; 'Full Intensity Spread'; 'Full Brightness'; 'Full Coverage'; 'Width (pixels)'; 'Height (pixels)'; 'Aspect Ratio'; 'COM_xcoord'; 'COM_ycoord'}, ...
    {selectedLocation; noiseSpread; brightness; coverage; dist_ratio; dist_ang; disp_from_maindiag; numel(img); noiseSpreadFull; fullbrightness; fullcoverage; cols; rows; aspectrat; COM_x; COM_y}, ...
    'VariableNames', {'Metric', 'Value'});

    
    disp('Image Analysis Results:');
    disp(metricsTable);
    
    % Create a new figure for displaying the table
    fig = figure('Visible', 'off');  % Create the figure without showing it
    
    % Create a uitable to display the metrics
    uit=uitable('Parent', fig, 'Data', table2cell(metricsTable), ...
            'ColumnName', metricsTable.Properties.VariableNames, ...
            'Units', 'Normalized', 'Position', [0, 0, 2, 1]);
    
   % % Increase the font size and adjust column width manually
    columnWidths = {350, 300};  % Example column width increase (in pixels)
    set(uit, 'ColumnWidth', columnWidths);
    
    
    %% Save Metrics Table as Image
    % Prepare data for uitable
    data = [metricsTable.Metric, cellfun(@num2str, metricsTable.Value, 'UniformOutput', false)];
    
    % Create a cell array for all metrics with proper individual rows for readability
    allMetrics = [ ...
        {"Selected Quadrant", selectedLocation}; 
        {"Quadrant Intensity Spread (std)", noiseSpread}; 
        {"Quadrant Brightness (mean)", brightness}; 
        {"Quadrant Coverage (%)", coverage}; 
        {"COM Distance Ratio", dist_ratio}; 
        {"COM Angle (for main quadrant, wrt centre) (radians)", dist_ang}; 
        {"COM Angle from the Centreline", disp_from_maindiag}; 
        {"Total Pixels", numel(quad_img)}; 
        {"Full Intensity Spread (std)", noiseSpreadFull}; 
        {"Full Brightness (mean)", fullbrightness}; 
        {"Full Coverage (%)", fullcoverage}
        {"Width (pixels)", cols}
        {"Height (pixels)", rows}
        {"Aspect Ratio", aspectrat}
        {"COM_xcoord", COM_x}
        {"COM_ycoord", COM_y}
    ];
    
    % Convert to table for display
    metricsTable = cell2table(allMetrics, 'VariableNames', {'Metric', 'Value'});
    
    %
    % Adjust figure size to accommodate table
    set(fig, 'Position', [100, 100, 800, 600]);  % Adjust dimensions if needed
    
    % Capture and save the table as an image
    frame = getframe(fig);
    
    % Convert frame data to uint8 if needed
    imData = frame.cdata;
    if ~isa(imData, 'uint8')
        imData = im2uint8(mat2gray(imData));  % Normalize and convert to uint8
    end
    
    
    % Sanitize filename to remove problematic characters
    %sanitizedImageName = regexprep(imageName, '[\/:*?"<>|]', '_');  % Replace invalid characters with '_'
    
    % Save using sanitized filename
    outputFilePath = fullfile(hostFolder, sprintf('image_analysis_%s_metrics.png', sanitizedImageName));
    imwrite(imData, outputFilePath);
% Save the metrics table as a .txt file (human-readable format)
outputTextPath = fullfile(hostFolder, sprintf('image_analysis_%s_metrics.txt', sanitizedImageName));
fid = fopen(outputTextPath, 'w');
for i = 1:size(metricsTable, 1)
    fprintf(fid, '%s: %s\n', metricsTable.Metric{i}, convertToPrintable(metricsTable.Value{i}));
end
fclose(fid);


    
    close(fig);
end

%% Function holder
% Function to create a line mask
function mask = drawLineMask(height, width, startPt, endPt)
    mask = false(height, width);
    x1 = startPt(1); y1 = startPt(2);
    x2 = endPt(1); y2 = endPt(2);
    numPoints = max(abs(x2 - x1), abs(y2 - y1)) + 1;
    xCoords = round(linspace(x1, x2, numPoints));
    yCoords = round(linspace(y1, y2, numPoints));
    for i = 1:numPoints
        if xCoords(i) > 0 && xCoords(i) <= width && yCoords(i) > 0 && yCoords(i) <= height
            mask(yCoords(i), xCoords(i)) = true;
        end
    end
end


%%
function outStr = convertToPrintable(val)
    if isnumeric(val)
        outStr = num2str(val);
    elseif ischar(val)
        outStr = val;
    elseif isstring(val)
        outStr = char(val);
    elseif islogical(val)
        outStr = mat2str(val);
    else
        try
            outStr = evalc('disp(val)');  % fallback for weird stuff
        catch
            outStr = '<unprintable>';
        end
    end
end
