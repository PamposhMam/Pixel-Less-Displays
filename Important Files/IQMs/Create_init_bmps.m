% Set the main directory containing 'Images'
mainDir = 'C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\classic-cgh-algorithms\IQM_Tests\Images4';

% Get a list of all folders in the directory
folders = dir(mainDir);
folders = folders([folders.isdir]); % Only keep directories
folders = folders(~ismember({folders.name}, {'.', '..'})); % Exclude '.' and '..'

% Loop through each folder
for i = 1:length(folders)
    % if folders(i).name=='J' %Uncomment if you need a specific folder
        folderPath = fullfile(mainDir, folders(i).name);
        fprintf('Processing folder: %s\n', folderPath);
        
        % Get all bitmap files in the folder
        bmpFiles = dir(fullfile(folderPath, '*.bmp'));
        
        % Preprocess: Convert all .bmp files to .jpg
        for j = 1:length(bmpFiles)
            bmpPath = fullfile(folderPath, bmpFiles(j).name);
            [~, name, ~] = fileparts(bmpFiles(j).name);
            jpgPath = fullfile(folderPath, [name, '.jpg']);
            
            try
                img = imread(bmpPath);
                imwrite(img, jpgPath, 'jpg');
                fprintf('Converted %s to %s\n', bmpPath, jpgPath);
            catch ME
                fprintf('Error converting %s to .jpg\nError message: %s\n', bmpPath, ME.message);
            end
        end
        
        % Get all image files in the folder (.jpg and .png)
        jpgFiles = dir(fullfile(folderPath, '*.jpg'));
        pngFiles = dir(fullfile(folderPath, '*.png'));
        imageFiles = [jpgFiles; pngFiles];
        
        if isempty(imageFiles)
            fprintf('No images found in folder: %s\n', folderPath);
            continue; % Skip to the next folder
        end
        
        % Process images
        for j = 1:length(imageFiles)
            % Read the image
            imgPath = fullfile(folderPath, imageFiles(j).name);
            fprintf('Attempting to read image: %s\n', imgPath);
            
            try
                img = imread(imgPath);
            catch ME
                fprintf('Error reading image: %s\nError message: %s\n', imgPath, ME.message);
                continue; % Skip to the next image
            end
            
            % Convert to grayscale if not already
            if size(img, 3) == 3
                imgGray = rgb2gray(img);
            else
                imgGray = img;
            end
            
            % Create a black square grid 4x the size of the image
            [rows, cols] = size(imgGray);
            gridSize = [2 * rows, 2 * cols];
            blackGrid = zeros(gridSize, 'uint8');
            
            % Randomly select a quadrant for the original image
            quadrant = randi(4);
            fprintf('Selected quadrant: %d\n', quadrant);
            
            % Determine the placement for the original image
            switch quadrant
                case 1 % Top-left
                    rowStart = 1;
                    colStart = 1;
                case 2 % Top-right
                    rowStart = 1;
                    colStart = cols + 1;
                case 3 % Bottom-left
                    rowStart = rows + 1;
                    colStart = 1;
                case 4 % Bottom-right
                    rowStart = rows + 1;
                    colStart = cols + 1;
            end
            blackGrid(rowStart:rowStart+rows-1, colStart:colStart+cols-1) = imgGray;
            
            % Rotate the original image 180 degrees
            imgRotated = rot90(imgGray, 2);
            
            % Place the rotated image in the opposite quadrant
            oppositeRowStart = max(1, gridSize(1) - rowStart - rows + 1);
            oppositeColStart = max(1, gridSize(2) - colStart - cols + 1);
            blackGrid(oppositeRowStart:oppositeRowStart+rows-1, ...
                      oppositeColStart:oppositeColStart+cols-1) = imgRotated;
            
            % Save the black grid with the images
            [~, name, ~] = fileparts(imageFiles(j).name);
            outputGridPath = fullfile(folderPath, [name, '_grid.bmp']);
            try
                imwrite(blackGrid, outputGridPath);
                fprintf('Saved grid BMP: %s\n', outputGridPath);
            catch ME
                fprintf('Error saving grid BMP: %s\nError message: %s\n', outputGridPath, ME.message);
            end
        end
    % end %uncomment if you need a specific folder
end

disp('Processing complete.');
