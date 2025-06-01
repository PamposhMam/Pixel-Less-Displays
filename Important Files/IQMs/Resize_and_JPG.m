% Get the current folder
currentFolder = pwd;

% Define the 'images' folder path
imagesFolder = fullfile(currentFolder, 'images');

% Recursively find all .bmp files
bmpFiles = dir(fullfile(imagesFolder, '**', '*.png'));

% Loop through each .bmp file
for k = 1:length(bmpFiles)
    % Get the full path of the current .bmp file
    bmpFilePath = fullfile(bmpFiles(k).folder, bmpFiles(k).name);
    
    % Read the .bmp image
    img = imread(bmpFilePath);
    
    % Get the size of the image and calculate the new size
    [rows, cols, ~] = size(img);
    newRows = round(rows / 2);
    newCols = round(cols / 2);
    
    % Resize the image to a quarter of its original size
    resizedImg = imresize(img, [newRows, newCols]);
    
    % Create a .jpg file name by replacing the extension
    [~, name, ~] = fileparts(bmpFiles(k).name);
    jpgFilePath = fullfile(bmpFiles(k).folder, [name, '.jpg']);
    
    % Write the resized image as a .jpg
    imwrite(resizedImg, jpgFilePath, 'Quality', 90); % Adjust quality as needed
end

disp('All .bmp files have been resized and converted to .jpg.');
