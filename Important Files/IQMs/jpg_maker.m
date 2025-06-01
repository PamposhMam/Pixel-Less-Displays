% Get the current folder
currentFolder = ("C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\LiveHologramsForProject\vector graphics images");

% Define the 'images' folder path
imagesFolder = fullfile(currentFolder);

% Recursively find all .bmp files
bmpFiles = dir(fullfile(imagesFolder, '**', '*.bmp'));

% Loop through each .bmp file
for k = 1:length(bmpFiles)
    % Get the full path of the current .bmp file
    bmpFilePath = fullfile(bmpFiles(k).folder, bmpFiles(k).name);
    
    % Read the .bmp image
    img = imread(bmpFilePath);
    
    % Create a .png file name by replacing the extension
    [~, name, ~] = fileparts(bmpFiles(k).name);
    pngFilePath = fullfile(bmpFiles(k).folder, [name, '.jpg']);
    
    % Write the image as .png
    imwrite(img, pngFilePath);
end

disp('All .bmp files have been converted to .jpg.');
