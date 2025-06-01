%% Montage Code with Labels
clear all
parentFolder = 'TMR';

% Get a list of all subfolders in the TMR folder
subfolders = dir(parentFolder);
subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));

% Define the border width in pixels
borderWidth = 3;

% Define font properties for the labels
fontSize = 50; % Font size for the labels
fontColor = [0, 0, 0]; % Black text

% Loop through each subfolder
for i = 1:length(subfolders)
    % Get the current subfolder path
    currentFolder = fullfile(parentFolder, subfolders(i).name);
    
    % Load images 'a', 'b', and 'c'
    imageA = imread(fullfile(currentFolder, 'a.jpg')); % Adjust extension as needed
    imageB = imread(fullfile(currentFolder, 'b.jpg'));
    imageC = imread(fullfile(currentFolder, 'c.jpg'));


    % %Convert if needed
    if size(imageA, 3) == 3 && size(imageB,3)==3 && size(imageC,3)==3 
            borderWidth=3;
    else
        imageA = gray2rgb(imageA);
        imageB = gray2rgb(imageB);
        imageC = gray2rgb(imageC);
     end

    
    % Find the minimum height and width across all images
    minHeight = min([size(imageA, 1), size(imageB, 1), size(imageC, 1)]);
    minWidth = min([size(imageA, 2), size(imageB, 2), size(imageC, 2)]);
    
    % Resize all images to the minimum dimensions
    imageA = imresize(imageA, [minHeight, minWidth]);
    imageB = imresize(imageB, [minHeight, minWidth]);
    imageC = imresize(imageC, [minHeight, minWidth]);
    
    % Create a white border
    whiteBorder = uint8(255 * ones(minHeight, borderWidth, 3)); % Assuming RGB images
    
    % Concatenate the images with borders in between
    montageImage = [imageA, whiteBorder, imageB, whiteBorder, imageC];
    %montageImage=imresize(montageImage,[400,700],"bilinear");

    minHeight = (size(montageImage,1));
    minWidth = (size(montageImage,2)-10)/3;
    
    % Add labels below the montage
    labelHeight = 110; % Height for the label area
    labelledCanvas = uint8(255 * ones(minHeight + labelHeight, size(montageImage, 2), 3));
    labelledCanvas(1:minHeight, :, :) = montageImage; % Add montage to the canvas
    labelledCanvas(minHeight:minHeight+2, :, :)=0;

    % Insert text labels '(a)', '(b)', '(c)'
    labelledCanvas = insertText(labelledCanvas, ...
        [minWidth/2, minHeight + 5; ... % Position for '(a)'
         minWidth + borderWidth + minWidth/2, minHeight + 5; ... % Position for '(b)'
         2 * (minWidth + borderWidth) + minWidth/2, minHeight + 5], ... % Position for '(c)'
        {'(a)', '(b)', '(c)'}, ...
        'FontSize', fontSize, ...
        'BoxColor', 'white', ...
        'BoxOpacity', 0, ...
        'TextColor', fontColor);
    
    % Display the labelled montage
    figure;
    imshow(labelledCanvas);
    title(['Montage of folder: ', subfolders(i).name]);
    
    % Save the labelled montage as a JPG file in the current subfolder
    outputFileName = fullfile(currentFolder, 'montage.jpg');
    imwrite(labelledCanvas, outputFileName);
    disp(['Montage saved as: ', outputFileName]);
end
