%% Compare All of the Images in the Comparison Folders Against Each Other:

% Define the host folder and subfolders
host_folder = 'Test3_comparisonqs\Comparisons'; % Replace with your folder path
subfolders = {'1', '3', '5', '8'};
output_folder = fullfile(host_folder, 'Test_Results');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Initialize a text file for the results
results_file = fullfile(output_folder, 'results.txt');
fileID = fopen(results_file, 'w');
fprintf(fileID, 'Folder\tImage\tMSE\tMAE\tSSIM\tMS-SSIM\tCableCheck\tPSNR-HVS-M\tPSNR-HVS\n');

% Loop through each subfolder
for i = 1:length(subfolders)
    folder_path = fullfile(host_folder, subfolders{i});
    
    % Get all image files of supported types
    image_files = [ ...
        dir(fullfile(folder_path, '*.bmp')); ...
        dir(fullfile(folder_path, '*.jpg')); ...
        dir(fullfile(folder_path, '*.png')) ...
    ];
    
    % Find the reference image ('grid' in filename, any format)
    ref_idx = find(contains(lower({image_files.name}), 'grid'), 1);
    if isempty(ref_idx)
        warning('No reference "grid" image found in %s', folder_path);
        continue;
    end

    % Load reference image and remove it from list
    og_image_path = fullfile(folder_path, image_files(ref_idx).name);
    original = im2double(imread(og_image_path));
    if size(original, 3) == 3
        original = rgb2gray(original);
    end
    image_files(ref_idx) = []; % remove reference from list

    % Loop through all other images in folder
    for j = 1:length(image_files)
        img_path = fullfile(folder_path, image_files(j).name);
        target = im2double(imread(img_path));
        if size(target, 3) == 3
            target = rgb2gray(target);
        end

        % Resize if needed
        if ~isequal(size(original), size(target))
            target = imresize(target, size(original));
        end

        % Metric calculations
        mse_val = immse(original, target);
        mae_val = mean(abs(original(:) - target(:)));
        ssim_val = ssim(target, original);
        ms_ssim_val = multissim(original, target);
        cable_val = cable_check(target, original);
        o8 = uint8(original * 255);
        t8 = uint8(target * 255);
        [p_hvs_m, p_hvs] = psnrhvsm(o8, t8);

        % Save to results file
        fprintf(fileID, '%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
            subfolders{i}, image_files(j).name, mse_val, mae_val, ssim_val, ...
            ms_ssim_val, cable_val, p_hvs_m, p_hvs);
    end
end

% Close file
fclose(fileID);
disp('Processing complete. Results saved.');
