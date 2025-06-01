%This code loops through every folder in the folder 'Images' looking for (in image/folder for e.g., looks for file: folder_grid.bmp)
%Then it takes it and performs the following hologram generation codes on
%it.
% Clear environment
close all;
clear;
clc;

% Define the main folder containing subfolders with images
mainFolder = "Images4";


% Get a list of all subfolders in the main folder
subfolders = dir(mainFolder);
subfolders = subfolders([subfolders.isdir]); % Keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Remove '.' and '..'

% Loop through each subfolder
for k = 10:length(subfolders)
    subfolderName = subfolders(k).name;
    subfolderPath = fullfile(mainFolder, subfolderName);

    % Construct the image filename "[folder]_grid.bmp"
    imageFileName = subfolderName + "_grid.bmp";
    imageFilePath = fullfile(subfolderPath, imageFileName);
    
    % % Check if the image file exists
    % if ~isfile(imageFilePath)
    %     warning("File not found: %s", imageFilePath);
    %     continue; % Skip to next iteration
    % end

    fprintf("Processing: %s\n", imageFilePath);
    
    %% ascend
    cd("C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\classic-cgh-algorithms\IQM_Tests")
    
    %% Load and Process Target Image
    ImageFileRead = imread(imageFilePath);

    % Convert to grayscale if needed
    if size(ImageFileRead, 3) == 1
        TargetImage = double(ImageFileRead) ./ 255;
    elseif size(ImageFileRead, 3) == 3
        warning("The target image has 3 channels, converting to grayscale by default")
        TargetImage = double(rgb2gray(ImageFileRead)) ./ 255;
    else
        error("Unexpected image format for: %s", imageFilePath);
    end

    TargetAmplitude = sqrt(TargetImage);
    %% descend
    cd(subfolderPath);
    %% Gerchberg-Saxton Algorithm
    N = 30; % Number of iterations
    E = TargetAmplitude .* exp(1i * rand(size(TargetAmplitude)) * 2 * pi);
    
    GS_NMSE_list = zeros(1, N);
    
    %  % Check if the final output file already exists
    % if isfile(fullfile(subfolderName + "_GS_Recon.bmp"))
    %     fprintf("Output file already exists: %s. Skipping the GS section.\n", subfolderName);
    %     continue; % Skip to the next iteration
    % end

    for iter = 1:N
        % Propagate and apply constraints
        A = ifftshift(ifft2(ifftshift(E)));
        A = exp(1i * angle(A));
        A = exp(1i * pi * double(angle(A) > pi/2 | angle(A) < -pi/2));

        % Propagate back
        E = fftshift(fft2(fftshift(A)));

        % Compute NMSE
        R = abs(E * sqrt(sum(TargetAmplitude(:).^2) / sum(abs(E(:)).^2))).^2;
        GS_NMSE_list(iter) = immse(R, TargetImage) / sum(TargetImage(:).^2);

        % Save final iteration images
        if iter == N
            imwrite(angle(A), fullfile(subfolderName + "_GS_Holo.bmp"), 'bmp');
            imwrite(R, fullfile(subfolderName + "_GS_Recon.bmp"), 'bmp');
        end
    end

    %% OSPR Algorithm
    N = 8;
    R_total = zeros(size(TargetAmplitude));
    FreemanHolo = zeros([size(TargetAmplitude), 3], 'uint8');
    

    output=(fullfile(subfolderName, "OSPR"));

    % if isfile(fullfile(subfolderName + + "_Freeman_OSPR_Holo_Frame_" + FreemanHoloChannelIndex + ".bmp"))
    %     fprintf("Output file already exists: %s. Skipping the OSPR section.\n", subfolderName);
    %     continue; % Skip to the next iteration
    % end

    mkdir(output) 
    
    for FreemanHoloChannelIndex = 1:3
        FreemanHoloPerChannel = 0;
        H_frame = zeros(size(TargetAmplitude));

        for PerChannelFrameIndex = 1:N
            holo_frame_i = (FreemanHoloChannelIndex - 1) * 8 + PerChannelFrameIndex;

            % Add random phase to the target
            E = TargetAmplitude .* exp(1i * 2 * pi * rand(size(TargetAmplitude)));

            % Compute the backward propagation
            A = ifft2(ifftshift(E));
            H = angle(A);
            H = double(H > 0);

            % Save individual binary hologram
            imwrite(H, fullfile(output, subfolderName + "_BinaryHolo_" + holo_frame_i + ".png"), 'png');

            FreemanHoloPerChannel = FreemanHoloPerChannel + H .* 2^(PerChannelFrameIndex-1);
            R = abs(fftshift(fft2(fftshift(exp(1i * H * pi))))) .^ 2;
            R = R .* sqrt(sum(abs(TargetImage(:)).^2) / sum(abs(R(:)).^2));

            % Save reconstructed frames
            imwrite(R, fullfile(output, subfolderName + "_Recon_" + holo_frame_i + ".png"), 'png');
            R_total = R_total + R;
            H_frame = H + H_frame;

            % Save combined H_frame for the current channel
            if PerChannelFrameIndex == N
                H_frame_norm = H_frame / max(H_frame(:));
                imwrite(uint8(H_frame_norm * 255), ...
                    fullfile(output, subfolderName + "_Freeman_OSPR_Holo_Frame_" + FreemanHoloChannelIndex + ".bmp"), 'bmp');
            end
        end

        % Store in final hologram
        FreemanHolo(:,:,FreemanHoloChannelIndex) = uint8(FreemanHoloPerChannel);
    end

    % Save averaged reconstruction
    R_avg = R_total / (N * 3);
    imwrite(R_avg, fullfile(output, subfolderName + "_ReconAvg.png"), 'png');

    % Save final hologram
    imwrite(FreemanHolo, fullfile(output, subfolderName + "_Freeman_OSPR_Holo.png"), 'png');
    imwrite(FreemanHolo, fullfile(output, subfolderName + "_Freeman_OSPR_Holo.bmp"), 'bmp');

    %% GS with Sinc Filter and Frequency Cut-off
     
    TargetAmplitude = sqrt(TargetImage);
    
    % if isfile(fullfile(subfolderName + "Holo_GSFiltered.png"))
    %     fprintf("Output file already exists: %s. Skipping the GS_Sinc section.\n", subfolderName);
    %     continue; % Skip to the next iteration
    % end
    % Define the sinc filter for phase modulation
    filter_size = 12; % Filter size (adjustable)
    cutoff_frequency = 550e9; % Normalized cutoff frequency (adjustable, 0 < cutoff_frequency < 0.5)
    
    % Create the 2D sinc filter
    [x, y] = meshgrid(linspace(-1, 1, filter_size), linspace(-1, 1, filter_size));
    r = sqrt(x.^2 + y.^2);
    sinc_filter = sinc(2 * cutoff_frequency * r);
    
    % Normalize the filter
    sinc_filter = sinc_filter / sum(sinc_filter(:));
    
    % Major loop for hologram generation
    N = 30; % Total number of iterations
    E = TargetAmplitude .* exp(1i * rand(size(TargetAmplitude)) * 2 * pi);
    
    for i = 1:N
        % Propagate from replay field to hologram aperture
        A = ifftshift(ifft2(ifftshift(E)));
    
        % Apply phase-only constraint with sinc-based phase modulation
        phase_modulation = angle(A);
        phase_modulation_filtered = imfilter(phase_modulation, sinc_filter, 'replicate');
        A = exp(1i * phase_modulation_filtered); % Re-apply amplitude constraint
    
        % Propagate from hologram aperture to replay field
        E = fftshift(fft2(fftshift(A)));
    
        % Save progress
        R = abs(E * sqrt(sum(TargetAmplitude(:).^2) / sum(abs(E(:)).^2))).^2;
        if i == N
            imwrite(angle(A), fullfile(subfolderName + "Recon_GSFiltered.png"));
            imwrite(R, fullfile(subfolderName + "Holo_GSFiltered.png"));
        end
    
        % Apply the target amplitude constraint at the replay field
        E = TargetAmplitude .* exp(1i * angle(E));
    end

    % %% Half-Dark GS, only for non-edge imgaes
    % 
    % 
    % % Modify target image to balance black & white pixels
    % white = numel(find(TargetImage == 1));
    % black = numel(find(TargetImage < 1));
    % wtb_ratio = white / (black + white);
    % 
    % wtb_threshold = 0.5; % Target black-white ratio (modifiable)
    % percentage_flip = 30; % Percentage of pixels to randomly flip
    % 
    % if white > black
    %     cond = 0;
    %     condother = 1;
    % else
    %     cond = 1;
    %     condother = 0;
    % end
    % 
    % if isfile("HalfBW_Target_Amplitude.bmp") || wtb_ratio<0.2
    %     fprintf("Output file already exists: %s. Skipping the OSPR section.\n", subfolderName);
    %     continue; % Skip to the next iteration
    % end
    % 
    % while wtb_ratio ~= wtb_threshold
    %     for i = 1:size(TargetImage, 1)
    %         for j = 1:size(TargetImage, 2)
    %             if TargetImage(i, j) == cond
    %                 if randi([1, 100]) <= percentage_flip
    %                     TargetImage(i, j) = condother;
    %                     white = numel(find(TargetImage == 1));
    %                     black = numel(find(TargetImage < 1));
    %                     wtb_ratio = white / (black + white);
    %                 end
    %             end
    %         end
    %     end
    % end
    % 
    % TargetAmplitude = sqrt(TargetImage);
    % 
    % %imwrite(TargetAmplitude,fullfile(subfolderName + "Recon_GSFiltered.png"))
    % imwrite(TargetAmplitude, "HalfBW_Target_Amplitude.bmp", 'bmp');
    % 
    % % Major loop for hologram generation
    % N = 30; % Total number of iterations
    % E = TargetAmplitude .* exp(1i * rand(size(TargetAmplitude)) * 2 * pi);
    % GS_NMSE_list = [];
    % 
    % for iter = 1:N
    %     % Propagate from replay field to hologram aperture
    %     A = ifftshift(ifft2(ifftshift(E)));
    % 
    %     % Apply phase-only constraint for hologram aperture
    %     A = exp(1i * angle(A));
    % 
    %     % Binary phase quantization (optional)
    %     A = exp(1i * pi * double(angle(A) > pi/2 | angle(A) < -pi/2)); % Rounds -pi/2 ~ pi/2 to 0 and otherwise to pi
    % 
    %     % Propagate from hologram aperture to replay field
    %     E = fftshift(fft2(fftshift(A)));
    % 
    %     % Save progress
    %     R = abs(E * sqrt(sum(TargetAmplitude(:).^2) / sum(abs(E(:)).^2))).^2;
    %     GS_NMSE_list(iter) = immse(R, TargetImage) / sum(TargetImage(:).^2);
    % 
    %     if iter == N
    %         imwrite(angle(A), "50BWGS_i_" + string(iter) + ".bmp");
    %         imwrite(R, "50BWGSReplay_i_" + string(iter) + ".bmp");
    %     end
    % 
    %     % Apply the target amplitude constraint at the replay field
    %     E = TargetAmplitude .* exp(1i * angle(E));
    % end
    % 
    % %% Half-Dark OSPR, only for non-edge images
    %  cd ..
    % 
    % output=(fullfile(subfolderName, "50BWOSPR"));
    % mkdir(output)
    % 
    % %  if isfile(fullfile("50BWOSPR/Freeman_OSPR_Holo_" + ImageFileName + ".png", 'png')) || wtb_ratio<0.2
    % %     fprintf("Output file already exists: %s. Skipping the OSPR section.\n", subfolderName);
    % %     continue; % Skip to the next iteration
    % % end
    % % Major loop for hologram generation
    % tic
    % N = 8; % Number of sub-frames to average per RGB channel
    % T = TargetAmplitude;
    % R_total = zeros(size(T)); % To store the accumulated reconstruction
    % H_frame = zeros(size(T));
    % 
    % cd(output);
    % % Freeman projector encoding loop
    % FreemanHolo = zeros([size(T), 3], 'uint8'); % Initialize 3-channel output
    % 
    % for FreemanHoloChannelIndex = 1:3
    %     FreemanHoloPerChannel = 0;
    %     H_frame = zeros(size(T));
    % 
    %     for PerChannelFrameIndex = 1:N
    %         holo_frame_i = (FreemanHoloChannelIndex - 1) * 8 + PerChannelFrameIndex;
    % 
    %         % Add random phase to the target
    %         E = T .* exp(1i * 2 * pi * rand(size(T)));
    % 
    %         % Compute backward propagation
    %         A = ifft2(ifftshift(E));
    % 
    %         % Phase-only constraint for SLM
    %         H = angle(A);
    % 
    %         % Binary phase quantization
    %         H = double(H > 0);
    % 
    %         % Save binary hologram frame
    %         imwrite(H, "50BWBinaryHologram_frame_" + holo_frame_i + ".png", 'png');
    % 
    %         % Encode Freeman hologram bit-plane by bit-plane
    %         FreemanHoloPerChannel = FreemanHoloPerChannel + H .* 2^(PerChannelFrameIndex - 1);
    % 
    %         % Compute reconstructed intensity
    %         R = abs(fftshift(fft2(fftshift(exp(1i * H * pi))))) .^ 2;
    %         R = R .* sqrt(sum(TargetImage(:).^2) / sum(R(:).^2));
    %         imwrite(R, "50BWOSPR/Recon_frame_" + holo_frame_i + ".png", 'png');
    % 
    %         % Accumulate reconstructions
    %         R_total = R_total + R;
    %         H_frame = H_frame + H;
    % 
    %         % Save combined frame for current channel
    %         if PerChannelFrameIndex == N
    %             H_frame_norm = H_frame / max(H_frame(:));
    %             imwrite(uint8(H_frame_norm * 255), ...
    %                     "50BWOSPR/Freeman_OSPR_Holo_Frame_" + ImageFileName + "_" + FreemanHoloChannelIndex + ".bmp", 'bmp');
    %         end
    %     end
    % 
    %     % Store encoded hologram channel
    %     FreemanHolo(:, :, FreemanHoloChannelIndex) = uint8(FreemanHoloPerChannel);
    % end
    % 
    % % Compute and save average reconstruction
    % R_avg = R_total / holo_frame_i;
    % imwrite(R_avg, "50BWOSPR/ReconAvg_" + holo_frame_i + "_" + ImageFileName + ".png", 'png');
    % 
    % % Compute image quality metrics
    % NMSE = immse(R_avg, TargetImage) / sum(TargetImage(:).^2);
    % SSIM_val = ssim(R_avg, TargetImage);
    % disp("NMSE: " + NMSE + ", SSIM: " + SSIM_val);
    % 
    % % Save Freeman hologram
    % imwrite(FreemanHolo, "50BWOSPR/Freeman_OSPR_Holo_" + ImageFileName + ".png", 'png');
    % imwrite(FreemanHolo, "50BWOSPR/Freeman_OSPR_Holo_" + ImageFileName + ".bmp", 'bmp');
  
end

fprintf("Processing complete!\n");





