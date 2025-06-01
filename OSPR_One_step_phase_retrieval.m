%%% Written for CMMPE by Daoming Dong, Youchao Wang and Jinze Sha
%%% Copyright 2018-2023

close all;
clear;
clc;

%% Set target image file
ImageFileName = "Arc_Hologram";
ImageFileExtension = ".bmp";
ImageFileNameFull = ImageFileName + ImageFileExtension;

%% Load and process target image file
ImageFileRead = imread(ImageFileNameFull);
% Check the number of channels of the target file, which should be 8 bit grayscale
if size(ImageFileRead, 3) == 1
    TargetImage = double(ImageFileRead) ./ 255;
elseif size(ImageFileRead, 3) == 3
    warning("The target image has 3 channels, converting to grayscale by default")
    TargetImage = double(rgb2gray(ImageFileRead)) ./ 255;
else
    error("The target image should have 1 channel, please check its bit depth")
end

%%% The following code might be useful, uncomment to use
% Image resize
TargetImage = imresize(TargetImage, [1024 1280], "bilinear");

% Horizontal flip
% TargetImage = flip(TargetImage, 2);

% Image normalization
% TargetImage = abs(ifft2(real(fft2(TargetImage))));
% TargetImage = TargetImage/max(TargetImage(:));

% Add upside down replica below the target
% T_rot = rot90(TargetImage, 2);
% TargetImage = [T_rot; TargetImage];
%%%


TargetAmplitude = sqrt(TargetImage); % target amplitude is the square root of intensity

%% Create output file folder for OSPR
if not(isfolder("OSPR_output"))
    mkdir("OSPR_output")
end

imwrite(TargetImage, "OSPR_output/TargetIntensity.png",'png');
imwrite(TargetAmplitude, "OSPR_output/TargetAmplitude.png",'png');

%% Major loop for hologram generation
tic
N = 8; % Number of sub-frames to average (per RGB encoding channel)
T = TargetAmplitude;
R_total = zeros(size(T)); % To get the averaged recon (for debug purpose)
H_frame=zeros(size(T));

% Assumming we want to output three channels in the Freeman projector
% Loop for Freeman hologram generation
for FreemanHoloChannelIndex = 1:3
    FreemanHoloPerChannel = 0;
    H_frame = zeros(size(T)); % Reset R_frame at the start of each channel
    for PerChannelFrameIndex = 1:N
        % Compute current subframe index
        holo_frame_i = (FreemanHoloChannelIndex-1) * 8 + PerChannelFrameIndex;

        % Add random phase to the target
        E = T .* exp(1i * 2 * pi * rand(size(T)));

        % Compute the backward propagation from the target to hologram plane
        A = (ifft2(ifftshift(E)));

        % Get the phase of hologram for phase-only SLM
        H = angle(A);

        % Binary phase quantization as constrained by the SLM
        H = double(H > 0);

        % Save individual frame of binary hologram (for debug purpose)
        imwrite(H, "OSPR_output/BinaryHologram_frame_" + holo_frame_i + ".png", 'png');

        % Encode the Freeman hologram bit-plane by bit-plane
        FreemanHoloPerChannel = FreemanHoloPerChannel + H .* 2^(PerChannelFrameIndex-1);

        % Calculate total reconstruction (for debug purpose)
        R = abs(fftshift(fft2(fftshift(exp(1i * H * pi))))) .^ 2;
        R = R .* sqrt(sum(abs(TargetImage(:)).^2) / sum(abs(R(:)).^2));
        imwrite(R, "OSPR_output/Recon_frame_" + holo_frame_i + ".png", 'png');
        R_total = R_total + R;
        H_frame = H + H_frame; % Accumulate subframe reconstructions

        % Save combined H_frame for the current channel
        if PerChannelFrameIndex == N
            % Normalize H_frame to [0, 255] before saving
            H_frame_norm = H_frame / max(H_frame(:));
            imwrite(uint8(H_frame_norm * 255), ...
                    "OSPR_output/Freeman_OSPR_Holo_Frame" + ImageFileName + "_ " + FreemanHoloChannelIndex + ".bmp", 'bmp');
        end
    end

    % Encode the final Freeman hologram
    FreemanHolo(:,:,FreemanHoloChannelIndex) = uint8(FreemanHoloPerChannel);
end


toc

R_avg = R_total/holo_frame_i;

%% Save average reconstruction
imwrite(R_avg, "OSPR_output/ReconAvg_" + holo_frame_i + "_" + ImageFileName + ".png",'png');
immse(R_avg, TargetImage) / sum(TargetImage(:).^2)
ssim(R_avg, TargetImage)


%% Save the hologram encoded for the Freeman projector
imwrite(FreemanHolo, "OSPR_output/Freeman_OSPR_Holo_" + ImageFileName + ".png",'png');

imwrite(FreemanHolo, "OSPR_output/Freeman_OSPR_Holo_" + ImageFileName + ".bmp",'bmp');






% %% Bitmap Combination- Pamposh
% output_folder = "GS_output";
% num_images = holo_frame_i; % Total number of iteration images generated
% 
% % Read the first image to get the size
% sample_image = imread(fullfile(output_folder, "B_holo_i_1.bmp"));
% [image_height, image_width] = size(sample_image);
% 
% % Initialize a matrix to store the cumulative sum
% combined_image = zeros(image_height, image_width);
% 
% % Sum all images
% for i = 1:num_images
%     img = double(imread(fullfile(output_folder, "GS_holo_i_" + string(i) + ".bmp")));
%     combined_image = combined_image + img;
% end
% 
% % Average the pixel values by dividing by the number of images
% combined_image = combined_image / num_images;
% 
% % Convert back to uint8 format and save
% combined_image = uint8(combined_image);
% imwrite(combined_image, fullfile(output_folder, "GS_combined_layered.bmp"), 'bmp');