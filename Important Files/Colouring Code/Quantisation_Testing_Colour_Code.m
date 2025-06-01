% QuantizerComparisonHolo.m
% Compare multiple quantizers in region-based CGH generation

clear; clc;

% --- Load grayscale image ---
filename = 'Biccy.jpg';
[~, name, ~] = fileparts(filename);
img = im2double(imread(filename));
if size(img, 3) > 1
    img = rgb2gray(img);
end

% % --- Embed into 4x4-style full canvas ---
% T_total = zeros(size(img) * 2);
% T_total(1:end/2, 1:end/2) = img;
% %T_total(end/2+1:end, end/2+1:end) = rot90(img, 2);
% img = T_total;


% --- Embed image onto top half of canvas with double the number of rows ---
T_total = zeros(size(img,1)*2, size(img,2));  % Double rows only
T_total(1:size(img,1), :) = img;              % Place original image at the top
img = T_total;



% --- Resize to 1024x1280 using bilinear interpolation ---
target_height = 1024;
target_width = 1280;
img = imresize(img, [target_height, target_width], 'bilinear');

% --- Image size AFTER resizing ---
[M, N] = size(img);
total_pixels = M * N;

% --- Frequency grid ---
[u, v] = meshgrid(linspace(-0.5, 0.5, N), linspace(-0.5, 0.5, M));
U = u;
V = v;

% --- Create output folder ---
outdir = 'Final Quantisation Tests_16 and 32';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

imwrite(img, "target.bmp")
% --- Set number of levels and define quantizers ---
% nLevels = 16;
% quantizer_names = {'linear', 'sqrt', 'gamma05','log9', 'log19', ...
%                    'log1', 'log4', 'log250', 'sinusoid_squared'};%'sinusoid', 'tangent', 'sqrt', 'gamma05',};
% quantizer_funcs = {
%     @(img) img;
%     @(img) sqrt(sqrt(img)); %sqrt
%     @(img) img .^ 0.5; %gamma05
%     @(img) log(1 + 9 * img) / log(10);
%     @(img) log(1 + img * 19) / log(100); %used to be called 'log 100' corrected that to 'log 19'
%     @(img) log(1 + img * 1) / log(100);
%     @(img) log(1 + img * 4) / log(100);
%     @(img) log(1 + img * 250) / log(100); %log 250
%     @(img) (sin(pi * img / 2)).^2;
%     %@(img) min(tan(pi * img * 0.49) / tan(pi * 0.49), 1)
% };
% 
% % --- Loop through each quantizer ---
% for q = 1:length(quantizer_names)
%     qname = quantizer_names{q};
%     qfunc = quantizer_funcs{q};
%     fprintf('Applying quantizer: %s\n', qname);
% 
%     % --- Apply quantizer and bin ---
%     mapped = qfunc(img);
%     quantized = floor(mapped * nLevels);
%     quantized = min(quantized, nLevels - 1);
% 
%     % --- Build hologram from quantized regions ---
%     H_partial = zeros(M, N);
%     count = 0;
% 
%     for i = 0:(nLevels - 1)
%         mask = quantized == i;
%         if ~any(mask(:)), continue; end
% 
%         CC = bwconncomp(mask, 8);
%         for regionIdx = 1:CC.NumObjects
%             pixelList = CC.PixelIdxList{regionIdx};
%             [rows, cols] = ind2sub([M, N], pixelList);
% 
%             if numel(pixelList) < 3
%                 continue;
%             end
% 
%             coords = [min(cols), min(rows); max(cols), max(rows)];
%             mean_intensity = mean(img(pixelList));
%             scale = 1 + numel(pixelList) / total_pixels;
%             amp = mean_intensity * scale;
% 
%             H = RECT_colour(amp, coords, U, V);
%             H_partial = H_partial + H;
%             count = count + 1;
%         end
%     end
%     Replay = (ifft2(H_partial));
% 
%     % % --- Compute replay ---
%     % H_total = H_partial;
%     % replay = ifft2(H_total);
%     % replay_abs = mat2gray(abs(replay));
%     % 
%     % % --- Save outputs ---
%     % prefix = sprintf('%s_%s_%dlevels', name, qname, nLevels);
%     % imwrite(replay_abs, fullfile(outdir, [prefix '_replay.bmp']));
%     % imwrite(mat2gray(abs(H_total)), fullfile(outdir, [prefix '_hologram.bmp']));
%     % 
%     % % --- Save binary phase version (0 or π) ---
%     % binary_phase = pi * (mod(angle(H_total), pi) > (pi / 2));
%     % imwrite(mat2gray(binary_phase), fullfile(outdir, [prefix '_binary_phase.bmp']));
% 
%     % --- Apply Jones encoding ---
% tau = pi;
% thetajones = pi/2;
% J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
%      -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];
% 
% Ex_mod = J(1,1) * H_partial + J(1,2) * H_partial;
% Ey_mod = J(2,1) * H_partial + J(2,2) * H_partial;
% Holo_mod = Ex_mod + Ey_mod;
% 
% % --- Convert to binary phase (0 or π) ---
% Holo_mod_bin = mod(angle(Holo_mod), pi) > (pi/2);
% Holo_mod_bin = Holo_mod_bin * pi;
% 
% % --- Compute replay field (FFT of original unbinarized H_partial) ---
% 
% 
% % --- Save outputs ---
% prefix = sprintf('%s_%s_%dlevels', name, qname, nLevels);
% imwrite(mat2gray(abs(Holo_mod_bin)), fullfile(outdir, [prefix '_binary_phase.bmp']));
% imwrite(mat2gray(abs(Replay)), fullfile(outdir, [prefix '_replay.bmp']));
% 
% % --- Save resized and flipped binary hologram ---
% Holo_mod_sized = imresize(Holo_mod_bin, [1024, 1280], 'bilinear');
% Holo_mod_sized = fliplr(Holo_mod_sized);
% imwrite(mat2gray(Holo_mod_sized), fullfile(outdir, [prefix '_Holo_Sized.bmp']));
% imwrite(mat2gray(Holo_mod_sized), fullfile(outdir, [prefix '_Holo_Sized.png']));
% 
% 
%     fprintf('Saved all outputs for quantizer: %s\n\n', qname);
% end
% 
% disp('🎉 All quantizer tests completed and saved to: quantizer tester');
% %%


%% version 2 for running efficiency:
quantization_levels = [4, 16];
quantizer_names = {'linear', 'log9', 'log19', ...
                   'log1', 'log250', 'sinusoid_squared'}; %sqrt', 'gamma05','log4'
quantizer_funcs = {
    @(img) img;
    %@(img) sqrt(sqrt(img));
    %@(img) img .^ 0.5;
    @(img) log(1 + 9 * img) / log(10);
    @(img) log(1 + img * 19) / log(100);
    @(img) log(1 + img * 1) / log(100);
    %@(img) log(1 + img * 4) / log(100);
    @(img) log(1 + img * 250) / log(100);
    @(img) (sin(pi * img / 2)).^2;
};

for nLevels = quantization_levels
    for q = 1:length(quantizer_names)
        ...
        % Everything in your quantizer loop remains the same
        % Just update the prefix line:
        
        
        qname = quantizer_names{q};
        qfunc = quantizer_funcs{q};
        fprintf('Applying quantizer: %s\n', qname);
    
        % --- Apply quantizer and bin ---
        mapped = qfunc(img);
        quantized = floor(mapped * nLevels);
        quantized = min(quantized, nLevels - 1);
    
        % --- Build hologram from quantized regions ---
        H_partial = zeros(M, N);
        count = 0;
    
        for i = 0:(nLevels - 1)
            mask = quantized == i;
            if ~any(mask(:)), continue; end
    
            CC = bwconncomp(mask, 8);
            for regionIdx = 1:CC.NumObjects
                pixelList = CC.PixelIdxList{regionIdx};
                [rows, cols] = ind2sub([M, N], pixelList);
    
                if numel(pixelList) < 3
                    continue;
                end
    
                coords = [min(cols), min(rows); max(cols), max(rows)];
                mean_intensity = mean(img(pixelList));
                scale = 1 + numel(pixelList) / total_pixels;
                amp = mean_intensity * scale;
    
                H = RECT_colour(amp, coords, U, V);
                H_partial = H_partial + H;
                count = count + 1;
            end
        end
        Replay = (ifft2(H_partial));
    
        % % --- Compute replay ---
        % H_total = H_partial;
        % replay = ifft2(H_total);
        % replay_abs = mat2gray(abs(replay));
        % 
        % % --- Save outputs ---
        % prefix = sprintf('%s_%s_%dlevels', name, qname, nLevels);
        % imwrite(replay_abs, fullfile(outdir, [prefix '_replay.bmp']));
        % imwrite(mat2gray(abs(H_total)), fullfile(outdir, [prefix '_hologram.bmp']));
        % 
        % % --- Save binary phase version (0 or π) ---
        % binary_phase = pi * (mod(angle(H_total), pi) > (pi / 2));
        % imwrite(mat2gray(binary_phase), fullfile(outdir, [prefix '_binary_phase.bmp']));
    
        % --- Apply Jones encoding ---
    tau = pi;
    thetajones = pi/2;
    J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
         -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];
    
    Ex_mod = J(1,1) * H_partial + J(1,2) * H_partial;
    Ey_mod = J(2,1) * H_partial + J(2,2) * H_partial;
    Holo_mod = Ex_mod + Ey_mod;
    
    % --- Convert to binary phase (0 or π) ---
    Holo_mod_bin = mod(angle(Holo_mod), pi) > (pi/2);
    Holo_mod_bin = Holo_mod_bin * pi;
    
    % --- Compute replay field (FFT of original unbinarized H_partial) ---
    
    
    % --- Save outputs ---
    prefix = sprintf('%s_%s_%dlevels', name, qname, nLevels);
    imwrite(mat2gray(abs(Holo_mod_bin)), fullfile(outdir, [prefix '_binary_phase.bmp']));
    imwrite(mat2gray(abs(Replay)), fullfile(outdir, [prefix '_replay.bmp']));
    
    % --- Save resized and flipped binary hologram ---
    Holo_mod_sized = imresize(Holo_mod_bin, [1024, 1280], 'bilinear');
    Holo_mod_sized = fliplr(Holo_mod_sized);
    imwrite(mat2gray(Holo_mod_sized), fullfile(outdir, [prefix '_Holo_Sized.bmp']));
    imwrite(mat2gray(Holo_mod_sized), fullfile(outdir, [prefix '_Holo_Sized.png']));





        % ADDITION: save BMP version of each output
        imwrite(mat2gray(abs(Replay)), fullfile(outdir, [prefix '_replay.bmp']));
        imwrite(mat2gray(abs(Holo_mod_bin)), fullfile(outdir, [prefix '_binary_phase.bmp']));

        Holo_mod_sized = imresize(Holo_mod_bin, [1024, 1280], 'bilinear');
        Holo_mod_sized = fliplr(Holo_mod_sized);
        imwrite(mat2gray(Holo_mod_sized), fullfile(outdir, [prefix '_Holo_Sized.bmp']));

        % === NEW: Random phase holograms ===
        num_rand_layers = 3;
        H_randoms = zeros(M, N);
       
        for k = 1:num_rand_layers
            rand_phase = exp(1j * 2 * pi * rand(M, N));
            H_rand = H_partial .* rand_phase;
            H_randoms = H_randoms + H_rand;
        end
        H_rand_phase = H_randoms / num_rand_layers;
        % --- Convert to binary phase (0 or π) ---
        Holo_rand_phase = mod(angle(H_rand_phase), pi) > (pi/2);
        Holo_rand_phase = Holo_rand_phase * pi;
        imwrite(mat2gray(abs(Holo_rand_phase)), fullfile(outdir, [prefix '_holo_rando_phase_avg.bmp']));

        fprintf('Saved all outputs for quantizer: %s, nLevels = %d\n\n', qname, nLevels);
    end
end

imwrite(img, fullfile(outdir, [prefix 'Actual Reconstructed Image.bmp']));