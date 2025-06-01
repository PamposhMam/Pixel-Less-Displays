% run_CableCheck.m
% Script to run perceptual noise test using manual file paths

clc;
clear;

% --- Define file paths manually ---
file_ref  = "C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\classic-cgh-algorithms\IQM_Tests\Test2_comparisonqs\1\p_grid.bmp";    % ← update this
file_test = "C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\classic-cgh-algorithms\IQM_Tests\Test2_comparisonqs\1\pHolo_GSFiltered.png";         % ← update this

% --- Load the images ---
RefIm  = imread(file_ref);
TestIm = imread(file_test);

% --- Convert to grayscale if needed ---
if size(RefIm, 3) == 3
    RefIm = rgb2gray(RefIm);
end
if size(TestIm, 3) == 3
    TestIm = rgb2gray(TestIm);
end

% --- Resize to match if needed ---
if ~isequal(size(RefIm), size(TestIm))
    warning('Image sizes do not match. Resizing test image to match reference.');
    TestIm = imresize(TestIm, size(RefIm));
end

% --- Run the noise metric function ---
metric = Cable_check(TestIm, RefIm);

% --- Show results ---
fprintf('Cable_check perceptual noise metric: %.4f\n', metric);

threshold = 0.0237;
if metric >= threshold
    disp('Perceptual noise is SIGNIFICANT.');
else
    disp('Perceptual noise is NOT significant.');
end

% %%
% function [metric] = Cable_check(recon, init)
% % CABLE_CHECK estimates perceptual noise using equal pixel energy assumption.
% %
% % Inputs:
% %   recon - Reconstructed/test image
% %   init  - Reference/original image
% %
% % Output:
% %   metric - Noise perceptibility score
% %
% % Formula:
% %   metric = coverage * (1 - p) / p
% %   where:
% %       p = fraction of white pixels in the reference image
% %       coverage = fraction of white pixels in the test image
% 
%     % Normalize input images to [0,1]
%     I_test = mat2gray(recon);
%     I_ref  = mat2gray(init);
% 
%     % Total pixels
%     num_pixels = numel(I_ref);
% 
%     % Proportion of image energy from reference
%     p = sum(I_ref(:) > 0.5) / num_pixels;
% 
%     % Coverage from reconstructed image
%     coverage = sum(I_test(:) > 0.5) / num_pixels;
% 
%     % Final metric (with safe division)
%     metric = coverage * (1 - p) / (p + eps);
% end


function [metric] = Cable_check(recon, init)
%CABLE_CHECK Summary of this function goes here
%Noise is perceptually significant if c(1-p)/p>=0.0237
%p= fraction of energy that forms the image (p=sumofall(alpha^2*T^2))---- i assume we
%can assume that each pixel has the same amount of energy? so then surely
%that would just be ((1-%white)/%white)*coverage-- also background energy
%mean=1-p so yeah... (u_e)



% I= recon; %transfer target image over
% pixel_energy = abs(I).^2; %energy of each pixel grid-wise
% p=0;    %energy of the non-noise
% all=0; %all the energy (for normalisation)
% numsofall=numel(recon); %find the total number of pixels 
% white=numel(find(I)>0.5);  %number of white pixels
% black=numel(find(I)<0.5);   %number of black pixels
% c=white/(white+black);      %covergage. This, of course, assumes a white picture on a black surface (the scope of this experiment)
% 
% [row, col]= size(I);
% 
% for i=1:row
%     for j=1:col
%         all=all+pixel_energy(i,j);
%         if I(i,j)==1
%             p=p+pixel_energy(i,j);
%         end
% 
%     end
% end
% 
% 
% metric= c*(1-p)/p;

    % CABLE_CHECK calculates a noise perceptibility metric for a reconstructed image.
    % Inputs:
    %   recon - Reconstructed image (grayscale or binary).
    % Outputs:
    %   metric - Noise metric based on energy distribution.

    % Normalize the input image to range [0,1]
    I = mat2gray(recon);

    % Compute pixel energy
    pixel_energy = abs(init).^2;

    % Total number of pixels
    num_pixels = numel(I);

    % Classify pixels as white or black
    white_mask = (I > 0.5);
    black_mask = (I <= 0.5);

    % Calculate coverage (proportion of white pixels)
    coverage = sum(white_mask(:)) / num_pixels;

    % Compute total energy and white pixel energy
    total_energy = sum(pixel_energy(:));
    white_energy = sum(pixel_energy(white_mask));

    % Fraction of energy forming the image (p)
    p = white_energy / total_energy;

    % Compute the noise metric
    metric = coverage * (1 - p) / (p + eps);
end
