% Specify file paths
file1 = "C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\Vector Graphics\Edge_Equation_Maker\BiccyFINAL250_scalefac_1\Var_Holo_.mat"; %line
file2 = "C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\Vector Graphics\Edge_Equation_Maker\SinSquared_Only_Output\Biccy_Holo_part.mat"; %colour

% Load the variables and extract the actual variable names
data1_struct = load(file1);
data2_struct = load(file2);

% Get field names
fields1 = fieldnames(data1_struct);
fields2 = fieldnames(data2_struct);

% Extract the first variable from each file
var1 = data1_struct.(fields1{1});
var2 = data2_struct.(fields2{1});

var2=flipud(var2);
var2=fliplr(var2);
%var1 = var1 * (norm(var2(:)) / norm(var1(:)));

% Add them
combined = var1 + var2;
combined = 5e2* var1_norm + 2e3 * var2_norm;

% Compute the 2D Fourier Transform
FT_result = fft2(combined);




% Display the magnitude
%imshow(mat2gray(abs(fft2(var1))));
%imshow(mat2gray(abs(fft2(var2))));
imshow(mat2gray(abs(FT_result)));
imshow(mat2gray(abs(combined)));

imwrite(combined, "Hologram_colouredandlined.jpg")
imwrite(mat2gray(abs(FT_result)), "ReplayFieldcolouredandlined.jpg")
% Normalize combined
combined = 5e2 * var1_norm + 2e3 * var2_norm;

% Binary phase hologram: 0 if phase < 0, pi if phase ≥ 0
phase_bin = pi * (angle(combined) >= 0);

% Save as image (optional: convert to uint8 for proper saving)
imwrite(uint8(phase_bin / pi * 255), 'BinaryPhaseHologram.jpg');

% Display
figure; imshow(phase_bin, []); title('Binary Phase Hologram (0 / π)');
