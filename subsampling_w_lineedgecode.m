clear; close all;

%% CONFIG
filename = "Easy.bmp";
Alpha = 150;
scalefac = 1;  % scaling the hologram resolution
subsample_factor = 2;  % downsample rate: 1 = all pixels, 2 = every other, etc.

%% OUTPUT FOLDER
[~, name, ~] = fileparts(filename);
name = name + "_Alpha" + num2str(Alpha) + "_scale" + num2str(scalefac) + "_sub" + num2str(subsample_factor);
outputFolder = fullfile(pwd, name);
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

%% STEP 1: Load Image and Prepare Binary Target
target = im2bw(imread(filename));
[orig_y, orig_x] = size(target);
grid_target = zeros(2*orig_y, 2*orig_x);
grid_target(1:orig_y, 1:orig_x) = target;
grid_target(end-orig_y+1:end, end-orig_x+1:end) = rot90(target, 2);
target = grid_target;
[num_y, num_x] = size(target);

%% STEP 2: Hologram Grid
scaled_x = round(num_x * scalefac);
scaled_y = round(num_y * scalefac);
[U, V] = meshgrid(linspace(0,1,scaled_x), linspace(0,1,scaled_y));
Holo = zeros(scaled_y, scaled_x);
ContributionCount = zeros(scaled_y, scaled_x);

%% STEP 3: Connected Components
CC = bwconncomp(target, 4);
StoredSegments = {};

for objIdx = 1:CC.NumObjects
    [yy, xx] = ind2sub(size(target), CC.PixelIdxList{objIdx});
    coords = [xx, yy];

    if size(coords,1) >= 3
        coords = ordered(coords);  % assumes you have this function

        % Subsample the coordinates
        subsampled = coords(1:subsample_factor:end, :);
        if size(subsampled,1) < 3
            subsampled = coords([1, round(end/2), end], :);  % fallback
        end

        A = subsampled(1, :);
        B = subsampled(round(end/2), :);
        C = subsampled(end, :);

        dAB = norm(A - B); dBC = norm(B - C); dAC = norm(A - C);
        if (dAB + dBC) < Alpha * dAC
            label = 'l';
            pts = [A(1), A(2), C(1), C(2), NaN, NaN];
        else
            label = 'c';
            pts = [A(1), A(2), B(1), B(2), C(1), C(2)];
        end

        StoredSegments{end+1} = struct('label', label, 'pts', pts);
    end
end

%% STEP 4: Apply Hologram Equations
for k = 1:length(StoredSegments)
    s = StoredSegments{k};
    pts = s.pts;

    if s.label == 'l' && ~isnan(pts(3))
        H = LineHolo(pts(1), pts(2), pts(3), pts(4), U, V);
    elseif s.label == 'c' && ~isnan(pts(5))
        H = ArcHolo(pts(1), pts(2), pts(3), pts(4), pts(5), pts(6), U, V, 10);
    else
        continue;
    end

    H(isnan(H)) = 0;
    Holo = Holo + H;
    ContributionCount(H ~= 0) = ContributionCount(H ~= 0) + 1;
end

Holo(ContributionCount > 0) = Holo(ContributionCount > 0) ./ ContributionCount(ContributionCount > 0);

%% STEP 5: Visualize and Save
Replay = ifft2(Holo);
Holo_mod = mod(angle(Holo), pi) > pi/2;
Holo_mod = pi * Holo_mod;

figure; imshow(mat2gray(abs(Holo))); title('Final Hologram');
figure; imshow(mat2gray(abs(Replay))); title('Replay Field');
figure; imshow(mat2gray(Holo_mod)); title('Binary Phase Hologram');

imwrite(mat2gray(Holo_mod), fullfile(outputFolder, "Binary_Hologram.bmp"));
imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "Replay_Field.bmp"));
