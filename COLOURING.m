%clear all
%rectftfunction
coords1 = [120, 650;
          120, 700;
          200, 700;
          200, 120]/2; % Rect function, requires 4 input-coordinates

coords2=[100, 650;
    100, 800;
    300, 800;
    300, 650]/2;





amp_list = [10, 8];

%------------------------------------------------------------------------

% Define grid size [1024x1280]
u = 1024;
v = 1280;


% Define spatial frequency grid (normalized)
[u_grid, v_grid] = meshgrid(...
    linspace(-0.5, 0.5, v), ...
    linspace(-0.5, 0.5, u));

% Initialise hologram
H_total = zeros(u, v);

% Pack into cell array and amplitude list
coords_list = {coords1, coords2};


% Loop to add each rectangle's FT
for i = 1:length(coords_list)
    coords = coords_list{i};

    % Convert to bounding box format
    bbox = [min(coords(:,1)), min(coords(:,2));
            max(coords(:,1)), max(coords(:,2))];

    % Call your RECT_colour function
    holo_i = RECT_colour(amp_list(i), bbox, u_grid, v_grid);

    % Superpose
    H_total = H_total + holo_i;
end

% Visualize
imagesc(abs(H_total)); colormap hot; axis image;
title('Superposed Rectangle Holograms');

replay=ifft2(H_total);
replay_save=mat2gray(abs(replay));
imshow(replay_save);
