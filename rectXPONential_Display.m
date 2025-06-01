% Display_Hologram_Overlap.m
% Visualizes overlapped Gaussian and Triangle Fourier-space holograms

% -------------------- USER INPUTS --------------------
% Coordinates: [x, y] pairs defining the rectangle
coords1 = [120 20; 210 100];  % Example: corners of a rectangle
coords2 = [120 120; 210 200];

% Resolution of u, v grid
N = 800;
[u, v] = meshgrid(linspace(-0.5, 0.5, N), linspace(-0.5, 0.5, N));

% Amplitude and envelope parameters
sigma_u = 50;
sigma_v = 10;
halfwidth_u = 10;
halfwidth_v = 10;

% -----------------------------------------------------

% Generate Gaussian hologram
%H_gauss = RectXPONential_Gaussian(coords1, u, v, sigma_u, sigma_v);

% Generate Triangle hologram (using rectangle center as default shift)
H_tri = RectXPONential_Triangle(coords2, u, v);

H=H_tri;

H = H / max(abs(H(:)));  % Normalize to peak magnitude = 1

RPF=ifft2(H);
imshow(mat2gray(abs(RPF)));







%%



% clear all
% % Heatmap of replay field brightness for σ_u and σ_v from 0 to 24 in steps of 8
% 
% % -------------------- GRID SETUP --------------------
% N = 600;
% [u, v] = meshgrid(linspace(-0.5, 0.5, N), linspace(-0.5, 0.5, N));
% 
% % -------------------- COORDINATES --------------------
% coords = [10 10; 250 10; 250 250 ; 10 250 ];  % Rectangle corners
% 
% % -------------------- SIGMA RANGES --------------------
% sigma_vals = 0:20:60;
% 
% % -------------------- FIGURE --------------------
% figure('Name', '|RPF| Heatmaps for σ_u and σ_v (0:20:60)');
% 
% for i = 1:length(sigma_vals)
%     for j = 1:length(sigma_vals)
%         sigma_u = sigma_vals(j);
%         sigma_v = sigma_vals(i);
% 
%         % Skip σ = 0 case to avoid division by zero
%         if sigma_u == 0 || sigma_v == 0
%             RPF = zeros(N);
%         else
%             % Generate hologram
%             H = RectXPONential_Gaussian(coords, u, v, sigma_u, sigma_v);
% 
%             % Normalize safely
%             if max(abs(H(:))) > 0
%                 H = H / max(abs(H(:)));
%             end
% 
%             RPF = abs(ifft2(H));
%         end
% 
%         % Plot heatmap
%         subplot(length(sigma_vals), length(sigma_vals), (i - 1) * length(sigma_vals) + j);
%         % imagesc(RPF);
%         %colormap hot;
%         imagesc(RPF);
%         colormap gray;
%         axis off; axis image;
%         title(sprintf('\\sigma_u = %d, \\sigma_v = %d', sigma_u, sigma_v), 'FontSize', 10);
% 
%     end
% end
% sgtitle('Replay Fields for \sigma_u and \sigma_v (0:20:60)', 'FontWeight', 'bold', 'FontSize', 30);
% 
% 
% % Add single colorbar
% %h = colorbar('Position', [0.93 0.1 0.015 0.8]);
% %h.Label.String = '|RPF| Intensity';
% 
