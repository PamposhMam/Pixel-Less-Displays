function [holo] = RectXPONential_Gaussian(coords, u, v, sigma_u, sigma_v)
% Fourier-space hologram using a 2D Gaussian envelope instead of a delta

% Rect dimensions
B2 = max(coords(:,1));
A2 = max(coords(:,2));
B1 = min(coords(:,1));
A1 = min(coords(:,2));
A = A2 - A1;
B = B2 - B1;
centre = [A2 - A/2, B2 - B/2];
decu=centre(2);
decv=centre(1);

gauss = exp( -2*pi^2*sigma_u^2 .* u.^2 - 2*pi^2*sigma_v^2 .* v.^2 ) ...
    .* exp( -2j*pi * (decu .* u + decv .* v) ); %*2*pi*sigma_u*sigma_v ...;

% Hologram function (smooth version)
holo = A * B ...
    .* sinc(pi*A*v) ...
    .* sinc(pi*B*u) ...
    .* exp(-2j*pi*(centre(1)*v + centre(2)*u)) ...
    .* gauss;
end
