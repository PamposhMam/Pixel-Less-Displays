function [holo] = RectXPONential_Gaussian(amp, coords, u, v, sigma_u, sigma_v)
% Fourier-space hologram using a 2D Gaussian envelope instead of a delta

% Rect dimensions
B2 = max(coords(:,1));
A2 = max(coords(:,2));
B1 = min(coords(:,1));
A1 = min(coords(:,2));
A = A2 - A1;
B = B2 - B1;
centre = [A2 - A/2, B2 - B/2];
decu, decv= centre;

% Gaussian envelope centered at (decu, decv)
gauss = amp * exp( -((u - decu).^2)/(2*sigma_u^2) - ((v - decv).^2)/(2*sigma_v^2) );

% Hologram function (smooth version)
holo = A * B ...
    .* sinc(pi*A*v) ...
    .* sinc(pi*B*u) ...
    .* exp(-2j*pi*(centre(1)*v + centre(2)*u)) ...
    .* gauss;
end
