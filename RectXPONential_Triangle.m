function [holo] = RectXPONential_Triangle(coords, u, v)
% Fourier-space hologram using a 2D triangle envelope instead of a delta

% Rect dimensions
B2 = max(coords(:,1));
A2 = max(coords(:,2));
B1 = min(coords(:,1));
A1 = min(coords(:,2));
A = A2 - A1;
B = B2 - B1;
centre = [A2 - A/2, B2 - B/2];

% Triangle envelope centered at (decu, decv)
tri_u = sinc(pi*u * B ).^2;
tri_v = sinc(pi*v * A ).^2;
triangle = A*B*tri_u .* tri_v.* exp(-2j * pi * (centre(1) * v + centre(2) * u));

% Hologram function (sharper version)
holo = (A * B)^2 ...
    .* sinc(pi * A * v) ...
    .* sinc(pi * B * u) ...
    .* exp(-2j * pi * (centre(1) * v + centre(2) * u)) ...
    .* triangle;
end
