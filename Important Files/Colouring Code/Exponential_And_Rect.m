% Multiple windowed exponentials (rect * exp) and their shifted sinc FTs
clear; clc;

% --- Spatial grid ---
N = 512;
L = 1;  % spatial width of rect
x = linspace(-1, 1, N);
y = linspace(-1, 1, N);
[XX, YY] = meshgrid(x, y);
dx = x(2) - x(1);  % spatial resolution

% --- Frequency grid ---
u = linspace(-1/(2*dx), 1/(2*dx), N);
v = linspace(-1/(2*dx), 1/(2*dx), N);
[UU, VV] = meshgrid(u, v);

% --- Define list of (alpha, beta) frequency centers ---
freq_coords = [ 0.0,  0.0;   % center
                200,  200;   % diagonal up-right
               -200,  200;   % diagonal up-left
                300, -300];  % diagonal down-right

% --- Initialize total function and FT ---
f_total = zeros(N, N);
F_total = zeros(N, N);

for k = 1:size(freq_coords, 1)
    alpha = freq_coords(k, 1);  % frequency shift in x
    beta  = freq_coords(k, 2);  % frequency shift in y

    % Spatial function: windowed plane wave (Rect × exp)
    window = double(abs(XX) <= L/2 & abs(YY) <= L/2);
    f_k = window .* exp(2j * pi * (alpha * XX + beta * YY));

    % Fourier transform (no shift, normalized)
    % F_k = fft2(f_k);
    % F_k = fftshift(F_k);  % for viewing centered frequencies
    sincfunc=
    expon=2*pi*(v-w0)

    % Combine
    f_total = f_total + f_k;
    F_total = F_total + F_k;
end


replay=mat2gray(abs(ifft2(f_total)));
% --- Plot result ---
figure;
subplot(2,2,1);
imshow(replay, []); title('Replay Field of f(x, y)');
subplot(2,2,2);
imshow(abs(f_total), []); title('|f(x, y)|');

subplot(2,2,3);
imshow(log(1 + abs(F_total)), []); title('log|F(u, v)|');
subplot(2,2,4);
imshow(angle(F_total), []); title('Phase of F(u, v)');
