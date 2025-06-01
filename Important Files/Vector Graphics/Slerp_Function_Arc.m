% %% Define Initial Params
% Holo=zeros(size(1280,1080));
% omega=pi/3; %define the angle of the arc
% x0= 10*cos(pi/3); y0= 10*sin(pi/3);
% x1= 2*cos(0); y1= 2*sin(0);
% 
% %when i do the actual code, omega will be calculated on it's own...
% 
% %% Calculation of parameters
% %alpha
% alpha= (sin((1-t)*omega)*x0 + sin(t*omega)*x1)/sin(omega);
% 
% %beta
% beta= (sin((1-t)*omega)*y0 + sin(t*omega)*y1)/sin(omega);
% v=linspace(numel(Holo(1,:)),);
% 
% 
% Holo_Deriv_Func= e^(2*pi*1j*(U.*alpha+ V.*beta));
% 
% 
% Holo_Func= Simpsons_Rule(Holo_Deriv_Funct, 0, 1, 10);

%%
% function integral = simpsons_rule(f, a, b, n)
%     % Ensure n is a positive integer
%     if n <= 0 || mod(n, 1) ~= 0
%         error('n must be a positive integer');
%     end
% 
%     % Use 2n intervals (ensuring even number of intervals)
%     num_intervals = 2 * n;
%     h = (b - a) / num_intervals;
% 
%     % Generate x values
%     x = linspace(a, b, num_intervals + 1);
% 
%     % Evaluate the function at x values
%     y = f(x);
% 
%     % Simpson's rule formula
%     integral = h / 3 * (y(1) + y(end) + 4 * sum(y(2:2:end-1)) + 2 * sum(y(3:2:end-2)));
% end

%% Define Initial Params
%Holo = zeros(1080, 1280);
clear all 


alias=1;
Holo=zeros(alias*1080,alias*1280);

[numrows, numcols]= size(Holo);
omega = pi-0.2 ; % define the angle of the arc
offset= 0.5;
r1=150;
r2=150;

x0 =  r1 * cos(0+omega+offset); 
y0 =  r1 * sin(0+omega+offset);
x1 =  r1 * cos(0);
y1 =  r1 * sin(0);

name= "splines_omega_" + num2str(omega) + "arcrad" + num2str(r1) + "_" + num2str(r2) + "_offset" + num2str(offset);

% Step 2: Create a Folder with the Same Name
outputFolder = fullfile(pwd, name); % Create folder path
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % Make directory if it doesn't exist
end



% when i do the actual code, omega will be calculated on its own...

%% Calculation of parameters {Me rehashing Jake}
% t values — evenly spaced points for Simpson’s rule
n = 90; % Number of subintervals
t = linspace(0, 1, 2 * n + 1); % 2n+1 points for Simpson’s Rule

% Numeric alpha and beta
alpha = @(t) (sin((1 - t) * omega) * x0 + sin(t * omega) * x1) / sin(omega);
beta = @(t) (sin((1 - t) * omega) * y0 + sin(t * omega) * y1) / sin(omega);

% %If want to use a for loop
% alpha = (sin((1 - t) * omega) * x0 + sin(t * omega) * x1) / sin(omega);
% beta = (sin((1 - t) * omega) * y0 + sin(t * omega) * y1) / sin(omega);



% U and V axes
U = linspace(-0.5, 0.5, size(Holo, 2)); %orignially -1 to 1
V = linspace(-0.5, 0.5, size(Holo, 1));
[U, V] = meshgrid(U, V);

% Holo Derivative Function — numeric and evaluated at each t
Holo_Deriv_Func = @(t) exp(2 * pi * 1j * (U .* alpha(t) + V .* beta(t)));

% Use Simpson's Rule to integrate over numeric t
Holo_Func = Simpsons_Rule(Holo_Deriv_Func, 0, 1, n);

Holo_Func=Holo_Func*r1*omega;


%(Changing the deriv function so i can use matlab inbuilt functions)
Holo_Deriv_Spare= @(x) exp(2 * pi * 1j * (U .* alpha(x) + V .* beta(x)));



% %Matlab Integrator
% Alt_holo = zeros(size(Holo)); % Initialize matrix
% 
% for i = 1:numrows
%     for j = 1:numcols
%         Holo_Deriv_Scalar = @(t) exp(2 * pi * 1j * (U(i,j) * alpha(t) + V(i,j) * beta(t)));
%         Alt_holo(i, j) = integral(Holo_Deriv_Scalar, 0, 1);
%     end
% end




% %% Calculation of parameters from FT of the thing
%  U = linspace(-1, 1, size(Holo, 2));
%  V = linspace(-1, 1, size(Holo, 1));
%  [U, V] = meshgrid(U, V);
% 
% %myHolo1 = (-1j/(2*sin(omega))*(x0*((exp(1j.*U)/(-omega+U))-(exp(1j*U)/(omega+U))-(exp(1j*omega)/(U-omega))+(exp(-1j*omega)/(omega+U))))+x1((exp(1j*(U+omega))/(U+omega))-(exp(-1j*(omega-U))/(U-omega)))-(1/(U+omega))+(1/(U-omega)));
% % Ensure omega and U are defined properly, and U does not have omega
% myHolo1 = (-1j / (2 * sin(omega))) * (x0 * (exp(1j * (U+V)) ./ (-omega + U+V) - exp(1j * (U+V)) ./ (omega + (U+V)) - exp(1j * omega) ./ (U+V - omega) + exp(-1j * omega) ./ (omega + U +V))) + ...
%            x0 * (exp(1j * (U+V + omega)) ./ (U+V + omega) - exp(-1j * (omega - U-V)) ./ (U+V - omega)) - ...
%            (1 ./ (V+U + omega)) + (1 ./ (V+U - omega));
% 
% 
% myHolo2 = (-1j / (2 * sin(omega))) * (y0 * (exp(1j * (V+U)) ./ (-omega + U+V) - exp(1j * (U+V)) ./ (omega + U+ V) - exp(1j * omega) ./ (V+U - omega) + exp(-1j * omega) ./ (omega + V+U))) + ...
%            y1 * (exp(1j * (V+U + omega)) ./ (V+U + omega) - exp(-1j * (omega - V-U)) ./ (V+U - omega)) - ...
%            (1 ./ (U+V + omega)) + (1 ./ (U+V - omega));
% 
% 
% 
% %Holo_Func = myHolo1 .* myHolo2;  % Combine multiplicatively in Fourier space
% %Holo_Func = sqrt(abs(myHolo1).^2 + abs(myHolo2).^2);
% 
% Holo_Func = exp(1j * angle(myHolo1 .* myHolo2)) .* sqrt(abs(myHolo1).^2 + abs(myHolo2).^2);
% 
% %Holo_Func=myHolo1+myHolo2;
% %Include a coordinate shift to wherever you want to start from

a= numcols/2; %coordinate along the x-axis
b= numrows/2; %coordinate along the y-axis

%%% save 

Holo = Holo + Holo_Func; %*exp(-2 * pi * 1j * (U .* a + V .* b));

% if ~isreal(HoloAlt)
%     HoloAlt = abs(HoloAlt); % Convert to real magnitude values
% end
% 
% 
% %HoloAlt=Holo+Alt_holo;
% imwrite(uint8(255*mat2gray(HoloAlt)), 'matlab_arc_integral.bmp')
% 
% replayAltArc=fft(fftshift(HoloAlt));
% if ~isreal(replayAltArc)
%     replayAltArc = abs(replayAltArc); % Convert to real magnitude values
% end
% 
% imwrite(uint8(255*mat2gray(replayAltArc)), 'AltArcReplay.bmp')

imwrite(mat2gray(abs(Holo)), 'Arc_Hologram.bmp') % If you want a BMP image

%% encoding for the hologram:
% Step 3: Define the Jones Matrix
tau = pi;
thetajones = pi/2;
J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
    -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];

% % Step 1: Apply the Jones Matrix Directly to Holo
Ex_mod = J(1,1) * Holo + J(1,2) * Holo;
Ey_mod = J(2,1) * Holo + J(2,2) * Holo;

% Step 2: Compute the Final Hologram (Before FFT)
Holo_mod = Ex_mod + Ey_mod;

% % Step 3: Convert to Binary Phase (0 or pi)
Holo_mod = mod(angle(Holo_mod), pi) > (pi/2);  % Threshold to 0 or 1
Holo_mod = Holo_mod * pi;  % Convert binary 0,1 → 0,π

% Step 4: Compute the Replay Field
Replay = fftshift(fft2(Holo)); % FFT of binary phase modulated hologram

% Step 5: Display Results
figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or π)');
figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT) from Binary Phase');

% Step 6: Save the Binary Hologram
imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "Binary_Phase_Hologram.bmp"));
imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "Binary_Replay.bmp"));

%Step 7: Also save a correct-sized hologram to the folder
 Holo_mod_sized= imresize(Holo_mod, [1024, 1280], 'bilinear');
% Holo_mod_

Holo_mod_sized=fliplr(Holo_mod_sized);
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.bmp"));
imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.png"));

%Holo=abs(Holo);

% Replay = fftshift(fft2(Holo));
% 
% imwrite(mat2gray(abs(Replay)), 'Replay_Arc.bmp') % If you want a BMP image





% %%
% % Step 1: Generate a Phase Mask for "JONES APERTURE"
% textSize = 50; % Adjust for visibility
% position = [round(numcols * 0.3), round(numrows * 0.5)]; % Center it
% 
% % Create a blank mask to overlay text
% textMask = ones(numrows, numcols); 
% textMask = insertText(textMask, position, 'JONES APERTURE', 'FontSize', textSize, 'BoxColor', 'black', 'TextColor', 'white');
% textMask = im2bw(rgb2gray(textMask)); % Convert to binary mask
% 
% % Convert binary mask into a phase modulation (0 or pi)
% phaseModulation = pi * textMask; 
% 
% %% random phase method
% 
% % % Step 2: Apply the Phase Modulation to the Hologram
% % Holo = Holo .* exp(1j * phaseModulation); 
% % 
% % %Apply random noise to every randomth pixel
% % 
% % % Generate 20 random phase values in the range [-pi, pi]
% % crayphase = (randi([1,3], 1, 20) * pi) - 2 * pi;
% % 
% % counter = 0;
% % randomness = randi([1,10]); % Single random number in [1,10]
% % nth = randomness;
% % 
% % % 
% % % for i = 1:numrows
% % %     for j = 1:numcols
% % %         counter = counter + 1;
% % %         if counter == nth
% % %             % Pick a random phase from crayphase
% % %             phase_index = randi([1, 20]); % Random index selection
% % %             Holo(i, j) = Holo(i, j) * exp(1j * crayphase(phase_index));
% % % 
% % %             % Reset counter and pick a new random interval
% % %             counter = 0;
% % %             nth = randi([1,10]); % New randomness for next iteration
% % %         end
% % %     end
% % end
% 
% %% radial cosine phase method: 
% % Step 2a: Apply the JONES APERTURE text-based phase modulation
% Holo = Holo .* exp(1j * phaseModulation); 
% 
% % Step 2b: Apply radial cosine phase modulation phi(u,v) = 1000π cos(sqrt(u² + v²))
% 
% % Create u, v coordinate grids centered at (0,0)
% [u, v] = meshgrid(linspace(-1,1,numcols), linspace(-1,1,numrows));  % normalized range
% 
% % Compute radial distance
% r = sqrt(u.^2 + v.^2);
% 
% % Compute phase
% phi = 1000 * pi * cos(r);  % This will oscillate rapidly
% 
% % Apply the additional phase modulation
% Holo = Holo .* exp(1j * phi);
% 
% 
% 
% % Step 3: Define the Jones Matrix
% tau = pi;
% thetajones = pi/2;
% J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
%     -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];
% 
% % % Step 1: Apply the Jones Matrix Directly to Holo
% Ex_mod = J(1,1) * Holo + J(1,2) * Holo;
% Ey_mod = J(2,1) * Holo + J(2,2) * Holo;
% 
% % Step 2: Compute the Final Hologram (Before FFT)
% Holo_mod = Ex_mod + Ey_mod;
% 
% % % Step 3: Convert to Binary Phase (0 or pi)
% Holo_mod = mod(angle(Holo_mod), pi) > (pi/2);  % Threshold to 0 or 1
% Holo_mod = Holo_mod * pi;  % Convert binary 0,1 → 0,π
% 
% % Step 4: Compute the Replay Field
% Replay = fftshift(fft2(Holo)); % FFT of binary phase modulated hologram
% 
% % Step 5: Display Results
% figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or π)');
% figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT) from Binary Phase');
% 
% % Step 6: Save the Binary Hologram
% imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "Binary_Phase_Hologram.bmp"));
% imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "Binary_Replay.bmp"));
% 
% %Step 7: Also save a correct-sized hologram to the folder
% Holo_mod_sized= imresize(Holo_mod, [1024, 1280], 'bilinear');
% % Holo_mod_
% 
% Holo_mod_sized=fliplr(Holo_mod_sized);
% imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.bmp"));
% imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "Holo_Sized.png"));
% 
% % %% and that sinc_convolved one
% % 
% % % Step 2: Apply the Phase Modulation to the Hologram
% % Holo_sinc = Holo_sinc .* exp(1j * phaseModulation); 
% % 
% % % Step 3: Define the Jones Matrix
% % tau = pi;
% % thetajones = pi/2;
% % J = [exp(-1j*tau/2)*(cos(thetajones))^2 + exp(1j*tau/2)*(sin(thetajones))^2, -1j*sin(tau/2)*sin(2*thetajones); ...
% %     -1j*sin(tau/2)*sin(2*thetajones), exp(1j*tau/2)*(cos(thetajones))^2 + exp(-1j*tau/2)*(sin(thetajones))^2];
% % 
% % % Step 1: Apply the Jones Matrix Directly to Holo_sinc
% % Ex_mod = J(1,1) * Holo_sinc + J(1,2) * Holo_sinc;
% % Ey_mod = J(2,1) * Holo_sinc + J(2,2) * Holo_sinc;
% % 
% % % Step 2: Compute the Final Hologram (Before FFT)
% % Holo_mod = Ex_mod + Ey_mod;
% % 
% % % Step 3: Convert to Binary Phase (0 or pi)
% % Holo_mod = mod(angle(Holo_mod), pi) > (pi/2);  % Threshold to 0 or 1
% % Holo_mod = Holo_mod * pi;  % Convert binary 0,1 → 0,π
% % 
% % % Step 4: Compute the Replay Field
% % Replay = fftshift(fft2(Holo_mod)); % FFT of binary phase modulated hologram
% % 
% % %%
% % Holo_mod=fliplr(Holo_mod);
% % 
% % % Step 5: Display Results
% % figure; imshow(mat2gray(abs(Holo_mod))); title('Binary Phase Hologram (0 or π)');
% % figure; imshow(mat2gray(abs(Replay))); title('Reconstructed Image (FFT) from Binary Phase');
% % 
% % % Step 6: Save the Binary Hologram
% % imwrite(mat2gray(abs(Holo_mod)), fullfile(outputFolder, "SincBinary_Phase_Hologram.bmp"));
% % imwrite(mat2gray(abs(Replay)), fullfile(outputFolder, "SincBinary_Replay.bmp"));
% % 
% % % Step 7: Also save a correct-sized hologram to the folder
% % Holo_mod_sized = fliplr(Holo_mod);
% % imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "SincHolo_Sized.bmp"));
% % imwrite(mat2gray(abs(Holo_mod_sized)), fullfile(outputFolder, "SincHolo_Sized.png"));
% % 
