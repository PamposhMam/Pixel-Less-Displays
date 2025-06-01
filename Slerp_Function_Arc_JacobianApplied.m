%% Define Initial Params
%Holo = zeros(1080, 1280);
clear all 

Holo=zeros(1280,1080);
[numrows, numcols]= size(Holo);
omega = pi/4 ; % define the angle of the arc
x0 =  10 * cos(0-omega); 
y0 =  10 * sin(0-omega);
x1 =  10 * cos(0);
y1 =  10 * sin(0);


% when i do the actual code, omega will be calculated on its own...

%% Calculation of parameters {Me rehashing Jake}

% U and V axes
U = linspace(-1, 1, size(Holo, 2));
V = linspace(-1, 1, size(Holo, 1));
[U, V] = meshgrid(U, V);

%Equation
gamma= (x0+x1)/sin(omega);
zeta= (y0+y1)/sin(omega);

Holo_Func= ((gamma*zeta)/(2*pi*1j))*exp(2*pi*1j*((x1-x0).*U + (y1-y0).*V));

a= numcols/2; %coordinate along the x-axis
b= numrows/2; %coordinate along the y-axis

%%% save 

Holo = Holo + Holo_Func; %*exp(-2 * pi * 1j * (U .* a + V .* b));

imwrite(mat2gray(abs(Holo)), 'Arc_Hologram.bmp') % If you want a BMP image

%Holo=angle(Holo);
%Holo=abs(Holo);

Replay = fftshift(fft(Holo));

imwrite(mat2gray(abs(Replay)), 'Replay_Arc.bmp') % If you want a BMP image
