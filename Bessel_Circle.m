% Parameters for aperture
grid_size = 512;             % Grid size for image
[x, y] = meshgrid(6:255, -256:255);  % Spatial grid
alpha=0.3; %radius
r = sqrt(x.^2 + y.^2);        % Radial distance from center

% Define aperture using first-order Bessel function
% Define aperture using first-order Bessel function divided by r
% bessel_aperture = alpha * besselj(1, 2*pi*alpha*r) ./ r;
% bessel_aperture= bessel_aperture -0.9*alpha * besselj(1, 2*pi*alpha*r) ./ r;
% 
% % Handle the singularity at r = 0
% bessel_aperture(r == 0) = alpha * pi; % Assign a finite value at r = 0
% 
% %Might not want angles outside of theta1-theta2. Add an eraser function
% %First define the angles i'd like to keep:
% theta0=pi/4; theta1=0.01;
% erasermod=2*alpha;
% %thetas= linspace(theta0, theta1, 1e6);
% x0=erasermod*cos(theta0); y0=erasermod*sin(theta0);
% x1=erasermod*cos(theta1); y1=erasermod*sin(theta1);
% centrex= grid_size/2;
% centrey=grid_size/2;

%Eraser function can be any shape that would get rid of what we want. I'm
%picking a rectangle
%eraser= a*b*sinc(pi*a*x)*sinc(pi*b*y);

% %% Erasing the unwanted sections of the circle:
% 
% if (0 <= theta0) && (theta0 < pi/2)
%     if (0 <= theta1) && (theta1 < pi/2) && theta1<theta0
%         a=x1-x0;
%         b=y0-y1;
%         remove=Rect_Eraser(a,b,x,y);
%         remove=remove*exp(-1j*2*pi*((x1+x0)*x/2+(y1+y0)*y/2));
%         %bessel_aperture=bessel_aperture-remove;
%     elseif (pi/2 >= theta1) && (theta1 < pi)
%         % Code for this case
%     elseif (pi >= theta1) && (theta1 < 3*pi/2)
%         % Code for this case
%     elseif (3*pi/2 >= theta1) && (theta1 < 2*pi)
%         % Code for this case
%     end
% elseif (pi/2 >= theta0) && (theta0 < pi)
%     if (0 >= theta1) && (theta1 < pi/2)
%         % Code for this case
%     elseif (pi/2 >= theta1) && (theta1 < pi)
%         % Code for this case
%     elseif (pi >= theta1) && (theta1 < 3*pi/2)
%         % Code for this case
%     elseif (3*pi/2 >= theta1) && (theta1 < 2*pi)
%         % Code for this case
%     end
% elseif (pi >= theta0) && (theta0 < 3*pi/2)
%     if (0 >= theta1) && (theta1 < pi/2)
%         % Code for this case
%     elseif (pi/2 >= theta1) && (theta1 < pi)
%         % Code for this case
%     elseif (pi >= theta1) && (theta1 < 3*pi/2)
%         % Code for this case
%     elseif (3*pi/2 >= theta1) && (theta1 < 2*pi)
%         % Code for this case
%     end
% elseif (3*pi/2 >= theta0) && (theta0 < 2*pi)
%     if (0 >= theta1) && (theta1 < pi/2)
%         % Code for this case
%     elseif (pi/2 >= theta1) && (theta1 < pi)
%         % Code for this case
%     elseif (pi >= theta1) && (theta1 < 3*pi/2)
%         % Code for this case
%     elseif (3*pi/2 >= theta1) && (theta1 < 2*pi)
%         % Code for this case
%     end
% end


%%
num_rows = 400; % Number of pixels along the vertical dimension
num_cols = 350;  % Number of pixels along the horizontal dimension

% Create the spatial coordinate arrays
u = linspace(-num_cols / 2, num_cols / 2, num_cols); % Horizontal pixel indices
v = linspace(-num_rows / 2, num_rows / 2, num_rows); % Vertical pixel indices
[X, Y] = meshgrid(u, v);

aperture=1*pi*cos(1*sqrt(X.^2+Y.^2));

% Fourier transform to compute the hologram
hologram = abs(fftshift(fft2(aperture)));

% Normalize for display
hologram = hologram / max(hologram(:));
% 
% % Display results
% figure;
% subplot(1, 2, 1);
% % imshow(bessel_aperture, []);
% % title('Aperture (Bessel Function)');

imshow(hologram, []);
title('Hologram (Arc from Fourier Transform)');
