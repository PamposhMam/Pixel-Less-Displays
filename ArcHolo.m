% function Holo_Func = ArcHolo(x0, y0, x1, y1, x2, y2, U, V, n)
%     % Compute omega (angle between vectors)
%     Ax = x1 - x0; Ay = y1 - y0;  
%     Bx = x2 - x1; By = y2 - y1;  
%     omega = acos((Ax * Bx + Ay * By) / (sqrt(Ax^2 + Ay^2) * sqrt(Bx^2 + By^2)));
% 
%     % Define t values — evenly spaced points for Simpson’s rule
%     t = linspace(0, 1, 2 * n + 1); % 2n+1 points for Simpson’s Rule
% 
%     % Compute α(t) and β(t) for arc trajectory
%     alpha = @(t) (sin((1 - t) * omega) * x0 + sin(t * omega) * x2) / sin(omega);
%     beta = @(t) (sin((1 - t) * omega) * y0 + sin(t * omega) * y2) / sin(omega);
% 
%     % Define the holographic function for integration
%     Holo_Deriv_Func = @(t) exp(2 * pi * 1j * (U .* alpha(t) + V .* beta(t)));
% 
%     % Perform numerical integration using Simpson’s Rule
%     Holo_Func = Simpsons_Rule(Holo_Deriv_Func, 0, 1, n);
% 
%         % Debugging: Check if Holo_Func has significant values
%     disp(['Max ArcHolo value for (', num2str(x0), ', ', num2str(y0), ') → ', num2str(max(abs(Holo_Func(:))))]);
% 
% end
%%
function Holo_Func = ArcHolo(x0, y0, x1, y1, x2, y2, U, V, n)
    % Define the upscaling factor
    scale_factor = 1;  

    % Expand grid resolution
    high_res_rows = size(U, 1) * scale_factor;
    high_res_cols = size(U, 2) * scale_factor;

    % Create high-res U and V grids
    U_high = linspace(min(U(:)), max(U(:)), high_res_cols);
    V_high = linspace(min(V(:)), max(V(:)), high_res_rows);
    [U_high, V_high] = meshgrid(U_high, V_high);

   % Scale arc coordinates directly
x0_high = x0 * scale_factor;
y0_high = y0 * scale_factor;
x1_high = x1 * scale_factor;
y1_high = y1 * scale_factor;
x2_high = x2 * scale_factor;
y2_high = y2 * scale_factor;

% Compute the perpendicular bisectors to find the circle center
mid1_x = (x0_high + x1_high) / 2;
mid1_y = (y0_high + y1_high) / 2;
mid2_x = (x1_high + x2_high) / 2;
mid2_y = (y1_high + y2_high) / 2;

% Slopes of the lines
dx1 = x1_high - x0_high;
dy1 = y1_high - y0_high;
dx2 = x2_high - x1_high;
dy2 = y2_high - y1_high;

% Perpendicular slopes
if dy1 == 0
    slope1 = Inf;
else
    slope1 = -dx1 / dy1;
end

if dy2 == 0
    slope2 = Inf;
else
    slope2 = -dx2 / dy2;
end

% Solve for intersection (center of circle)
if isinf(slope1) % Vertical case
    xc = mid1_x;
    yc = slope2 * (xc - mid2_x) + mid2_y;
elseif isinf(slope2)
    xc = mid2_x;
    yc = slope1 * (xc - mid1_x) + mid1_y;
else
    A = [1, -slope1; 1, -slope2];
    B = [mid1_y - slope1 * mid1_x; mid2_y - slope2 * mid2_x];
    center = A \ B;
    xc = center(1);
    yc = center(2);
end

% Compute angles
theta0 = atan2(y0_high - yc, x0_high - xc);
theta1 = atan2(y1_high - yc, x1_high - xc);
theta2 = atan2(y2_high - yc, x2_high - xc);

% Compute the angle turned through
omega = mod(theta2 - theta0, 2 * pi);

% Ensure direction is correct
if (theta1 < theta0 && theta2 > theta1) || (theta1 > theta0 && theta2 < theta1)
    omega = 2 * pi - omega;
end



    % Define t values — evenly spaced points for Simpson’s rule
    t = linspace(0, 1, 2 * n + 1); % 2n+1 points for Simpson’s Rule
    
    % Compute α(t) and β(t) for arc trajectory (now in high-res grid)
    alpha = @(t) (sin((1 - t) * omega) * x0_high + sin(t * omega) * x2_high) / sin(omega);
    beta = @(t) (sin((1 - t) * omega) * y0_high + sin(t * omega) * y2_high) / sin(omega);
    
    % Define the holographic function for integration
    Holo_Deriv_Func = @(t) exp(2 * pi * 1j * (U_high .* alpha(t) + V_high .* beta(t)));
    

    % Perform numerical integration using Simpson’s Rule
    Holo_HighRes = Simpsons_Rule(Holo_Deriv_Func, 0, 1, n);
    Holo_HighRes(isnan(Holo_HighRes)) = 0;

    % Downsample back to original grid size
    Holo_Func = imresize(Holo_HighRes, [size(U,1), size(U,2)], 'bilinear');
    
       % Compute arc radius (rho)
    rho = sqrt((x0_high - xc)^2 + (y0_high - yc)^2);
    
    % Compute arc length
    scale_mag = rho * omega;
    % Chord length between start and end
    chord_length = sqrt((x2_high - x0_high)^2 + (y2_high - y0_high)^2);

    % Clamp arc length if it exceeds 3× chord length
    if scale_mag > 3 * chord_length
        disp(['⚠️ Clamped scale_mag from ', num2str(scale_mag), ' to ', num2str(3 * chord_length)]);
        scale_mag = 3 * chord_length;
    end



    %Normalize amplitude:
    Holo_Func=Holo_Func*scale_mag; %Delete this line if it performs worse...

    % Debugging: Check if Holo_Func has significant values
    disp(['Max ArcHolo value for (', num2str(x0), ', ', num2str(y0), ') → ', num2str(max(abs(Holo_Func(:))))]);
end
