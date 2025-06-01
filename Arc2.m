% function Holo_Func = Arc2(x0, y0, x1, y1, x2, y2, U, V, n)
%     % Compute the circle center
%     [xc, yc, omega] = ComputeCircleCenter(x0, y0, x1, y1, x2, y2);
% 
%     % If omega > π/2, split into multiple sub-arcs
%     if omega > pi/2
%         num_splits = ceil(omega / (pi/2)); % Number of sub-arcs
%         theta_step = omega / num_splits; % Angle per sub-arc
% 
%         % Initialize total hologram
%         Holo_Func = zeros(size(U));
% 
%         % Iterate over sub-arcs and accumulate hologram
%         theta_start = atan2(y0 - yc, x0 - xc);
%         for k = 1:num_splits
%             theta_end = theta_start + theta_step;
% 
%             % Compute new intermediate points along the arc
%             x_mid = xc + cos(theta_end) * (x0 - xc);
%             y_mid = yc + sin(theta_end) * (y0 - yc);
% 
%             % Compute the hologram for this arc segment
%             Holo_Segment = ArcHoloSegment(x0, y0, x_mid, y_mid, xc, yc, U, V, n);
% 
%             % Accumulate into total hologram
%             Holo_Func = Holo_Func + Holo_Segment;
% 
%             % Update for next segment
%             x0 = x_mid;
%             y0 = y_mid;
%             theta_start = theta_end;
%         end
%     else
%         % Compute hologram normally if omega ≤ π/2
%         Holo_Func = ArcHoloSegment(x0, y0, x1, y1, x2, y2, U, V, n);
%     end
% end
% 
% %% Function to compute arc segment hologram
% function Holo_Segment = ArcHoloSegment(x0, y0, x1, y1, x2, y2, U, V, n)
%     % Compute omega for this arc
%     [xc, yc, omega] = ComputeCircleCenter(x0, y0, x1, y1, x2, y2);
% 
%     % Define t values — evenly spaced points for Simpson’s rule
%     t = linspace(0, 1, 2 * n + 1);
% 
%     % Compute α(t) and β(t) for arc trajectory
%     alpha = @(t) (sin((1 - t) * omega) * x0 + sin(t * omega) * x2) / sin(omega);
%     beta = @(t) (sin((1 - t) * omega) * y0 + sin(t * omega) * y2) / sin(omega);
% 
%     % Define the holographic function for integration
%     Holo_Deriv_Func = @(t) exp(2 * pi * 1j * (U .* alpha(t) + V .* beta(t)));
% 
%     % Perform numerical integration using Simpson’s Rule
%     Holo_Segment = Simpsons_Rule(Holo_Deriv_Func, 0, 1, n);
% end
% 
% %% Function to compute circle center and angle omega
% function [xc, yc, omega] = ComputeCircleCenter(x0, y0, x1, y1, x2, y2)
%     % Midpoints
%     mid1_x = (x0 + x1) / 2;
%     mid1_y = (y0 + y1) / 2;
%     mid2_x = (x1 + x2) / 2;
%     mid2_y = (y1 + y2) / 2;
% 
%     % Slopes of perpendicular bisectors
%     dx1 = x1 - x0; dy1 = y1 - y0;
%     dx2 = x2 - x1; dy2 = y2 - y1;
% 
%     if dy1 == 0
%         slope1 = Inf;
%     else
%         slope1 = -dx1 / dy1;
%     end
% 
%     if dy2 == 0
%         slope2 = Inf;
%     else
%         slope2 = -dx2 / dy2;
%     end
% 
%     % Solve for intersection (circle center)
%     if isinf(slope1)
%         xc = mid1_x;
%         yc = slope2 * (xc - mid2_x) + mid2_y;
%     elseif isinf(slope2)
%         xc = mid2_x;
%         yc = slope1 * (xc - mid1_x) + mid1_y;
%     else
%         A = [1, -slope1; 1, -slope2];
%         B = [mid1_y - slope1 * mid1_x; mid2_y - slope2 * mid2_x];
%         center = A \ B;
%         xc = center(1);
%         yc = center(2);
%     end
% 
%     % Compute angles
%     theta0 = atan2(y0 - yc, x0 - xc);
%     theta2 = atan2(y2 - yc, x2 - xc);
% 
%     % Compute angle omega (ensuring positive direction)
%     omega = mod(theta2 - theta0, 2 * pi);
%     if omega > pi
%         omega = 2 * pi - omega;
%     end
% end



function Holo_Func = ArcHolo(x0, y0, x1, y1, x2, y2, U, V, n)
    % Step 1: Compute perpendicular bisector gradients using your formulas
    % Slopes from endpoints to the midpoint
    dx01 = x1 - x0; dy01 = y1 - y0;
    dx12 = x2 - x1; dy12 = y2 - y1;

    % Handle division by zero by using Inf
    if dy01 == 0
        m1_perp = Inf;
    else
        m1_perp = (x0 - x1) / (y0 - y1);
    end

    if dy12 == 0
        m2_perp = Inf;
    else
        m2_perp = (x1 - x2) / (y1 - y2);
    end

    % Y-intercepts using your equations
    if isinf(m1_perp)
        x_c = x0;
        y_c = m2_perp * (x_c - x1) + y1 - m2_perp * x1;
    elseif isinf(m2_perp)
        x_c = x2;
        y_c = m1_perp * (x_c - x0) + y0 - m1_perp * x0;
    else
        % Solve y1(x) = y2(x) to find intersection
        % Line 1: y = m1*x + c1
        c1 = y0 - m1_perp * x0;
        % Line 2: y = m2*x + c2
        c2 = y1 - m2_perp * x1;

        x_c = (c2 - c1) / (m1_perp - m2_perp);
        y_c = m1_perp * x_c + c1;
    end

    % Step 2: Compute angles from (xc, yc) to start and end points
    theta_start = atan2(y0 - y_c, x0 - x_c);
    theta_end   = atan2(y2 - y_c, x2 - x_c);

    % Step 3: Compute total swept angle
    omega = mod(theta_end - theta_start, 2 * pi);
    
    % Correct for direction if the middle point suggests clockwise
    theta_mid = atan2(y1 - y_c, x1 - x_c);
    if ~isBetween(theta_mid, theta_start, theta_end)
        omega = 2 * pi - omega;
    end

    % Step 4: Define trajectory functions
    t = linspace(0, 1, 2 * n + 1); % 2n+1 points for Simpson’s Rule

    radius = sqrt((x0 - x_c)^2 + (y0 - y_c)^2); % constant radius

    alpha = @(t) x_c + radius * cos(theta_start + t * omega);
    beta  = @(t) y_c + radius * sin(theta_start + t * omega);

    % Step 5: Hologram integration
    Holo_Deriv_Func = @(t) exp(2 * pi * 1j * (U .* alpha(t) + V .* beta(t)));
    Holo_Func = Simpsons_Rule(Holo_Deriv_Func, 0, 1, n);

    % Optional: normalize by arc length
    arc_length = radius * omega;
    chord_length = sqrt((x2 - x0)^2 + (y2 - y0)^2);
    if arc_length > 3 * chord_length
        disp(['Clamped arc length from ', num2str(arc_length), ' to ', num2str(3 * chord_length)]);
        arc_length = 3 * chord_length;
    end

    Holo_Func = Holo_Func * arc_length;

    % Debug output
    disp(['→ Max |Holo_Func| = ', num2str(max(abs(Holo_Func(:))))]);
end

% Utility function to check if angle lies between two others (counter-clockwise)
function result = isBetween(theta, theta_start, theta_end)
    theta = mod(theta, 2*pi);
    theta_start = mod(theta_start, 2*pi);
    theta_end = mod(theta_end, 2*pi);

    if theta_start < theta_end
        result = theta_start < theta && theta < theta_end;
    else
        result = theta > theta_start || theta < theta_end;
    end
end
