function HoloFunct = LineHolo(x0, y0, x1, y1, U, V)
    % Debugging: Print Inputs
    %U=U/2; %Uncomment these two lines if i'm using -1,1 on the meshgrid
    %V=V/2;
    disp(['LineHolo called with x0 = ' num2str(x0) ', y0 = ' num2str(y0) ', x1 = ' num2str(x1) ', y1 = ' num2str(y1)]);

    % Compute L_num1 using complex exponential
    L_num1 = 1j * exp(1j * 2 * pi * (U * x0 + V * y0));

    % Define the denominator term L_denom
    L_denom = 2 * pi * ((x1 - x0) * U + (y1 - y0) * V);

    % Avoid division by zero (set very small denominators to a safe value)
    L_denom(abs(L_denom) < 1e-2) = 1e-2;

    % Compute L_num2
    L_num2 = 1 - exp(1j * 2 * pi * ((x1 - x0) * U + (y1 - y0) * V));

    %I shall include the bit i THINK Jake left out!
    scale= ((x1-x0)^2+(y1-y0)^2)^(0.5);

    % Compute Hologram Function
    HoloFunct = ((L_num1 .* L_num2) ./ L_denom)*scale; %*5 try adjusting the scale???

    % Debugging: Check if result is valid
    if isempty(HoloFunct) || all(isnan(HoloFunct(:))) || all(abs(HoloFunct(:)) < 1e-10)
        warning('LineHolo produced invalid output at x0=%f, y0=%f, x1=%f, y1=%f', x0, y0, x1, y1);
    end
end
