function integral = Simpsons_Rule(f, a, b, n)
    % Ensure n is a positive integer
    if n <= 0 || mod(n, 1) ~= 0
        error('n must be a positive integer');
    end

    % Use 2n intervals (ensuring even number of intervals)
    num_intervals = 2 * n;
    h = (b - a) / num_intervals;

    % Generate t values
    t = linspace(a, b, num_intervals + 1);

    % Evaluate the function at each t value
    y = arrayfun(f, t, 'UniformOutput', false);

    % Convert cell array to 3D numeric array if function outputs matrices
    y = cat(3, y{:});

    % Apply Simpson's rule across the third dimension (t axis)
    integral = h / 3 * (y(:, :, 1) + y(:, :, end) + 4 * sum(y(:, :, 2:2:end-1), 3) + 2 * sum(y(:, :, 3:2:end-2), 3));
end
