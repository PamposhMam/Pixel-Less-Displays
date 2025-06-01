%Square
% Define coordinates for multiple points
% Fixed starting points (bottom-left corner of the square)
a0 = num_cols / 2; % Starting horizontal point
b0 = num_rows / 2; % Starting vertical point

% Define the ending points to form a square
x1 = [a0 + 100, a0 + 100, a0, a0]; % Horizontal points of the square (clockwise or counterclockwise)
y1 = [b0, b0 + 100, b0 + 100, b0]; % Vertical points of the square (clockwise or counterclockwise)


%5-Pointed Star
a0= num_cols/2;
b0= num_ros/2;
x1= [a0 + 6*cos(pi/10), a0 + 2*cos(3*pi/10), a0 + 6*cos(pi/2), a0 + 2*cos(7*pi/10), a0 + 6*cos(9*pi/10), a0 + 2*cos(11*pi/10), a0 + 6*cos(13*pi/10),0, a0 + 6*cos(17*pi/10), a0 + 2*cos(19*pi/10)];
y1= [b0 + 6*sin(pi/10), b0 + 2*sin(3*pi/10), b0 + 6*sin(pi/2), b0 + 2*sin(7*pi/10), b0 + 6*sin(9*pi/10), b0 + 2*sin(11*pi/10), b0 + 6*sin(13*pi/10),b0 -2, b0 + 6*sin(17*pi/10), b0 + 2*sin(19*pi/10)];