function [holo] = RECT_colour(amp, coords,u, v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

B2=max(coords(:,1));
A2=max(coords(:,2));
B1=min(coords(:,1));
A1=min(coords(:,2));
A=A2-A1; %height
B=B2-B1; %width
centre=[A2-A/2,B2-B/2];

holo = 2e-1*amp*A*B .* sinc(pi*A*v) .* sinc(pi*B*u) .* exp(-2j*pi*(centre(1)*v + centre(2)*u));



end