function eraser = Rect_Eraser(a,b, x, y)
%RECT_ERASER Summary of this function goes here
%   Detailed explanation goes here
eraser= a*b*sinc(pi*a*x).*sinc(pi*b*y);

end

