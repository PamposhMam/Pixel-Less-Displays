function [remove] = Angle_Eraser(theta1,theta0, y0, y1, x0, x1)
%ANGLE_ERASER: bec it would be MADNESS to pop that if statement into every
%section of code!!!!
%   Detailed explanation goes here
if (0 <= theta0) && (theta0 < pi/2)
    if (0 <= theta1) && (theta1 < pi/2) && theta1<theta0        %case 1
        a=x1-x0;
        b=y0-y1;
        remove=Rect_Eraser(a,b,x,y);
        remove=remove*exp(-1j*2*pi*((x1+x0)*x/2+(y1+y0)*y/2));
        %bessel_aperture=bessel_aperture-remove;

    elseif (0 <= theta1) && (theta1 < pi/2) && theta1>theta0         %case 2
        a=x1-x0;
        b=y0-y1;
        remove=Rect_Eraser(a,b,x,y);
        remove=remove*exp(-1j*2*pi*((x1+x0)*x/2+(y1+y0)*y/2));

    elseif (pi/2 >= theta1) && (theta1 < pi) && theta1<theta0   %case 3
        % Code for this case

     elseif (pi/2 >= theta1) && (theta1 < pi) && theta1>theta0   %case 4
        % Code for this case

    elseif (pi >= theta1) && (theta1 < 3*pi/2)
        
        % Code for this case
    elseif (3*pi/2 >= theta1) && (theta1 < 2*pi)
        
        % Code for this case
    end
elseif (pi/2 >= theta0) && (theta0 < pi)
    if (0 >= theta1) && (theta1 < pi/2)
        
        % Code for this case
    elseif (pi/2 >= theta1) && (theta1 < pi)
        
        % Code for this case
    elseif (pi >= theta1) && (theta1 < 3*pi/2)
        
        % Code for this case
    elseif (3*pi/2 >= theta1) && (theta1 < 2*pi)
        
        % Code for this case
    end
elseif (pi >= theta0) && (theta0 < 3*pi/2)
    if (0 >= theta1) && (theta1 < pi/2)
        
        % Code for this case
    elseif (pi/2 >= theta1) && (theta1 < pi)
        
        % Code for this case
    elseif (pi >= theta1) && (theta1 < 3*pi/2)
        
        % Code for this case
    elseif (3*pi/2 >= theta1) && (theta1 < 2*pi)
        
        % Code for this case
    end
elseif (3*pi/2 >= theta0) && (theta0 < 2*pi)
    if (0 >= theta1) && (theta1 < pi/2)
        % Code for this case
    elseif (pi/2 >= theta1) && (theta1 < pi)
        % Code for this case
    elseif (pi >= theta1) && (theta1 < 3*pi/2)
        % Code for this case
    elseif (3*pi/2 >= theta1) && (theta1 < 2*pi)
        % Code for this case
    end
end
end

