% ==============================================================

% jacobian motion model
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%theta: robot orientation
%v: speed 1 x 3
%a: steering 1
%t: time step

%J: jacobian of motion model

% =============================================================
function [ J ] = jacobian_motion(theta, v, a, t)

    %car parameters
    H=0.74;
    L=2.83;
    a1= 0.95+L; b1=0.5;
    
    J = [1 0 -t*(v*sin(theta) + v/L*tan(a)*(a1*cos(theta)-b1*sin(theta)));
         0 1 t*(v*cos(theta)-v/L*tan(a)*(a1*sin(theta)+b1*cos(theta)));
         0 0 1];

end