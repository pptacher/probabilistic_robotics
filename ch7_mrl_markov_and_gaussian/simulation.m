% ==============================================================

%simulate motion and sensing of the simplistic underwater robot
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%x0: start position
%t: time step
%v: speed commands n x 3

%X: n x 3 positions
%Z: n x k beacons sensing

% ==============================================================
function [ X, Z ] = simulation(x0, t, v)
    n = size(v,1);
    B = [-10 -10 -10;10 10 -10;2 3 5;2 -4 8];%beacons

    X = zeros(n+1,3);
    X(1,:) = x0;
    k = size(B,1);
    Z = zeros(n+1,4);
    s=0.01;%motion noise variance
    
    for ii=2:n+1
        X(ii,:) = X(ii-1,:) + v(ii-1,:)*t + mvnrnd([0 0 0], s*eye(3));
        Z(ii,:) = sqrt(sum((bsxfun(@minus, X(ii,:), B)).^2, 2))';
    end
    
    X = X(2:n+1,:);
    Z = Z(2:n+1,:);
end