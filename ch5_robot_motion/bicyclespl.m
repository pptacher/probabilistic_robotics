% ==============================================================

%sample bicyle motion model
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%X1: current state
%alpha, v: commanded steering and forward velocity
%n: number of samples
%t: sample time interval

%X: n x 3 samples matrix.

% ==============================================================
function [ X ] = bicyclespl(X1, alpha, v, t, n)
    X = zeros(n,3);
    [x, y, theta] = deal(X1(1,1), X1(1,2), X1(1,3));
    
    sigma1 = deg2rad(15);
    sigma2 = 0.2;
    l = 1; 
    
    alphas = max(min(normrnd(deg2rad(alpha), sigma1, [n,1]),pi*80/180),-pi*80/180);
    alphas = tan(alphas);
    vs = normrnd(v, sigma2*v, [n,1]);
    
    X = bsxfun(@plus, X, alphas);
    X  = bsxfun(@times, X, [-sin(theta) cos(theta) 1]);
    X = bsxfun(@plus, X, [cos(theta) sin(theta) 0]);
    X =  X .* bsxfun(@times,vs, t*[1 1 1/l]);
    X = bsxfun(@plus, X, [x, y, theta]);
    
end

