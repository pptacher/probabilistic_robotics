% ==============================================================

%sample robot pose from landmark
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%Xi: coordinates of image 
%Xm= coordinates of landmark
%n: number of samples
%s: variance noise
%r: h/f ratio

%X: n x 3 samples matrix.

% ==============================================================
function [ X ] = robotposepl(Xi, Xm, s, r, n)
    X = zeros(n,3);
    [xi, yi, thetai] = deal(Xi(1,1), Xi(1,2), Xi(1,3));
    [xm, ym, thetam] = deal(Xm(1,1), Xm(1,2), Xm(1,3));
    
    thetas = thetam-thetai-unidrnd(4,n,1)*pi/2;
    C = cos(thetas);S=sin(thetas);
    X(:,1:2) = -r*[xi*C-yi*S xi*S+yi*C];
    X = bsxfun(@plus, X, [xm ym 0]);
    X(:,1:2) = X(:,1:2) + mvnrnd([0 0], [s 0;0 s], n);
    X(:,3) = thetas;

end