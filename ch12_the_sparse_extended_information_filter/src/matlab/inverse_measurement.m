% ==============================================================

%inverse of measurement function
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%m:current robot position estimate
%c: vector of sensings new landmarks

%X: initialization of new landmarks positions

% =============================================================
function [ X ] = inverse_measurement(m,z)

    C = cos(m(3,1)+z(2,:));
    S = sin(m(3,1)+z(2,:));
    
    X = bsxfun(@plus,m(1:2,1), [z(1,:).*S;-z(1,:).*C]);
    
end