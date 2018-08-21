% ==============================================================

% jacobian measurement model
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%m: current mean estimate
%c: correspondence N x 1

%J: jacobian of measurement model 2 x 5 x N

% =============================================================
function [ J ] = jacobian_measurement(m,c)

    N = size(c,2);
    c1 = 2 + 2*c';
    x = m(c1);y = m(c1+1);
    
    delta = bsxfun(@minus ,[x';y'], m(1:2));
    q = sum(delta.^2);
    
    J = [1 0;0 -1;0 1;+1 0]*delta;
    J = [-J;zeros(1,N);-ones(1,N);J];
    J(1:2:end,:,:) = bsxfun(@times, J(1:2:end,:,:), 1./sqrt(q));
    J([2 4 8 10],:,:) = bsxfun(@times, J([2 4 8 10],:,:), 1./(q));
    
    J = reshape(J,2,5,N);
end