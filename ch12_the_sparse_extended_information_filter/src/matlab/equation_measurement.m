% ==============================================================

%equations of measurement model
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%m:current mean estimate
%c: correspondance vector

%F: measurement model 2 x n

% =============================================================
function [ F ] = equation_measurement(m,c)

    n=size(c,2);
    c1 = 2 + 2*c';
    x = m(c1);y = m(c1+1);

    delta = bsxfun(@minus ,[x';y'], m(1:2));
    F = [sqrt(sum(delta.^2));zeros(1,n)];
    F(2,:) = atan2(delta(2,:),delta(1,:))-m(3)+pi/2;
    
end