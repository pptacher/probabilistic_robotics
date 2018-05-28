% ==============================================================

%bivariate Epanechnikov kernel density estimate
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%X: data whose underlying density is to be estimated. 
%xi: defines grid where the density is to be evaluated. 
%h: window width

%f: values of density estimate.

% ==============================================================

function [ f ] = epankde(X,xi,h)
    m = size(xi, 2);
    n = size(X, 1);
    f = zeros(m,m);
    
    [x1,x2] = meshgrid(xi);

    [B,I] = sort(X(:,1));
    X = [B X(I,2)];
    
    for ii=1:m
        for jj=1:m
        Y = bsxfun(@minus, [x1(ii,jj) x2(ii,jj)],X);
        first = find(Y(:,1)>=-h,1,'first');
        last = find (Y(:,1)<=h, 1, 'last');
        if isempty(first) || isempty(last)
            f(ii,jj)=0;
        else
            Y = Y(first(1):last(1), :);
            N = (1-1/h^2.*dot(Y,Y,2));
            N = N(N >=0);
            f(ii,jj) = 2/(pi*h^2*n)*sum(N,1);
        end
    end
end

