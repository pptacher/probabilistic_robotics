% ==============================================================
%Implementation of histogram filter algorithm.
%written by Pierre-Paul TACHER (ptacher@gmail.com).
% ==============================================================

%initial , mean, covariance matrix
S = [0.01 0;0 0.01];
mu=[0 0];

cell = 0.015;
[xmin ymin] = deal(-1.5, -1.5);

%size of square grid.
n = 200;

prob1 = zeros(n, n);
prob2 = zeros(n, n);

%precompute transition probabilities.
%requires (2n-1)^2 probabilities. needs to be done only once.
[X1,Y1] = meshgrid(0:cell:(n-1)*cell);

tmp = arrayfun(@(x, y) detect_circle(x, y, cell), X1, Y1);
D = zeros(2*n-1, 2*n-1);
tmp = [tmp(n:-1:2,:);tmp];
D = [tmp(:,n:-1:2) tmp];

%generate coordinates of cell centers.
[X,Y]=meshgrid(xmin:cell:xmin+(n-1)*cell);

X = reshape(X,n^2,1);
Y = reshape(Y,n^2,1);
Z= [X Y];

%compute initial distribution.
D1 = pdist2(Z, mu, 'mahalanobis', S);

prob1 = exp(-1/2*reshape(D1.^2,n,n));
eta1 = 1./sum(sum(prob1,1),2);
prob1 = eta1*prob1;

figure;
%plotting takes a while.
plotdistrib3(xmin:cell:xmin+(n-1)*cell,xmin:cell:xmin+(n-1)*cell,prob1, cell);

for ee=1:1
    
    %this should be vectorized probably.or using blockproc-like function.
    for ii=1:n
        for jj=1:n
            prob2(ii,jj) = sum(sum((D(n+1-jj:2*n-jj , n+1-ii:2*n-ii )).*prob1,1),2);
        end
    end
    
    eta = 1./sum(sum(prob2,1),2);
    prob2 = eta.* prob2;

    figure;
    plotdistrib3(xmin:cell:xmin+(n-1)*cell,xmin:cell:xmin+(n-1)*cell,prob2, cell);
    prob1 = prob2;
end

%measurement update.
z = 0.65;
Q = 0.01;
N = exp(-1/(2*Q)*((xmin:cell:xmin+(n-1)*cell)-z).^2);

prob1 = bsxfun(@times, N', prob1);
eta = 1./sum(sum(prob1,1),2);
prob1 = eta.* prob1;

figure;
plotdistrib3(xmin:cell:xmin+(n-1)*cell,xmin:cell:xmin+(n-1)*cell,prob1, cell);
    