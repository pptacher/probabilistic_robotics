% ==============================================================
%Implementation of histogram filter algorithm.
%written by Pierre-Paul TACHER (ptacher@gmail.com).
% ==============================================================

%motion model matrix.
A = [1 1;0 1];
 
%initial , mean, covariance matrix
S = [2.5 2;2 2];
mu=[0 0];
 
%noise covariance. tweaked to avoid singularity.
R = [1/4 1/2;1/2 1.001];
R1 = inv(R);
d = det(R);

cell = 0.2;
[xmin ymin] = deal(-20, -20);

%size of square grid.
n = 200;

prob1 = zeros(n, n);
prob2 = zeros(n, n);

%precompute transition probabilities for gaussian noise.
%requires (2n-1)^2 probabilities. needs to be done only once.
[X1,Y1] = meshgrid(-(n-1)*cell:cell:(n-1)*cell);
X1 = reshape(X1,(2*n-1)^2,1);
Y1 = reshape(Y1,(2*n-1)^2,1);
Z1= [X1 Y1];

[P,D]=eig(R1);
U = sqrt(D)*P';
DD=sum((Z1*U').^2, 2);
D = zeros(2*n-1, 2*n-1);
D = exp(-1/2*(reshape(DD,2*n-1,2*n-1)));

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

for ee=1:3
    
    %motion step

    indexmap = @(l) max(min(l-floor(l/n)+1-ymin/cell,n^2), 1);
    
    %move probability masses according to X -> AX linear map.
    prob1 = prob1(indexmap(reshape(1:n^2, n, n,1)'));

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
z = 5;
Q = 10;
N = exp(-1/(2*Q)*((xmin:cell:xmin+(n-1)*cell)-z).^2);

prob1 = bsxfun(@times, N', prob1);
eta = 1./sum(sum(prob1,1),2);
prob1 = eta.* prob1;

figure;
plotdistrib3(xmin:cell:xmin+(n-1)*cell,xmin:cell:xmin+(n-1)*cell,prob1, cell);
    