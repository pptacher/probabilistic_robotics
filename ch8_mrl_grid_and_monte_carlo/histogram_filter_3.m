% ==============================================================
%Implementation of histogram filter algorithm.
%written by Pierre-Paul TACHER (ptacher@gmail.com).
% ==============================================================

%motion model matrix.
A = eye(3);

deltat = 0.1;%time step
t=0;
B = [-10 -10 -10;10 10 -10;2 3 5;2 -4 8];%beacons

epsilon = 1e-7;
 
%initial , mean, covariance matrix
S = 0.001*eye(3);
mu=[0 0 0];
 
%noise covariance. tweaked to avoid singularity.
R = 0.1*eye(3); 
R1 = inv(R);
d = det(R);

cell = 0.025;
[xmin ymin] = deal(-0.25, -0.25);

%size of square grid.
n = 100;

prob1 = zeros(n, n, n);
prob2 = zeros(n, n, n);
prob3 = zeros(n,n);

%precompute transition probabilities for gaussian noise.
%requires (2n-1)^3 probabilities. needs to be done only once.
[X1,Y1,Z1] = ndgrid(-(n-1)*cell:cell:(n-1)*cell);
X1 = reshape(X1,(2*n-1)^3,1);
Y1 = reshape(Y1,(2*n-1)^3,1);
Z1 = reshape(Z1,(2*n-1)^3,1);
T1= [X1 Y1 Z1];

[P,D]=eig(R1);
U = sqrt(D)*P';
DD=sum((T1*U').^2, 2);
D = zeros(2*n-1, 2*n-1, 2*n-1);
D = exp(-1/2*(reshape(DD,2*n-1,2*n-1, 2*n-1)));

%generate coordinates of cell centers.
[X,Y,Z]=ndgrid(xmin:cell:xmin+(n-1)*cell);

X = reshape(X,n^3,1);
Y = reshape(Y,n^3,1);
Z = reshape(Z,n^3,1);
T= [X Y Z];

%compute initial distribution.
D1 = pdist2(T, mu, 'mahalanobis', S);

prob1 = exp(-1/2*reshape(D1.^2,n,n,n));
eta1 = 1./sum(sum(sum(prob1,1),2), 3);
prob1 = eta1*prob1;

%figure;
%plotting takes a while.
%plotdistrib3(xmin:cell:xmin+(n-1)*cell,xmin:cell:xmin+(n-1)*cell,prob1, cell);

for ee=1:100
    %motion step
    t = t + deltat;
    v = [cos(0.6*t) sin(0.6*t) 0.1];
    indexmap = @(l) max(min(l-sum(floor(deltat/cell.*v).*[1 n n^2],2),n^3), 1);
    
    %move probability masses according to X -> AX + v*deltat affine map.
    prob1 = prob1(indexmap(reshape(1:n^3, n, n, n)));

    eta = 0;
    %this should be vectorized probably.or using blockproc-like function.
    for ii=1:n
        for jj=1:n
            for kk=1:n
                if prob1(ii,jj,kk) >= epsilon
                    prob2(ii,jj,kk) = sum(sum(sum((D(n+1-ii:2*n-ii , n+1-jj:2*n-jj , n+1-kk:2*n-kk)).*prob1)));
                else
                    prob2(ii,jj,kk)= prob1(ii,jj,kk);
                end
            end
        end
    end

    eta = sum(sum(sum(prob2)));
    prob2 = 1/eta.* prob2;

    if mod(ee,10)==0 || ee==1
        figure;
        prob3 = sum(prob2,3);
        plotdistrib3(xmin:cell:xmin+(n-1)*cell,xmin:cell:xmin+(n-1)*cell,prob3, cell);
    end
    prob1 = prob2;
end

%measurement update.
%todo

%figure;
%plotdistrib3(xmin:cell:xmin+(n-1)*cell,xmin:cell:xmin+(n-1)*cell,prob1, cell);
    