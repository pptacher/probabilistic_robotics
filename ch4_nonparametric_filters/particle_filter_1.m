% ==============================================================
%Implementation of particle filter algorithm.
%written by Pierre-Paul TACHER (ptacher@gmail.com).
% ==============================================================

%motion model matrix.
A = [1 1;0 1];

%size of square grid.
n = 200;

%noise covariance. tweaked to avoid singularity.
R = [1/4 1/2;1/2 1.001];
mu=[ 0 0];

cell = 0.2;
[xmin ymin] = deal(-20, -20);

[x1,x2] = meshgrid(xmin:cell:xmin+(n-1)*cell);
f = zeros(n, n);

%number of particles
M = 10000;

%particles
X = zeros(M, 2);

%weights (importance factors)
w = ones(M, 1);

for ii=1:5
    %sampling from conditional probabilities N(A*xt-1, R)
    X = X*A' + mvnrnd(mu, R, M);
    
     figure;
    axis([-20,20,-15,15]);
    hold on;
    set(gca,'XLimMode','manual','YLimMode','manual','DataAspectRatio',[1,1,1]);
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'Box', 'off');
    plot(gca, X(:,1), X(:,2), '.','MarkerSize', 2);
    hold off;
    
    %Epanechnikov kernel density estimation
    %f = epankde(X, xmin:cell:xmin+(n-1)*cell, 3);
    
    %figure;
    %surf(gca, x1, x2, f);
   
    drawnow;

end

%low variance resampling.
z = 5;
Q = 10;
    
N = exp(-1/(2*Q)*((X(:,1)-z).^2));
w= w.*N;

bins = cumsum(w'*1/sum(w,1));
[~,~,index] = histcounts(1/M*rand + (0:1/M:1-1/M), bins);
X = X(index+1,:);
w = w(index+1);
    
figure;
axis([-20,20,-15,15]);
hold on;
set(gca,'XLimMode','manual','YLimMode','manual','DataAspectRatio',[1,1,1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'Box', 'off');
plot(gca, X(:,1), X(:,2), '.', 'MarkerSize', 2);
hold off;
    
f = epankde(X, xmin:cell:xmin+(n-1)*cell, 3);
    
figure;
surf(gca, x1, x2, f);
   
    