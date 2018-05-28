% ==============================================================
%Implementation of particle filter algorithm.
%written by Pierre-Paul TACHER (ptacher@gmail.com).
% ==============================================================

cell = 0.015;
[xmin ymin] = deal(-1.5, -1.5);

%size of square grid.
n = 200;

[x1,x2] = meshgrid(xmin:cell:xmin+(n-1)*cell);
f = zeros(n, n);

%initial , mean, covariance matrix
S = [0.01 0;0 0.01];
mu=[0 0];

%variance of theta.
s1 = 10000;

%number of particles
M = 10000;

%particles
X = zeros(M, 2);

%weights (importance factors)
w = ones(M, 1);

%sampling from conditional probabilities p(xt | xt-1)
X =  simmov(S,s1,M);

figure;
axis([-1.5,1.5,-1.5,1.5]);
hold on;
set(gca,'XLimMode','manual','YLimMode','manual','DataAspectRatio',[1,1,1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'Box', 'off');
plot(gca, X(:,1), X(:,2), '.','MarkerSize', 2);
hold off;

%Epanechnikov kernel density estimation
f = epankde(X, xmin:cell:xmin+(n-1)*cell, 0.3);

figure;
axis([-1.5,1.5,-1.5,1.5,0 inf]);
hold on;
set(gca,'XLimMode','manual','YLimMode','manual','ZLimMode','auto','DataAspectRatio',[1,1,1]);
surf(gca, x1, x2, f);
drawnow;

%low variance resampling.
z = 0.65;
Q = 0.01;
    
N = exp(-1/(2*Q)*((X(:,1)-z).^2));
w= w.*N;

bins = cumsum(w'*1/sum(w,1));
[~,~,index] = histcounts(1/M*rand + (0:1/M:1-1/M), bins);
X = X(index+1,:);
w = w(index+1);
    
figure;
axis([-1.5,1.5,-1.5,1.5]);
hold on;
set(gca,'XLimMode','manual','YLimMode','manual','DataAspectRatio',[1,1,1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'Box', 'off');
plot(gca, X(:,1), X(:,2), '.', 'MarkerSize', 2);
hold off;
    
f = epankde(X, xmin:cell:xmin+(n-1)*cell, 0.3);
    
figure;
axis([-1.5,1.5,-1.5,1.5,0 inf]);
hold on;
set(gca,'XLimMode','manual','YLimMode','manual','ZLimMode','auto','DataAspectRatio',[1,1,1]);
%surf(gca, x1, x2, f);
   
    