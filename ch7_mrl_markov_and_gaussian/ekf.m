% ==============================================================

% E Kalman filter algorithm localization implementation
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%mu: initial mean
%S: initial covariance
%v: commands n x 3
%z: measurements n x k
%t: time step

%X: n x 3 means
%Sigma: 3 x 3 x n covariance matrices

% ==============================================================
function [ X, Sigma ] = EKF(mu, S, v, z, t)
    n = size(v,1);%number of steps

    X = zeros(n+1,3);
    X(1,:) = mu;
    k = size(z,2);
    Sigma = zeros(3,3,n+1);
    Sigma(:,:,1) = S;
    s=0.01*eye(3);%motion noise variance
    K = zeros(3,1);
    B = [-10 -10 -10;10 10 -10;2 3 5;2 -4 8];%beacons
    
    V = eye(3);
    G = eye(3);
    Q = 0.001;
    
    H = zeros(k,3);
    
    for ii=2:n+1
        %motion update
        X(ii,:) = X(ii-1,:) + v(ii-1,:)*t;
        Sigma(:,:,ii) = Sigma(:,:,ii-1) + s;
        
        %measurement update
        H = bsxfun(@minus, X(ii,:), B);
        q = sqrt(sum(H.^2,2));
        H = bsxfun(@times, 1./q, H); 
        for jj = 1:k
            K = Sigma(:,:,ii) *H(jj,:)'* inv(H(jj,:)*Sigma(:,:,ii)*H(jj,:)'+Q);
            X(ii,:) = X(ii,:) + K'*(z(ii-1,jj)-q(jj));
            Sigma(:,:,ii) = (eye(3) - K*H(jj,:))*Sigma(:,:,ii);
        end
    end
    
    X = X(2:n+1,:);
    Sigma = Sigma(:,:,2:n+1);
end