% ==============================================================

% SEIF measurement update
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%xi: information vector
%O: information matrix
%m: current mean estimate
%z: N x 2 measurements
%c: correspondence N x 1

% =============================================================
function [] = SEIF_measurement( z, c)
    global G;global Q;
    global m xi O
    
    n = (size(xi,1)-3)/2;
    new = max(c)-n;

    if new > 0
        O1 = zeros(size(xi,1)+2*new,size(xi,1)+2*new);
        O1(1:size(xi,1),1:size(xi,1))= O;
        xi1 = zeros(size(xi,1)+2*new,1);
        xi1(1:size(xi,1),1) = xi;
        m1 = zeros(size(xi,1)+2*new,1);
        m1(1:size(xi,1),1) = m;
        m1(size(xi,1)+1:end,1) = reshape(inverse_measurement(m(1:3),z(:,c>n)),[2*new,1]);
        s=size(G,1);
        G(s+new,s+new)=0;
    else
        O1 = O;
        xi1 = xi;
        m1=m;
    end

    J = jacobian_measurement(m1,c);
    zhat = equation_measurement(m1,c);

    for ii=1:size(c,2)

        jj = c(ii);
        H = J(:,:,ii);
        dz = z(:,ii)-zhat(:,ii);
        dz(2) = measure(dz(2));
        xi1([1:3 2*jj+2 2*jj+3],1) = xi1([1:3 2*jj+2 2*jj+3],1) + H'*inv(Q)*(dz+H*m1([1:3 2*jj+2 2*jj+3],1));
        O1([1:3 2*jj+2 2*jj+3],[1:3 2*jj+2 2*jj+3]) =  O1([1:3 2*jj+2 2*jj+3],[1:3 2*jj+2 2*jj+3])+ H'*inv(Q)*H; 
        %change inv Q \

    end
   m=m1;
   xi=xi1;
   %O=sparse(O1);
   O=O1;
    %unneeded in theory
    c = unique(c);
    G=G | sparse(ones(1,size(c,2)), c+1,ones(1,size(c,2)),(size(xi,1)-3)/2+1,(size(xi,1)-3)/2+1);
   
   
end

function alpha = measure(theta)
    tmp =mod(theta,2*pi);
    %there is a bug in matlab mod 
    alpha = tmp-min(floor(tmp/pi),1)*2*pi;
end