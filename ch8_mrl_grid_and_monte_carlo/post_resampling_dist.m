% ==============================================================

%compute probability distribution of particle location after resampling
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%N: number of particles 
%p: 1 by 4 multinomial probabilities (k bins). must sum to 1.
%z: 1 by 4 measurement model probabilities p(z | x ).

%P: 1 by 4 by 2 vectors containing distributions P(X | z ) and P( X | not z).

% ==============================================================
function [ P1 ] = post_resampling_dist(N,p,z)

    [X1,X2,X3] = ndgrid(0:N);
    
    X4 = N - (X1+X2+X3);
    
    X1 = X1(:);X2=X2(:);X3=X3(:);X4=X4(:);
    %filter out vectors not in DN 
    ind = X4 >= 0;
    X1 = X1(ind);X2=X2(ind);X3=X3(ind);X4=X4(ind);
    
    Y = mnpdf([X1,X2,X3,X4],repmat(p,size(X4,1),1));
   
    
    Z=zeros(1,4,2);
    Z(1,:,1)= z;
    Z(1,:,2)= 1-z;
    
    X = bsxfun(@times, [X1, X2, X3, X4], Z);
    eta = sum(X,2);
    eta = 1./eta;

    X = bsxfun(@times, bsxfun(@times,Y,eta), X);
    P1 = sum(X,1);
end