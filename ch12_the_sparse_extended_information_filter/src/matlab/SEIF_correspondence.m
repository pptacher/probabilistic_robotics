% ==============================================================

% SEIF correspondence determination
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%z: N x 2 measurements

%c: indexes of landmarks corresponding to observations

% =============================================================

function [ c ] = SEIF_correspondence( z,m0)


global G;global Q;

global m xi O
k = size(z,2);
if k==0
    c=zeros(1,0);
    return;
end
c= zeros(1,k);
n = (size(xi,1)-3)/2;

if n > 0
    
    %rough landmarks selection
    t = m(4:end)';t1=t(1:2:end);t2=t(2:2:end);
    lst = find(abs(t1-m(1))<=75 ...
        & abs(t2-m(2))<=75 ...
        & (t1-m(1))*cos(m(3))+(t2-m(2))*sin(m(3))>=-5);
    a=size(lst,2);
    A = zeros(2*a,2);
    
    J = jacobian_measurement(m,lst);
    zhat = reshape(equation_measurement(m,lst),2,1,a);
    
    for ii=1:a
        
        ind = markov_blanket(lst(ii),m0);
        
        ind=ind(2:end);
        c1 = 2 + 2*ind;c2=c1+1;ind=reshape([c1;c2],[1,2*size(ind,2)]);
        
        ind1 = find(ind==2+2*lst(ii));
        
        F=Fxm(ind1,3+size(ind,2));
        S = J(:,:,ii)*F/O([1:3 ind],[1:3 ind]);
        
        [U,M1]=eig(S*F'*J(:,:,ii)'+Q);
        U = diag(sqrt(1./diag(M1)))*U';
        A(2*ii-1:2*ii,:) = U;
        
        clear F;
        
    end
    
    F = reshape(permute(bsxfun(@minus,zhat,z),[1,3,2]),2*a,k);
    F(2:2:end,:)=measure(F(2:2:end,:));
    %first pass over active landmarks
    N = zeros(2*a,k);
    [C,ia,~] = intersect(lst,m0,'stable');
    
    for j=1:size(C,2)
        jj=ia(j);
        N(2*jj-1:2*jj,:)= A(2*jj-1:2*jj,:)*F(2*jj-1:2*jj,:);
    end
    N1 = N(reshape([2*ia'-1;2*ia'],1,[]),:);
    N1 = N1.^2;M2=N1(2:2:end,:);N1 = N1(1:2:end,:)+M2;
    [m1, c] = min(N1,[],1);%search min along first dim, even if row vector.
    
    new = find(m1>chi2inv(0.5,2));
    new1 = find(m1(new)<chi2inv(0.95,2));ia=ia(setdiff(1:size(ia,1),c(new(new1)),'stable'));
    c = C(1,c);
    %second pass over left over landmarks
    if (~isempty(new) && size(C,2)<a) || isempty(C)
        if isempty(C)
            new=1:k;
        end
        ib=setdiff(1:a,ia,'stable');%ib is a row
        lst = lst(1,ib);
        for j=1:size(lst,2)
            jj=ib(j);
            N(2*jj-1:2*jj,new)= A(2*jj-1:2*jj,:)*F(2*jj-1:2*jj,new);
        end
        if ~isempty(ib)
            N1 = N(reshape([2*ib-1;2*ib],1,[]),new);
            N1 = N1.^2;M2=N1(2:2:end,:);N1 = N1(1:2:end,:)+M2;
            [m2, c(new)] = min(N1,[],1);%search min along first dim, even if row vector.
            c(new) = lst(1,c(new));
            new = new(1,find(m2>chi2inv(0.95,2)));
            new1=zeros(1,0);
        end
    end
else
    new = 1:k;
    new1=zeros(1,0);
end
if ~isempty(new1)
    new = new(setdiff(1:size(new,2),new1));
end

if ~isempty(new)
    c(1,new) = n + (1:(size(new,2)));
end

end

function alpha = measure(theta)
tmp =mod(theta,2*pi);
%matlab mod may return 2pi due to precision limitation
alpha = tmp-min(floor(tmp/pi),1)*2*pi;
end

function [ F ] =Fxm(p,n)

if p>0
    F=sparse(1:5, [1 2 3 3+p 4+p],ones(1,5),5,n);
else if n==0
        F=sparse(1:3,1:3,ones(1,3),3,n);
    end
end
end
