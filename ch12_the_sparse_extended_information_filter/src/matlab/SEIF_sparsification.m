% ==============================================================

% SEIF sparsification
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%xi: information vector
%O: information matrix
%m0:  index active landmarks
%m1:  index landmarks to be deactivated
%m2:  index passive landmarks

% =============================================================
function [  ] = SEIF_sparsification( m0,m1,m2)

    global m xi O
    
    global G;
    sg=size(G,1);
    mask=ones(sg,sg)-sparse(ones(1,size(m1,2)),m1+1,ones(1,size(m1,2)),sg,sg);
    G= G & mask;
    [s,t]=meshgrid(union(m1,m0),m1);
    s=s(:)'+1;t=t(:)'+1;
    edges=unique(sort([s;t],1)','rows')';
    edges=edges(:,edges(1,:)-edges(2,:)~=0);
    G=G | sparse(edges(1,:),edges(2,:) ,size(edges,2),sg,sg);

    s0=size(m0,2);s1=size(m1,2);s2=size(m2,2);s=size(m,1);
    m0=sort(m0,2);m1=sort(m1,2);m2=sort(m2,2);
    m0 = 2 + 2*m0;
    m0 = reshape([m0;m0+1],1,[]);
    m1 = 2 + 2*m1;
    m1 = reshape([m1;m1+1],1,[]);
    s0=size(m0,2);s1=size(m1,2);s2=size(m2,2);s=size(m,1);

    O2 = O([1:3 m0 m1],m1);
    L1 = O2/O(m1,m1)*O2';
    O3=O([1:3 m0 m1],[1:3 m1]);
    L2 = O3/O([1:3 m1],[1:3 m1])*O3';
    O4=O(:,1:3);
    L3 = O4/O(1:3,1:3)*O4';

    O1 = O;
    O1([1:3 m0 m1],[1:3 m0 m1]) = O1([1:3 m0 m1],[1:3 m0 m1])-L1+L2;
    %O1([1:3 m0 m1],[1:3 m0 m1]) = O1([1:3 m0 m1],[1:3 m0 m1])+L2;
    O1 = O1 - L3;
    
    %reset to 0 for sparse conversion
    O1(1:3,m1)=zeros(3,size(m1,2));
    O1(m1 ,1:3)=zeros(size(m1,2),3);
    

      xi = xi + (O1-O)*m;
      O=O1;
    
end