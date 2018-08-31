% ==============================================================

% SEIF motion update
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%xi: initial information vector
%O: initial information matrix
%m: initial mean
%v: speed 1 x 3
%a: steering 1
%t: time step
%m0: active landmarks

% =============================================================
function [ ] = SEIF_motion(v, a, t,m0)
   
    global G;
    
    global m xi O
    
     n = size(xi,1);
    
    R= diag([8e-4 8e-4 1e-6]);
    J = jacobian_motion(m(3,1), v, a, t);
    
    psi = (inv(J)-eye(3));
    
    OX=sparse(O(:, 1:3));
    psi1 = OX*psi;

    lambda=sparse(n,n);
    lambda(:, 1:3) = psi1;lambda(1:3,:) = lambda(1:3,:) + psi1';
    lambda(1:3, 1:3)= lambda(1:3, 1:3) + psi'*O(1:3,1:3)*psi;
    
    phi = O + lambda;
    phiX=sparse(phi(:,1:3));
    k = phiX*inv(inv(R)+ phi(1:3,1:3))*phiX';
    O = phi-k;

    dx = equation_motion(m(3,1), v, a, t);
    OX=sparse(O(:, 1:3));
    xi = xi + lambda*m-k*m+OX*dx;
    m = m + [dx;zeros(n-3,1)];
    
    clear OX
  
    if ~isempty(m0)
        [s,t]=meshgrid(m0);
        s=s(:)'+1;t=t(:)'+1;
        edges=unique(sort([s;t],1)','rows')';
        edges=edges(:,edges(1,:)-edges(2,:)~=0);
        if ~isempty(edges)
             G=G | sparse(edges(1,:),edges(2,:) ,size(edges,2),size(G,1),size(G,1));
        end
    end
    
end