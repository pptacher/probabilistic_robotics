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

    lambda(:, 1:3) = psi1;lambda(n,n)=0;lambda(1:3,:) = lambda(1:3,:) + psi1';
    lambda(1:3, 1:3)= lambda(1:3, 1:3) + psi'*O(1:3,1:3)*psi;
    
    phi = O + lambda;
    k = phi(:,1:3)*inv(inv(R)+ phi(1:3,1:3))*phi(1:3, :);
    O = phi-k;

    dx = equation_motion(m(3,1), v, a, t);
    xi = xi + (lambda-k)*m+O(:,1:3)*dx;
    m = m + [dx;zeros(n-3,1)];
    
    clear OX
  
    if ~isempty(m0)
        [s,t]=meshgrid(m0);
        s=s(:)'+1;t=t(:)'+1;
        edges=unique(sort([s;t],1)','rows')';
        edges=edges(:,edges(1,:)-edges(2,:)~=0);
        if ~isempty(edges)
            G=rmedge(G,edges(1,:),edges(2,:));
            G=addedge(G,edges(1,:),edges(2,:));
        end
    end
    
end