% ==============================================================

%SEIF update state mean
%written by Pierre-Paul TACHER (ptacher@gmail.com)


%m0:  index active landmarks

% =============================================================
function [  ] = SEIF_estimate(m0)

   global xi O m
   n=(size(xi,1)-3)/2;
   
   if n < 100
       m=O\xi;
   else
 
        for i=1:size(m0,2)
           m(2*m0(i)+2:2*m0(i)+3) = (O(2*m0(i)+2:2*m0(i)+3,2*m0(i)+2:2*m0(i)+3))...
               \(xi(2*m0(i)+2:2*m0(i)+3)...
               -O(2*m0(i)+2:2*m0(i)+3,:)*m...
               +O(2*m0(i)+2:2*m0(i)+3,2*m0(i)+2:2*m0(i)+3)*m(2*m0(i)+2:2*m0(i)+3));
        end
        m(1:3) = (O(1:3,1:3))\(xi(1:3)-O(1:3,:)*m+O(1:3,1:3)*m(1:3));

   end
end

function [ F ] =Fm(p,n)
    if p>0
        F=sparse(1:2, [2*p+2 2*p+3],ones(1,2),2,n);
    else if p==0
          F=sparse(1:3, 1:3,ones(1,3),3,n);  
        end
    end
end