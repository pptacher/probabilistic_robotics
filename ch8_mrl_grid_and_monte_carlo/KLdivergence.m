% ==============================================================

%compute Kullback?Leibler divergence from
%       discrete probability distribution q to p. 
%written by Pierre-Paul TACHER (ptacher@gmail.com)

% ==============================================================
function [ KL ] = KLdivergence(p,q)

    KL = sum(p.*log2(p./q),2);

end