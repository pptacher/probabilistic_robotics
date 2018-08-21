% ==============================================================

%SEIF update state estimate
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%xi: information vector
%O: information matrix
%m0:  index active landmarks

% =============================================================
function [  ] = SEIF_estimate2()

    global m xi O
    m=O\xi;
    
end