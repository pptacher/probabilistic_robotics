% ==============================================================

%compute combined markov blanket of robot position and landmark
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%n: index of landmark (starts from 1)
%D: pair distances of adjacency information graph 

%mb: list of landmarks in combined markov blanket

% =============================================================
function [ mb ] = markov_blanket(n,m0)
    
    global G;
    
    mb1 = m0;
    mb2 = neighbors(G,n+1)'-1;

    if isempty(intersect(mb1,mb2))
        mb = shortestpath(G,1,n+1);
        if ~isempty(mb)%should never be empty
            mb2 = union( mb2,mb-1);
        end
    end
    mb = union([0 n],union(mb1,mb2),'stable');
    
end