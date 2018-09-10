% ==============================================================

%compute combined markov blanket of robot position and landmark
%written by Pierre-Paul TACHER (ptacher@gmail.com)

%n: index of landmark (starts from 1)

%mb: list of landmarks in combined markov blanket

% =============================================================
function [ mb ] = markov_blanket(n,m0)
    
    global G;
    
    mb1 = m0+1;
    mb2 = find(G(n+1,:) | G(:,n+1)');
    
    if isempty(intersect_int(mb1,mb2))
         mb = bfs(n+1);%%%%%%%%%%%%%%%this shoul not be empty
         if ~isempty(mb)
             mb2 = union_int( mb2,mb);
         end
    end

    mb2 = [1 n+1 mb2];
    mb = union_int(mb1,mb2);
    mb=mb-1;
    
end


function [ z ] = union_int(m1,m2)
global xi;

s=(size(xi,1)-3)/2+1;
x=zeros(1,s);
y=zeros(1,s);
x(m1) = 1;
y(m2) = 1;
z = find(x | y);


end

function [ z ] = intersect_int(m1,m2)
global xi;
s=(size(xi,1)-3)/2+1;
x=zeros(1,s);
y=zeros(1,s);
x(m1) = 1;
y(m2) = 1;
z = find(x & y);

end

function [ p ] = bfs(n)

global G;
sg=size(G,1);
d=-ones(1,sg);
q = zeros(1,sg);i=1;j=2;q(1)=1;
found=false;
p=zeros(1,0);
curr=1;
while ~found && i<j
    curr=q(i);
    i=i+1;
    next=find((G(curr,:)|G(:,curr)') & d==-1);
    s=size(next,2);
    q(j:j+s-1)=next;j=j+s;d(next)=curr*ones(1,s);
    found=d(n)==curr;
end
if found
    p  = curr;
    prev=d(curr);
    while prev~=1
        p=[p prev];
        prev=d(prev);
    end
end
end
