% ==============================================================

%plot 3D histogram representing joint probability distibution 
%inspired from references:

%scatterbar3
%<https://fr.mathworks.com/matlabcentral/fileexchange/1420-scatterbar3>

%Objects creation performance
%<https://blogs.mathworks.com/graphics/2015/06/09/object-creation-performance/>

%Note 'patch' function is called only once.

% ==============================================================

function plotdistrib3(X, Y, Z, width)

i=0;
[r,c]=size(Z);
verts=zeros(8*r,3);
col=zeros(8*r,1);
faces=zeros(6*r,4);
for j=1:r,
    x = X(j);
    for k=1:c,
        if ~isnan(Z(j,k))
            y=Y(k);z=Z(j,k);
            faces(1+6*i:6*i+6,:) = [1 2 3 4;2 5 6 3;5 7 8 6;7 1 4 8;1 2 5 7;4 3 6 8]+8*i;
            w = width/2;
            verts(1+8*i:8*i+8,:) = [x-w y-w 0; x+w y-w 0;x+w y-w z; x-w y-w z;x+w y+w 0;x+w y+w z;x-w y+w 0;x-w y+w z];
            col(1+8*i:8*i+8,:)= [0;0;z;z;0;z;0;z];
            i=i+1;
        end
    end
end

patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData',col,'FaceColor','interp', 'LineWidth', 0.1);
zlim=[min(Z(:)) min(max(Z(:))*1.5,1)];
axis([min(X(:))-width max(X(:))+width min(Y(:))-width max(Y(:))+width zlim])
caxis([min(Z(:)) max(Z(:))])
xlabel('x'); ylabel('y'); zlabel('probability density');
view(15,40);
drawnow;