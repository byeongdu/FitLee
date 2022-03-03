function [F, V, Rmat, Rmatfaces, Rmatvertices, Vbase] = saxsdodecahedron( qx, qy, qz, L)
sizeqx = size(qx);
Qall = [qx(:), qy(:), qz(:)];
F = zeros(size(qx));


%% calculate Rotmatrix.
V = 10*(3+sqrt(5))/3*(L/2)^3;
% obj = polyhedra('dodecahedron');
% [Rmat, Rmatinv, Rmatfaces, Rmatvertices, Vbase, H] = ...
%     getRotMat_Subunit2Face(obj.vertices, obj.faces);
% H = H*L;

% use already computed Rmat. 
p = (1+sqrt(5))/2;
H = sqrt((2+2*p+1/p)^2+(2+p)^2)/5;
t = tan(pi/5);
edgeL = 2/p;
Rmat{1} = [[0,1,0];[-((p^2+2*p)*t)/5,0,-((3*p^2-2*p-1)*t)/5];[(p+2)/(5*H),0,-(2*p^2+2*p+1)/(5*p*H)]]';
Rmat{2} = [[0,1,0];[((p^2+2*p)*t)/5,0,-((3*p^2-2*p-1)*t)/5];[-(p+2)/(5*H),0,-(2*p^2+2*p+1)/(5*p*H)]]';
Rmat{3} = [[p/2,(p-1)/2,(p^2-p)/2];[(p*t)/2,-((4*p^2-p-3)*t)/10,-((3*p^2+p)*t)/10];[0,(2*p^2+2*p+1)/(5*p*H),-(p+2)/(5*H)]]';
Rmat{4} = [[(p-1)/2,(p^2-p)/2,p/2];[((4*p^2-p-3)*t)/10,((3*p^2+p)*t)/10,-(p*t)/2];[-(2*p^2+2*p+1)/(5*p*H),(p+2)/(5*H),0]]';
Rmat{5} = [[0,0,1];[-((3*p^2-2*p-1)*t)/5,((p^2+2*p)*t)/5,0];[-(2*p^2+2*p+1)/(5*p*H),-(p+2)/(5*H),0]]';
Rmat{6} = [[(p-1)/2,-(p^2-p)/2,p/2];[-((p+1)*t)/2,-((p^2+p-2)*t)/10,((2*p^2-p)*t)/10];[0,-(2*p^2+2*p+1)/(5*p*H),-(p+2)/(5*H)]]';
H = H/edgeL*L;
%% Calculate formfactor.
j = sqrt(-1);
for i=1:numel(Rmat)
    M = Rmat{i};
    Q = Qall*M;
    F0 = saxspolygonalcon(Q(:,1), Q(:,2), Q(:,3), [5, L, H]);
    F0 = F0.*exp(j*Q(:,3)*H);
    F = F+2*real(F0);
end
F = reshape(F, sizeqx);
end