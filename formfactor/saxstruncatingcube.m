function [Fq, Vtc] = saxstruncatingcube(qx, qy, qz, Lb)
% function F = saxstruncatedcube(qx, qy, qz, Lb)
% L is the edge length of the original cube.
% b is the edge length of the truncated tetrahedron.
if numel(Lb) == 1
    L = Lb;
    b = L/(1+sqrt(2));% ideal truncatedcube.
else
    b = Lb(:,2);
    L = Lb(:,1);
end
Q0 = [qx(:), qy(:), qz(:)];

Fq = zeros(length(qx), numel(L));
F = Fq;
La = L;
ba = b;
[Fcube, Vc] = saxscube(Q0(:,1),Q0(:,2),Q0(:,3),L/2);
for iL=1:numel(La)
    L = La(iL);
    b = ba(iL);
    if b~=0
        [F(:,iL), Vt] = compute(L, b, Q0);
    else 
        F(:,iL) = 0;
        Vt = 0;
    end
    Fq(:,iL) = Fcube - F(:, iL);
end
%Vtc = L^3-sqrt(2)*b^3/3;
Vtc = Vc-Vt;
k = abs(Fq) > Vtc;
if sum(k)>1
    Fq(k) = Vtc;
end

%Ft = Fcube - F;
%Fq = reshape(Ft, size(qx));

function [Fspine, Vt] = compute(L, b, Q0)
Fspine = zeros(size(Q0, 1), 1);
a = polyhedra('truncatedcube', L, L*0.5);
if b == 0
    Rmat = {};
else
    Nf = cellfun(@numel, a.faces);
    Ni = find(Nf == 3);
    Htc = sqrt(3*(L/2-b/sqrt(2)/3)^2);
    % Vicosceless(3, b, Htc)
%     Vicos = [b*sqrt(2)/2, -b*sqrt(2)/2/tan(pi/3), -Htc;
%         -b*sqrt(2)/2, -b*sqrt(2)/2/tan(pi/3), -Htc;
%         0, b*sqrt(2)/tan(pi/3), -Htc];

    tri_tet = [
        1  -1/sqrt(3)  0
      -1 -1/sqrt(3) 0
      0 2/sqrt(3)  0 ];
    Vicos = tri_tet*b/2;
    Vicos(:,3) = -Htc;

    VicosT = inv(Vicos);

    Ht = sqrt(3)*L/2-Htc;
    j = sqrt(-1);
    Vt = 0;

    Rmat = {};
    Rmatinv = {};
    Rmatfaces = {};
    Rmatvertices = {};

    
    for i=1:numel(Ni)
        Fi = a.faces{Ni(i)};
        Mi = round(a.vertices(Fi', :)*100000)/100000;
        
        M = (VicosT*Mi)';
        invM = inv(M);
        
        %invM = VicosT*inv(Mi)';
        if isempty(Rmat)
            Rmat{1} = M;
            Rmatinv{1} = invM;
            Rmatvertices{1} = sort(Mi);
            Rmatfaces{1} = Fi;
        else
            newrot = 1;
            for mati = 1:numel(Rmat)
                if isequal(Rmatvertices{mati}, sort(-Mi))
                    newrot = 0;
                end
            end
            if newrot
                n = numel(Rmat)+1;
                Rmat{n} = M;
                Rmatinv{n} = invM;
                Rmatvertices{n} = sort(Mi);
                Rmatfaces{n} = Fi;
            end
        end
    
    end
end

for i=1:numel(Rmat)
    R = Rmat{i};
    Q = Q0*R;
    [F0,V] = saxstetrahedron(Q(:,1),Q(:,2),-Q(:,3), [b, Ht]);
    Vt = Vt+2*V;
    Fspine = Fspine + 2*real(F0.*exp(j*Q(:,3)*Htc));
%    Fspine = Fspine + F0.*exp(j*Q(:,3)*Htc);
end
