function [F, Ftet, Foct] = saxsstellateoctahedron(qx, qy, qz, RH)
% calculate F of saxsstellateoctahedron by tranforming coordinates 
% instead of transforming particles.
% This uses analytical formulae for calculating F of tetrahdedron..
% see saxsstellateoctahedron3.m for comparison. It still stransforming
% coordinates (as oppose to saxsstellateoctahedron2.m), but calculate F
% numerically using a model.
% Stellate Octahedron is composed of an octahedron and 8 tetrahedra on each
% face of the ocahedron. 4 at z>0 and 4 at z<0.

if numel(RH) == 1
    R = RH;
    H = sqrt(2/3)*RH;    % Lengths of all edges of this tetrahedron is the same.
else
    R = RH(1);H = RH(2);
end

% The following will calculate F of the octahedron which orients so that
% its center of mass at [0,0,0], 4 vertices are on xy plane, and 2 vertices
% are on z axis.
% each of 4 vertices on xy plane points [1,1], [1,-1], [-1,-1],[-1,1]
% direction.
Foct = saxsoctahedron(qx,qy,qz,R);


% The following is to calculate F of a tetrahedron.
% The orientation and location of the tetrahedron: 
%   without any rotatio or translation.
% one of its 4 triangle faces is on xy plane. Center of a mass of the
% triangle is on [0,0,0] and one vertice of the triangle is on y axis.
% As a result, a vertice of the tetrahedron is on +z axis.
sinth = sqrt(2/3);
ang = asin(sinth);
angRot = [ang, ang+pi];
%angRot = [ang];
angIn = pi/2*[0,1,2,3];
if numel(size(qx))==2
    Ftet = zeros(size(qx,1), size(qx,2), 8);
else
    Ftet = zeros(size(qx));
end

j = sqrt(-1);
Qall = [qx(:), qy(:), qz(:)];
for l=1:numel(angRot)
    ang = angRot(l);
    for m=1:numel(angIn)
        angxy = angIn(m);
        
        MRP = [1,0,0;0,cos(ang), sin(ang);0, -sin(ang), cos(ang)]; % for Rotation for Particle.
        MRPxy = [cos(angxy), sin(angxy) 0;-sin(angxy), cos(angxy) 0; 0,0,1];
        M = MRP*MRPxy;
        %MRP = [1,0,0;0,cos(ang), sin(ang);0, -sin(ang), cos(ang)]; % for Rotation for Particle.
        %MRPxy = [cos(angxy), sin(angxy) 0;-sin(angxy), cos(angxy) 0; 0,0,1];
        %M = inv(MRPxy*MRP);
        Q = (M*Qall')';
        q1 = reshape(Q(:,1), size(qx));
        q2 = reshape(Q(:,2), size(qx));
        q3 = reshape(Q(:,3), size(qx));
        %Zt = abs(R/sqrt(3)*sin(ang));
        %Yt = abs((1-cos(ang)/sqrt(3))*R);
        Zt = R/sqrt(3)*sin(ang);
        Yt = -(1-cos(ang)/sqrt(3))*R;
        %fprintf('ang, angxy, Zt and Yt = %0.3f, %0.3f, %0.3f and %0.3f\n',ang, angxy, Zt,Yt)
        switch angxy
            case 0
                Zm = exp(-j*qz*Zt);
                Ym = exp(-j*qy*(-1*abs(Yt)));
            case pi/2
                Zm = exp(-j*qz*Zt);
                Ym = exp(-j*qx*abs(Yt));
            case pi
                Zm = exp(-j*qz*Zt);
                Ym = exp(-j*qy*abs(Yt));
            case pi/2*3
                Zm = exp(-j*qz*Zt);
                Ym = exp(-j*qx*(-1*abs(Yt)));
        end
        
        tt = 4*(l-1) + m;
        F0 = saxstetrahedron2(q1,q2,q3,[2*R,H]);
        %if (l==1) & (m==1)
        %    tet = evalin('base', 'tet');
        %    Vtet = 1.1785e+05;
        %    F0 = saxsPolyhedronAmp(q1,q2,q3,tet)*Vtet;
        %end
        if numel(size(q1)) == 2
            Ftet(:,:,tt) = F0.*Zm.*Ym;
        else
            Ftet = Ftet + F0.*Zm.*Ym;
        end

    end
end

%Vtet = 1.1785e+05;
%Ftet = Ftet*Vtet;
if numel(size(Ftet))==2
    F = Foct + sum(Ftet, 3);
else
    F = Foct + Ftet;
end