function [F, V] = saxstetrahedron(qx,qy,qz, L)
% saxs form factor of a tetrahedron with the edge length L
% Position: 
%  1. A fact of triangle is on the xy plane.
%  2. The mass center of the bottom triangle is at [0, 0, 0], therefore,
%  the 4th tip is on the z axis.
%  3. A tip of the bottom triangle is on the y axis.
% 
tol = 1E-15;
tol2 = tol*1E-10;
% [F, V] = saxstetrahedron2(qx,qy,qz, L);
% return
% [F, V] = saxstetrahedron(qx,qy,qz, L)
% Derived by B. Lee
% 11/26/2012
%
% Updates:
% 12/03/2018 Formula for singularities are derived.

if numel(L) == 2
    H = L(2);
    L = L(1);
else
    H = sqrt(2/3)*L;    % Lengthes of all edges of this tetrahedron is the same.
end
R = L/2;
tanalp = H/R*sqrt(3);
%V = 1/3*tanalp*(R^3-(R-sqrt(3)*H/tanalp)^3);
% S = 4*R*sqrt(H^2+R^2/3);
V = L(1)^2*H/(4*sqrt(3));
% tri_tet = [
%         1  -1/sqrt(3)  0
%       -1 -1/sqrt(3) 0
%       0 2/sqrt(3)  0 ];
% obj.vertices = [tri_tet*R;0,0, H];
% obj.faces = [1, 2, 3; 1, 2, 4; 2, 3, 4; 3, 1, 4];
%   obj = polyhedron_facecheck(obj);
j = sqrt(-1);

D = tanalp*R/sqrt(3) - H/2;
B1 = (qx.^2-3*qy.^2);
B = (qx.*B1);
%D(abs(D)<eps) = eps;

%td1 = abs(D1) < 1E-8;
%td2 = abs(qx) < 1E-8;
con1 = abs(qx) < tol;
con2 = abs(qy) < tol;
con3 = abs(qz) < tol;
con4 = abs(B1) < tol2;
con5 = sign(qx)==sign(qy);
t{1} = ~con1 & ~con4;
t{2} = con1 & ~con2;
t{3} = ~con1 & ~con2 & con4 & con5;
t{4} = ~con1 & ~con2 & con4 & ~con5;
t{5} = con1 & con2 & ~con3;
t{6} = con1 & con2 & con3;
con = zeros(size(con1));
for i=1:numel(t)
    con = con + t{i}*10^(i-1);
end
F = zeros(size(qx));
T = F;
za = linspace(0, H, 50);
za = diff(za)/2+za(1:end-1);
dz = za(2)-za(1);
Npnt = 0;
% B(abs(B)<eps) = eps;

for i=1:6
    qxn = qx(t{i});
    qyn = qy(t{i});
    qzn = qz(t{i});
%     q1n = q1(t{i});
%     q2n = q2(t{i});
%     q3n = q3(t{i});
    Bn = B(t{i});
    Npnt = Npnt + sum(t{i});

    if sum(t{i})<1
        continue
    end
    switch i
        case {1} % qx~=0 and D ~=0
            q1 = (qyn - sqrt(3)*qxn)/tanalp + qzn;
            q2 = (qyn + sqrt(3)*qxn)/tanalp + qzn;
            q3 = 2*qyn/tanalp - qzn;

            F(t{i}) = sqrt(3)*H./Bn.*exp(-j*qzn*R*tanalp/sqrt(3));
            T(t{i}) = (2*qxn.*sinc(q3*H/2).*exp(-j*q3*D) +...
                -(qxn+sqrt(3)*qyn).*sinc(q1*H/2).*exp(j*q1*D) + ...
                -(qxn-sqrt(3)*qyn).*sinc(q2*H/2).*exp(j*q2*D));
            F(t{i}) = F(t{i}).*T(t{i});
            
        case 2 %qx == 0 and qy ~=0
%             q2a = 2*R/(sqrt(3)*H)*qyn-qzn;
%             q1a = sqrt(3)*R/H*qyn-q2a;
            
            q0 = R*qyn/sqrt(3);
            c = 2./(3*qyn.^2);
            q1 = 1/sqrt(3)*R/H*qyn+qzn;
            q2 = -2/sqrt(3)*R/H*qyn+qzn;

            c2 = sqrt(3);
            F1 = (3*q0+j).*exp(j*q0)./q1.*(exp(-j*H*q1)-1);
            F2 = 3*j*q0/H.*exp(j*q0)./q1.^2.*(-1+exp(-j*H*q1).*(1+j*H*q1));
            F3 = j./q2.*(exp(j*H*q2)-1).*exp(-j*(2*q0+H*q2));
            F0 = (F1+F2+F3);
            
            k = abs(q1) < tol;
            if sum(k)>0
                q0 = q0(k);
                q2 = q2(k);
                F1 = (1-3/2*j*q0).*exp(j*q0)*H;
                F3 = j./q2.*(exp(j*H*q2)-1).*exp(-j*(2*q0+H*q2));
                F0(k) = (F1+F3);
            end

            k = abs(q2) < tol;
            if sum(k)>0
                q0 = q0(k);
                q1 = q1(k);
                F1 = (3*q0+j).*exp(j*q0)./q1.*(exp(-j*H*q1)-1);
                F2 = 3*j*q0/H.*exp(j*q0)./q1.^2.*(-1+exp(-j*H*q1).*(1+j*H*q1));
                F3 = -exp(-2*j*q0)*H;
                F0(k) = (F1+F2+F3);
            end
            F(t{i}) = c*c2.*F0;

%            F(t{i}) = 2./(3*qyn.^2).*(F1b+F1a+F2);
%             for k=1:numel(za)
%                 z = za(k);
%                 Rz = R/H*(H-z);
%                 F(t{i}) = F(t{i}) + 2*((-3*j*Rz*qyn+sqrt(3)).*exp(j*sqrt(3)*Rz*qyn)-sqrt(3)).*exp(-2/sqrt(3)*j*Rz*qyn-j*z*qzn)./(3*qyn.^2)*dz;
%             end
        case {3, 4} % D == 0 
            F(t{i}) = 0;
%             for k=1:numel(za)
%                 z = za(k);
%                 ss1 = -(j*q0+1)*exp(-j*q0/3)*exp(j*q1*z);
%                 ss2 = q0*j/H*exp(-j*q0/3)*z*exp(j*q1*z);
%                 ss3 = exp(j*q0*2/3)*exp(-j*q2*z);
% 
%                 S2 = -c*(ss1+ss2+ss3);
%                 FF = FF+S2;
%             end
            q0 = 2*sqrt(3)*R*qyn;
            c = sqrt(3)./(6*qyn.^2);
            q1 = 2/sqrt(3)*R/H*qyn-qzn;
            q2 = 4/sqrt(3)*R/H*qyn+qzn;

            F1 = (q0-j).*exp(-j*q0/3)./q1.*(1-exp(j*q1*H));
            F2 = q0/H./q1.^2.*exp(-j*q0/3).*(exp(j*H*q1).*(j+H*q1)-j);
            F3 = j./q2.*exp(j*q0*2/3-j*q2*H).*(1-exp(j*q2*H));
            F0 = (F1+F2+F3);
            
            k = abs(q1) < tol;
            if sum(k)>0
                q0 = q0(k);
                q2 = q2(k);
                F1 = (3/2*j*q0 + 1)*H.*exp(-j*q0/3);
                F3 = j./q2.*exp(j*q0*2/3-j*q2*H).*(1-exp(j*q2*H));
                F0(k) = (F1+F3);
            end

            k = abs(q2) < tol;
            if sum(k)>0
                q0 = q0(k);
                q2 = q2(k);
                F1 = (q0-j).*exp(-j*q0/3)./q1.*(1-exp(j*q1*H));
                F2 = q0/H./q1.^2.*exp(-j*q0/3).*(exp(j*H*q1).*(j+H*q1)-j);
                F3 = exp(j*q0*2/3)*H;
                F0(k) = (F1+F2+F3);
            end

            F(t{i}) = -c.*F0;            
        case 5 % qx=0, qy = 0, qz~=0
%             %F(t{i}) = 0;
%             for k=1:numel(za)
%                 z = za(k);
%                 Rz = R/H*(H-z);
%                 F(t{i}) = F(t{i})-sqrt(3)*Rz^2*exp(-j*qzn*z)*dz;
%             end
            F(t{i}) = -2*sqrt(3)*j*R^2/H^2./qzn.^3.*(exp(-j*qzn*H)+qzn*j*H+qzn.^2*H^2/2-1);
        case 6
            F(t{i}) = V;
    end
end
mm = abs(F)>V;
if sum(mm)
    F(mm) = V;
end

% mm = abs(F)>V;
% if sum(mm)
%     k = find(mm>0);
%     for i=1:numel(k)
%         fprintf('qx %0.2e, qy %0.2e, qz %0.2e is chosen to the type %i for singularity evaluation, which is not OK.\n', ...
%             qx(k(i)), qy(k(i)), qz(k(i)), con(k(i)));
%         F(k(i)) = V;
%     end
% end

