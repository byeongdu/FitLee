function [F, V] = saxsoctant0(qx0,qy0,qz0,RLH)
% octant with arbitary dimensions, R, L and H.
% Definition of this octant
% Bottom plane (isoless triangle)
%            ^ y   
%            |   
%            |   
%            | 
%   ----------------------->x
%          . | .
%         .  |  . R
%        .   |   .
%       .....|.....            
%            | L    
%
% H is the height of the octant.
% whose coordinate is at [0,0,H]
% 
% Coordinates of vertices
% vert = [0,0,0;
%       L/2, -sqrt(R^2-(L/2)^2), 0;
%       -L/2, -sqrt(R^2-(L/2)^2), 0;
%       0, 0, H]
% faces = [1,2,3;1,2,4;1,3,4;2,3,4]
%
% Byeongdu Lee
% Dec. 1. 2012
%
% Updates:
% Dec. 3. 2018, Calculation for singularities is added. 

R = RLH(1);
L = RLH(2);
H = RLH(3);

tol = 1E-15;
qx0 = round(qx0/(tol*0.1))*tol*0.1;
qy0 = round(qy0/(tol*0.1))*tol*0.1;
qz0 = round(qz0/(tol*0.1))*tol*0.1;

t = sqrt(R^2-(L/2)^2);
tanth = L/2/t;
costh = t/R;

j = sqrt(-1);

sinth = L/2/R;

C = [0,0,0; 0,0,H; R*sinth, -R*costh, 0; -R*sinth, -R*costh, 0];
V = 1/6*det([C(1,:)-C(2,:);C(1,:)-C(3,:);C(1,:)-C(4,:)]);

B1 = (qy0.^2-qx0.^2*tanth^2);

con1 = abs(qx0) < tol;
con2 = abs(qy0) < tol;
con3 = abs(qz0) < tol;
con4 = abs(B1) < tol;
con5 = sign(qx0)==sign(qy0);
t = [];
t{1} = ~con1 & ~con4;
t{2} = con1 & ~con2;
t{3} = ~con1 & ~con2 & con4 & con5;
t{4} = ~con1 & ~con2 & con4 & ~con5;
t{5} = con1 & con2 & ~con3;
t{6} = con1 & con2 & con3;

F = zeros(size(qx0));
Npnt = 0;

for i=1:6
    qx = qx0(t{i});
    qy = qy0(t{i});
    qz = qz0(t{i});
    Npnt = Npnt + sum(t{i});

    if sum(t{i})<1
        continue
    end
    switch i
        case {1} % qx~=0 and D ~=0
            
            q1 = R/H*(costh*qy-sinth*qx)+qz;
            q2 = R/H*(costh*qy+sinth*qx)+qz;

            D = qx.*(qy.^2-qx.^2*tanth^2);
            Sz = sinc(qz*H/2);
            Ez = exp(-j*qz*H/2);
            F2 = qx.*Sz - ...
                1/2*Ez.*((qy/tanth+qx).*exp(j*q1*H/2).*sinc(q1*H/2) - ...
                (qy/tanth-qx).*exp(j*q2*H/2).*sinc(q2*H/2));
            F(t{i}) = -2*H*tanth*Ez./D.*F2;

        case 2 %qx == 0 and qy ~=0
            q0 = qy*R*costh;
            q1 = q0+H*qz;
            c = 2*tanth./qy.^2;
            F1 = j*(exp(-j*qz*H)-1)./qz;
            F2 = (j*q0-1).*exp(j*q0)./q1*H.*(sin(q1)+j*(cos(q1)-1));
            F3 = -j*q0.*exp(j*q0)./q1.^2*H.*(-1+exp(-j*q1).*(1+j*q1));
            F = -c.*(F1+F2+F3);
        case {3, 4} % D == 0 
            q0 = qx*R*tanth*costh;
            q1 = q0+H*qz;
            c = 2./qx.^2/tanth;
            F1 = (1+j*q0).*j./qz.*(exp(-j*qz*H)-1);
            F2 = -exp(j*q0)*H./q1.*(sin(q1)+j*(cos(q1)-1));
            F3 = -j/H*q0./qz.^2.*(exp(-j*qz*H).*(1+j*qz*H)-1);
            F = c.*(F1+F2+F3);
        case 5 % qx=0, qy = 0, qz~=0
            Ez = exp(-j*qz*H);
            c = R^2*costh^2*tanth;
            F1 = c*j./qz.*(Ez-1);
            F2 = -2*c/H./qz.^2.*(Ez.*(j*qz*H+1)-1);
            F3 = c/H^2./qz.^3.*(Ez.*(H*qz.*(2+j*H*qz)-2*j)+2j);
            F(t{i}) = F1+F2+F3;
        case 6
            F(t{i}) = V;
    end
end