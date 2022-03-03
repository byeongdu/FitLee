function [F, Vpy] = saxspyramid(qx0, qy0, qz0, RH)
% function [F, Vpy] = saxspyramid(qx, qy, qz, RH)
% this is for a pyramid but not truncated pyramid...
% for truncated pyramid, another input alpha is required....
% In this pyramid, bottom square is on xy plane and its centr of mass is at
% [0,0,0]. Vertices of the bottom square is [R, R; -R, R; R, -R; -R, -R]
% Lazzari's paper define F(q) = Int: exp(iqr)dr instead of 
%   F(q) = Int: exp(-iqr)dr. Thus, I multiply -1 to the vector q.
% Example
% [F, Vpy] = saxspyramid([1;2], [2;3], [0,1], [5, 10])
% [F, Vpy] = saxspyramid([1;2], [2;3], [0,1], [5, 10; 6, 11;7,22])
%

if size(RH,2) == 1
    Rm = RH(:,1);
    Hm = sqrt(2)*Rm;   % equal edge length pyramid...
else
    Rm = RH(:,1);Hm = RH(:,2);
end
% if size(Rm > 1)
%     [R, qx0] = meshgrid(Rm, qx0);
%     [H, qy0] = meshgrid(Hm, qy0);
%     [~, qz0] = meshgrid(Rm, qz0);
% else
    R = Rm;
    H = Hm;
% end

tanalpha = H./R;
%Vpy = 4/3*tanalpha.*(R.^3-(R-H./tanalpha).^3)
Vpy = 4*R.^2.*H/3;

%F = saxspyramid0(qx0, qy0, qz0, [R, H]);
%return

j = sqrt(-1);
% qx0 = -qx0;%+eps; % don't add eps, then it will not calculate correctly.
% qy0 = -qy0;%+eps; % don't add eps, then it will not calculate correctly.
% qz0 = -qz0;%+eps; % don't add eps, then it will not calculate correctly.

tol = 1E-14;
con1 = abs(qx0) < tol;
con2 = abs(qy0) < tol;
con3 = abs(qz0) < tol;


% tqx = find(abs(qx) < eps);
% qxy0default = 0.00001;
% if ~isempty(tqx)
%     qx(tqx) = qx(tqx)+qxy0default;
% end
% tqy = find(abs(qy) < eps);
% if ~isempty(tqy)
%     qy(tqy) = qy(tqy)+qxy0default;
% end

t{1} = con1 & con2 & con3;
t{2} = con1 & con2 & ~con3; % qx=0 qy=0, qz ~=0
t{3} = con1 & ~con2; % qx=0, qy~=0
t{4} = ~con1 & con2; % qx~=0, qy==0
t{5} = ~(t{1} | t{2} | t{3} | t{4});

F = zeros(size(qx0));

for i=1:numel(t)
    qx = qx0(t{i});
    qy = qy0(t{i});
    qz = qz0(t{i});
    switch i
        case 2
            Ez = exp(-j*qz*H);
            F0 = j*(Ez-1)./qz + 2/H*(1-Ez.*(1+j*H*qz))./qz.^2 + 1/H^2*(Ez.*(H*qz.*(2+j*H*qz)-2*j)+2*j)./qz.^3;
            F(t{i}) = F0*4*R^2;
        case 1
            F(t{i}) = Vpy;
        case 3
            F(t{i}) = qz_RH_qy(R, H, qy, qz);
            
        case 4
            F(t{i}) = qz_RH_qy(R, H, qx, qz);
            
        case 5
            
            q1 = -1/2*((qx-qy)./tanalpha + qz);
            q2 = -1/2*((qx-qy)./tanalpha - qz);
            q3 = -1/2*((qx+qy)./tanalpha + qz);
            q4 = -1/2*((qx+qy)./tanalpha - qz);
            K1 = sinc(q1.*H).*exp(j*q1.*H) + sinc(q2.*H).*exp(-j*q2.*H);
            K2 = -j*sinc(q1.*H).*exp(j*q1.*H) + j*sinc(q2.*H).*exp(-j*q2.*H);
            K3 = sinc(q3.*H).*exp(j*q3.*H) + sinc(q4.*H).*exp(-j*q4.*H);
            K4 = -j*sinc(q3.*H).*exp(j*q3.*H) + j*sinc(q4.*H).*exp(-j*q4.*H);
            F(t{i}) = H./(qx.*qy).*(cos((qy-qx).*R).*K1 + sin((qy-qx).*R).*K2 - cos((qx+qy).*R).*K3 + sin((qx+qy).*R).*K4);
    end
end
end

function F = qz_RH_qy(R, H, qy, qz)
    tol = 1E-14;
    j = sqrt(-1);
    c = 2*R*j./qy.*exp(-j*qz*H);
    q1 = -qy*R/H+qz;
    q2 = qy*R/H+qz;
    F0 = -j*(exp(j*q1*H))./q1...
        + j*(exp(j*q2*H))./q2...
        -(1-exp(j*q1*H))./q1.^2/H...
        + (1-exp(j*q2*H))./q2.^2/H;

    k = abs(q1) < tol;
    if sum(k)>0
        q2n = q2(k);
        F0(k) = H/2 + j*exp(j*q2n*H)./q2n ...
            +(1-exp(j*q2n*H))./q2n.^2/H;
    end
    k = abs(q2) < tol;
    if sum(k)>0
        q1n = q1(k);
        F0(k) = -H/2 - j*exp(j*q1n*H)./q1n ...
            - (1-exp(j*q1n*H))./q1n.^2/H;
    end

    F = c.*F0;
end
% F = F./Vpy;
% t = isinf(F); F(t) = 1;
% t = find(abs(F) > 1); F(t) = 1;
% F = reshape(F, size(qx));
% F = F.*Vpy;