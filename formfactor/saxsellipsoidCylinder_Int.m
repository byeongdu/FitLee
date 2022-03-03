function [f, name, pnames, pin] = saxsellipsoidCylinder_Int(q, p, flag)
% cylinder with ellipsoid cross section.
% eccentricity axis is qx direction.
%R = [Rin2, Rin1, Rout];
%rho = [rho1, rho2, rho3, rho4];
% qp and qz should have a same size....

if nargin == 2
    
    qpORqz = abs(p(1));
    Rp1 = abs(p(2));
    Rp2 = abs(p(3));
    Rp3 = abs(p(4));
    Rz1 = abs(p(5));
    Rz2 = abs(p(6));
    Rz3 = abs(p(7));
    rho1 = p(8);
    rho2 = p(9);
    rho3 = p(10);
    rho4 = 0;
    I0 = p(11);
    dR = abs(p(12));

    if (dR ~= 0)
        [gDistR, Rp3out] = schultzdist99(Rp3, dR, 10);%gDistR=gDistR*(Rout(2)-Rout(1));
    else
        gDistR = 1;
        Rp3out = Rp3;
    end
    
    if qpORqz == 1   % qp
        qp = q;qz = zeros(size(qp));
    elseif qpORqz == 0  %qz
        qz = q;qp = zeros(size(qz));
    elseif qpORqz == 2 % 2D data.
        [qp, qz] = meshgrid(q, q);
    end
%    onex = ones(size(qp));
%    zerox = zeros(size(qz));
    alpha = acos(qp./sqrt(qp.^2 + qz.^2+eps));
    q = sqrt(qp.^2 + qz.^2);
    rho = [rho1, rho2, rho3, rho4];
    Vall = 0;
    f = zeros(size(qp));
    
    for pd=1:length(gDistR)
%f = besseljc(q.*Re);
        Rp = [Rp1, Rp2, Rp3];
        Rz = [Rz1, Rz2, Rz3];

        V = 0;
        Rp = Rp/Rp3*Rp3out(pd);
        Rz = Rz/Rp3*Rp3out(pd);
        
        fq = zeros(size(qp));
        for i=1:length(Rp)
            Re = sqrt(Rz(i)^2*sin(alpha).^2 + Rp(i)^2.*cos(alpha).^2);
            fq = fq + 2*pi*(rho(i) - rho(i+1))*Rp(i)*Rz(i)*besseljc(q.*Re);
            V = V + pi*(rho(i) - rho(i+1))*Rp(i)*Rz(i);
        %    fq = fq + besseljc(q.*Re);
        end
        
        f = f + gDistR(pd)*abs(fq).^2*I0;
        Vall = Vall + V^2*gDistR(pd);
    end
    f = f/Vall;
else
    
   f=[];
   name='Hollow Ellipsoid Cylinder';
   
    pnames=str2mat('qpORqz', 'Rp1', 'Rp2', 'Rp3', 'Rz1', 'Rz2', 'Rz3', 'rho1', 'rho2', 'rho3', 'I0', 'dR');
	if flag==1, pin=[1, 15, 68, 85, 22, 56, 73.5, 0.1, 0.5, 0.7, 1, 0]; else pin = p; end
end           