function [P, fq] = Cylinderlayered(q, p, flag)
% Formfactor of multiple layered cylinder
% [P, name, pnames, pin] = Cylinderlayered(q, p, flag)
% currently, only fq are calculated, avFq2 = abs(fq).^2;
%   p = [R1, ... R5], dR, [rho1, ... rho5], L, I0, bkg
%
% L. Ramos et al. Eur. Phys. J. B 1, 319(1998)
%n = length(rho);
%sumBottom = 0;
%for i=1:n-1
%    sumBottom = sumBottom + (rho(i) - rho(i+1))*R(i)^2;
%end

[xq, yq] = size(q);

if xq < yq 
    q = q';     % q should be column vector
end

if nargin == 2
    
    Rinner = abs(p(1));
    R2 = abs(p(2));
    R3 = abs(p(3));
    R4 = abs(p(4));
    Rout = abs(p(5));
    dR = abs(p(6));
    rho1 = p(7);
    rho2 = p(8);
    rho3 = p(9);
    rho4 = p(10);
    rho5 = p(11);
    L = abs(p(12));
    I0 = p(13);
%    bkg = p(14);
bkg = 0;
    rhoout = 0;
    rho6 = 0;
    RR = Rout;
%    fl = cylindertype3(q, L);   
    if (dR ~= 0)
        [gDistR, Rout] = schultzdist99(Rout, dR, 10);gDistR=gDistR*(Rout(2)-Rout(1));
    end
    
    P = zeros(size(q));
%    dalpha = pi/1000;
    %rho4 = rho3;
    %rho5 = rho3;
    %rho6 = 0;
    rho = [rho1, rho2, rho3, rho4, rho5, rho6];

    V = 0;
    for j = 1:length(Rout)
        R = [Rinner, R2, R3, R4, RR];
        R = R/RR*Rout(j);
        fq = 0;
        for i=1:length(R)
            fq = fq + 2*pi*L*(rho(i) - rho(i+1))*R(i)^2*besseljc(q*R(i));
        end
        Y = abs(fq).^2;
        Vr = 0;
        for i=1:length(R)
            Vr = Vr + rho(i)*R(i);
        end
        
        if (dR ~= 0)
            P = P + Y*gDistR(j);
        else
            P = P + Y;
        end
        V = V + Vr;
    end
    fl = 1./q;
%    fl = cylindertype3(q, L);
%    fl = ones(size(q));
%    P = P.*fl/(pi*(p(1)^2-p(2)^2)*L)*I0 + bkg;
 %   P = P/(L).*q.^2*I0 + bkg;
        P = P/L*I0.*fl + bkg;
else
    
    P=[];
    name='Multi-Cylinder Model';
    pnames=str2mat('Rinner', 'R2', 'R3', 'R4', 'Rout', 'dR', 'rho1', 'rho2', 'rho3',  'rho4', 'rho5', 'L', 'I0', 'bkg');
    if flag==1, pin=[5, 10, 20, 50, 83, 4, 0, 0, 0, 0, 1, 1, 0]; else pin = p; end
end           