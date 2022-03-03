function [Formfactor, name, pnames, pin] = cylindertype(q, p, flag)
% this is for cylinder type particle : cylindertype(q, p)
% p = [R, dR, H/R, orientationtype, I0) 
% 0: random, 1: along radius: 2: along Height

[xq, yq] = size(q);

if xq < yq 
    q = q';     % q should be column vector
end

if nargin == 2
    
    R = p(1);dR = p(2);
    HoverR = p(3);
    type = p(4); % 0: random, 1: along radius: 2: along Height
    I0 = p(5);
    if (dR ~= 0)
        [gDistR, R] = schultzdist99(R, dR, 10);gDistR=gDistR*(R(2)-R(1));
    else
        gDistR = 1;
    end
    L = HoverR*R;
    Rall = R;Lall = L;
    Y = zeros(size(q));

for j = 1:length(Rall)
    R = Rall(j);
    L = Lall(j);
    switch type
    case 0
    
%        alpha = 0:pi/50:pi/2; % alpha is low vector 
            % Hashimoto, (0 to pi) integral range
%        xx = 1:length(q);
%       [alpha, A] = meshgrid(0:pi/1000:pi, q');
         Formfactor = zeros(size(q));
        dalpha = pi/1000;
        for alpha = 0:dalpha:pi
%            P = (2*pi*R.^2*L*besseljc(R.*A.*sin(alpha)).*sinc(L*A.*cos(alpha)/2)).^2.*sin(alpha);
             Formfactor = Formfactor+(2*pi*R.^2*L*besseljc(R*q*sin(alpha)).*sinc(L*q*cos(alpha)/2)).^2*sin(alpha)*dalpha;
%            P = P+(2*pi*R.^2*L*sinc(L*q*cos(alpha)/2)).^2*sin(alpha);
%            P = (2*besseljc(R*q*sin(alpha))).^2.*sin(alpha);
        end
%            Formfactor = sum(P, 2);
            Formfactor = Formfactor;
            %Formfactor = Formfactor*I0;
    case 1
        qp = q;qz=0.00001;
        Formfactor = abs(corn_vert_type(qp, qz, [R, L, 90])).^2';
    case 2
        qz = q;qp=0.00001;
        Formfactor = abs(corn_vert_type(qp, qz, [R, L, 90])).^2;
    end
    if dR ~= 0
        Y = Y + reshape(Formfactor, size(Y))*gDistR(j);
    else
        Y = Y + Formfactor(:);
    end
end
    Formfactor = Y*I0;
else
    
   Formfactor=[];
   name='Cylinder Model';
   
    pnames=str2mat('R','dR', 'HoverR','type','I0');
	if flag==1, pin=[10 0.1 2 0 1]; else pin = p; end
end       