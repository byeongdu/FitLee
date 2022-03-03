function y = sq2Gr(q, sq, r, np)
% this pair distribution calculation is from the atomic PDF..
% g(r) = 1 + 1/(2*pi^2)/n_p * int q^2[S(q)-1]sinc(qr)dq;
% G(r) = 4*pi*n_p*(g(r)-1)*r;
% Gr = sq2Gr(q, sq, r)
% return g(r) when np is provided
% return G(r) when np is not provided.

% Byeongdu Lee

% if q min is not 0, extrapolate Sq to q = 0. This is an important step to
% reduce Fourier ripple.
if q(1) > 0
    dq = q(2) - q(1);
    qext = 0:dq:(q(1)-dq);
    %Sqext = qext*sq(1)/q(1);
    Sqext = sq(1)*ones(size(qext));
    q = [qext(:); q(:)];
    sq = [Sqext(:); sq(:)];
end

t = isinf(sq);
sq(t) = [];q(t)=[];
dq = min(diff(q));
q2 = min(q):dq:max(q);
sq = interp1(q, sq, q2); q = q2;
%dq = q(2)-q(1);
%dsp = 0.2732;
%R = dsp/2;

%volumefractionofparticle = 0.68;
%volumeofaparticle = (4/3*pi*R).^3;
%particlenumberdensity = volumefractionofparticle/volumeofaparticle; 
    % closed packed disordered fluid
    % PRB 59, 22, 14191, ref 35 threrein
%y = sum(r2.*q2.*pi.*(sq2-ones(size(sq2))).*sin(q2.*r2)*dq, 2);
y = zeros(size(r));
for i=1:length(r)
    %1/2/pi^2*q.*(sq-1).*sin(q*r(i))*dq;
    %y(i) = sum(1./(r(i)+eps)*1/2/pi^2*q.*(sq-1).*sin(q*r(i))*dq);
    y(i) = 1/(2*pi^2)*sum(q.^2.*(sq-1).*sinc(q*r(i)))*dq;
    if nargin==3
        y(i) = y(i)*4*pi*r(i);
    elseif nargin == 4
        y(i) = y(i)/np+1;
    end
        
end
%y = 1 + y/particlenumberdensity;