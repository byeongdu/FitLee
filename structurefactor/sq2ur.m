function [ur, gr, cr] = sq2ur(q, sq, r, numdensity)
% calculate potentional funtion using HNC approximation
% PRB. 59, 14191 (1999)
% smaller q region should be far more important to get correct u(r)
% truncation effect is obviously seen in case not enough small q region.
% Byeongdu Lee
% 2009. 4. 3
%
t = sq == 0;
sq(t) = [];
q(t) = [];

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



pdf = pairdistf(q, sq, r);
r = r(:);
gr = 1 + 1./(numdensity*r).*pdf;
%gr(1:385) = ones(size(gr(1:385)))*gr(386);
k = 1;
while min(gr) < 0
    k = k/2;
    gr = 1 + k./(numdensity*r).*pdf;
end

q = q(:)';   
Sq = sq(:)';
r = r(:);
dr = r(2)-r(1);
temp = (r*((Sq-1)./Sq.*q)).*sin(r*q);
Pdf_r = sum(temp, 2)*dr;
cr = k/(4*pi*numdensity)./r.*Pdf_r;
%cr(1:385) = ones(size(cr(1:385)))*cr(386);

ur = gr - 1 - cr - log(gr);