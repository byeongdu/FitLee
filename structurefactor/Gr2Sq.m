function Sq = Gr2Sq(r, gGr, q, np)
% Sq = Gr2Sq(r, gGr, q, np)
% If np is provided,gGr is G(r) that is 4*pi*np*r*(g(r)-1);
% otherwise, gGr is assumed as g(r).
% see pdf2sq.m

dr = min(diff(r));
r2 = 0:dr:r(end);
gGr = interp1(r, gGr, r2);
r = r2(:)';%dr = r(2)-r(1);
gGr = gGr(:)';
q = q(:);
tq0 = q==0;q(tq0) = [];
tq = ones(size(q));
if nargin<4
    gr = gGr/(4*pi*np)./r+1; % if gGr is G(r)
    np = 1;
elseif nargin == 4
    gr = gGr; % if gGr is g(r)
end
gr = (gr-1)*np;
temp = (tq*(gr.*r.^2)).*sinc(q*r)*dr;
temp(isnan(temp))=0;
Sq = 4*pi*sum(temp, 2)+1;
