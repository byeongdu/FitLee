function [Pq, contrastV] = multilayersphere2(q, radius, density, distribution)
% [y, y1] = multilayersphere(q, radius, density, sig)
% radius = [r11, r12, r13, r14, ..., r1n,...;
%    r11, r12, r13, r14, ..., r1n,...;
%    ...;
%    rm1, rm2, rm3, rm4, ..., rmn];
% density = [d1, d2, ..., d(n+1)];
% d(n+1) : solvent density
% n : number of interfaces.
%
%y = zeros(size(q));
%V = 0;
if (nargin == 3)
	distribution = 1;
end

%Aq = y;
%Iq = y;
m = size(radius, 1);
n = size(radius, 2);
if m ~= numel(distribution)
    error('number of distribution should be the same with row of radius')
end
if n ~= (numel(density)-1)
    error('number of density should be the same with col of radius')
end
numq = numel(q);
density = density(:)';
eden = density(1:end-1) - density(2:end);
Aq = sphereamp(q, radius);
Aq = reshape(Aq, [numq*m, n]);
Aq = Aq*eden(:);
Pqdistr = abs(reshape(Aq, [numq, m])).^2;
Pq = Pqdistr*distribution(:)/sum(distribution);

V = 4*pi/3*eden*transpose(radius).^3;
contrastV = sqrt(V.^2*distribution(:)/sum(distribution));

%for i=1:length(radius)
%    Aq = Aq + (density(i)-density(i+1))*sphereAmp(q, radius(i));
%    V = V + (density(i)-density(i+1))*radius(i).^3/3*4*pi;
%end
