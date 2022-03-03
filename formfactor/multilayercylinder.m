function [Pq, Aq] = multilayercylinder(q, radius, density, distribution, L)
% [y, y1] = multilayersphere(q, radius, density)
% radius = [r11, r12, r13, r14, ..., r1n,...;
%    r11, r12, r13, r14, ..., r1n,...;
%    ...;
%    rm1, rm2, rm3, rm4, ..., rmn];
% density = [d1, d2, ..., d(n+1)];
% density(n+1) : solvent density
% n : number of interfaces.
% distribution: its length should be m.

if (nargin == 3)
	distribution = 1;
    L = 1;
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
density = density(:)';
eden = density(1:end-1) - density(2:end);
q = q(:);
q_ones = ones(size(q));
Pq = zeros(size(q));
for i=1:m
    R = radius(i, :); % n shell
    Rmat = q_ones*R;
    Eden_ones = q_ones*eden;
    Aq = sum(2*pi*L*Eden_ones.*Rmat.^2.*besseljc(q*R), 2);
    Pq = Pq + distribution(i)*abs(Aq).^2;
end
Pq = Pq./q.^2;
