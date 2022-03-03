function y=debye(q, atm, formf)
% Debye scattering formula for a single type scatter
% y=debye(q, R, formf)
% R is atomic coordinates and kind[x, y, z, kind];
% formfactor
N = numel(atm)/4;
x = atm(:,1);
y = atm(:,2);
z = atm(:,3);
c = atm(:,4);
q = vect2row(q)';
data = zeros(size(q));

for m=1:N
    data=data+formf(:, c(m)).*formf(:, c(m));%*atm(m)*atm(m);
    for k=(m+1):N
        r_jk = sqrt((x(m)-x(k)).^2+(y(m)-y(k)).^2+(z(m)-z(k)).^2);
        data=data + formf(:, c(m)).*formf(:, c(k)).*sinc(q*r_jk)*2;%*atm(m)*atm(k);
    end
end
y = data;

function u = sinc(x)
if x==0
    u = 1;
else
    u = sin(x)./x;
end