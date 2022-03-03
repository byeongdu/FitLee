function [qx, qy, qz] = spiral_on_sphere(q, N)
if nargin < 2
    N=21;
end
N = fix(N);
k=1:N;
k=fliplr(k);
hk=-1+(2.*(k-1))./(N-1);
Thetak=acos(hk);
Thetak=fliplr(Thetak);
Phik(1)=0;
for n=2:N-1;
   Phik(n)=mod((Phik(n-1)+(3.6./sqrt(N)).*1/sqrt(1-hk(n).^2)),2*pi);
end
Phik(1)=0;
Phik(N)=2*pi;
%R=ones(size(Thetak));
%[x1 y1 z1]=sph2cart(Thetak,Phik,R);
x = cos(Phik).*sin(Thetak);
y = sin(Phik).*sin(Thetak);
z = cos(Thetak);
x = x(:);y=y(:);z=z(:);
qx = x*q;
qy = y*q;
qz = z*q;