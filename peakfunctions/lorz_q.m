function I=lorz_q(varargin)
% I=lorz_q(x,y,p)
% p = [1/Q.size(1), qpz(i,1), qpz(i,2), orientw, w]
% I=lorz_q(x,y,z,p)
% p = [1/Q.size(1), qvalue(i,1), qvalue(i,2), qvalue(i,3), orientw, w]
% lorentz function along the q direction.
% gaussian function for the orientational distribution.
% y=lorz_q(x,y,p)
% http://en.wikipedia.org/wiki/Lorentzian_function
% p = [ Area Centrex centery Width width_q ]

% Author: B. Lee
% when p(1) =1, the area should be 1.
cs = 0;
if numel(varargin) == 3
    x = varargin{1};
    y = varargin{2};
    p = varargin{3};
    z = zeros(size(x));
elseif numel(varargin) == 4
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    p = varargin{4};
    cs = 1;
end 
switch cs
    case 0
        p0 = sqrt(sum(p(2:3).^2));
        xt = sqrt(x.^2+y.^2);               % qz direction = 1
        I1 = lorza(xt, [p(1), p0, p(5), 0]); % qxy direction
        Ig = lorza2D(x, y, [1, p(2), p(3), p(4)]); % azimuthal direction.
    case 1
        p0 = sqrt(sum(p(2:4).^2));
        xt = sqrt(x.^2+y.^2+z.^2);
        I1 = lorza(xt, [p(1), p0, p(6), 0]); % q direction
        Ig = lorza3D(x, y, z, [1, p(2:5)]);
end
%Ig=gauss2d(x, y, [p(1), p(2), p(3), p(4), p(4), 0]);
I = I1.*Ig;