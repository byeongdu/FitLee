function [yall, avgF2qall] = saxs_average_cylindrical(qx, qz, funcname, parameters, nr, Rotmat)
% [yall, avgF2qall] = saxs_average_cylindrical(qx, qz, funcname, parameters, nr, Rotmat)
% 2D average version of saxs_average
% How to use, take a look at saxs_average.m
% qx and qz should be equal in size.

if iscell(parameters)
    nrp = size(parameters{1}, 1);
else
    nrp = size(parameters, 1);
end
yall = [];
avgF2qall = [];

if nargin < 5
    sizeaverage = 0; % when average over polydisperse particles..
    nr = ones(nrp, 1);
else
    sizeaverage = 1; % when average over monodisperse particles..
end
if nargout < 2
    averageoverF = 0;
else
    averageoverF = 1;
end

if nargin <6
    for i=1:numel(nr)
        Rotmat{i} = eye(3);
    end
end

n = [];
if iscell(nr)
    for i = 1:numel(nr)
        %xxr = [xxr, xr{i}];
        n = [n, nr{i}];
    end
    np = numel(n);
else
    n = nr/sum(nr);
    np = numel(nr);
end

%strcmd = ['F = @', funcname, ';'];
%eval(strcmd);
N = 2^5;
MAXpntat1compute = 25600*2; % This number provides the fast calculation.

%if N^2*nrp*numel(q) > MAXpntat1compute
[inx, iny] = size(qx);
qx = qx(:); qz = qz(:);
QQ1 = qx;
QQ2 = qz;
numq = numel(qx);

qseg = fix(N^2*nrp*numq/MAXpntat1compute);
NumQseg = round(numq/qseg);
if qseg == 0
    NumQseg = numq;
    qseg = 1;
else
    qseg = round(numq/NumQseg);
end

if NumQseg == 0
    NumQseg = 5;
    qseg = round(numq/NumQseg);
end

dtheta = pi/(2*N);
theta = 0:dtheta:pi;

for seg=1:qseg
    if seg < qseg
        q1 = QQ1((seg-1)*NumQseg+1:NumQseg*seg);
        q2 = QQ2((seg-1)*NumQseg+1:NumQseg*seg);
    else
        q1 = QQ1((seg-1)*NumQseg+1:end);
        q2 = QQ2((seg-1)*NumQseg+1:end);
    end
    %runcal;
    [y, avgF2q] = runcal(theta, dtheta, np, q1, q2, nr, funcname, parameters, Rotmat, sizeaverage, averageoverF);
    yall = [yall;y];
    avgF2qall = [avgF2qall;avgF2q];
    %yall = reshape(yall, inx,iny);
    %avgF2qall = reshape(avgF2qall, inx,iny);
    fprintf('Running seg %i/%i, numel of q = %i\n', seg, qseg, numel(q1))
end
    
function [y, avgF2q] = runcal(theta, dtheta, np, q1, q2, nr, funcname, parameters, Rotmat, sizeaverage, averageoverF)
    avgF2q = [];
    sinth = 1;
    dphi = 1;
    
%    case 3 %2D average
        %theta = 0:dtheta:pi;
        %phi = 0:dphi:2*pi;
    [r, theta] = ndgrid(q1, theta);
    qx = r.*cos(theta);
    qy = r.*sin(theta);
    qz = repmat(q2, [1, size(theta, 2)]);
    sizeqx = size(qx);

qx = qx(:);
qy = qy(:);
qz = qz(:);


if iscell(nr)
    n = [];
    F1 = zeros(numel(qx), numel(nr));
    for inr=1:numel(nr)
        if strcmp(funcname{inr}, 'sphere')
            Q = sqrt(qx.^2+qy.^2+qz.^2);
            [~, ~,~,~, ~, mF] = SchultzsphereFun(Q, parameters{inr}(1), parameters{inr}(2));
            Ftemp = mF;
            %Ftemp = feval(funcname{inr}, qx, qy, qz, parameters{inr});
        else
            Q = [qx, qy, qz]*Rotmat{inr};
            qx=Q(:,1);
            qy=Q(:,2);
            qz=Q(:,3);
            Ftemp = feval(funcname{inr}, qx, qy, qz, parameters{inr});
        end
        %nn = nn+numel(nr{inr});
        %F1(:,nnold:nn) = Ftemp;
        F1(:,inr) = Ftemp;
        n(inr) = sum(nr{inr});
        %n = [n, nr{i}];
    end
else
    n = nr;
    Q = [qx, qy, qz]*Rotmat{1};
    qx=Q(:,1);
    qy=Q(:,2);
    qz=Q(:,3);
    F1 = feval(funcname, qx, qy, qz, parameters);
end


if np == 1
    F1 = reshape(F1, sizeqx);
    y = abs(F1).^2.*sinth;
    y = trapz(y, 2)*dtheta*dphi/pi;
else
    if sizeaverage == 0
        y = zeros(numel(F1(:,1)), np);
        for k = 1:np
            Fq2 = abs(F1(:,k)).^2;
            Iq = reshape(Fq2, sizeqx).*sinth;
            y(:,k) = trapz(Iq, 2)*dtheta*dphi/pi;
        end
    end
    if sizeaverage == 1
        if size(F1, 2)==numel(n)
            Fq2 = sum(abs(F1).^2*n(:), 2);
        else
            Fq2 = sum(abs(F1).^2*n(:)', 2);
        end
        Iq = reshape(Fq2, sizeqx).*sinth;
        y = trapz(Iq, 2)*dtheta*dphi/pi;
    end
    if averageoverF == 1
        if size(F1, 2)==numel(n)
            Fq = sum(F1*n(:),2);
        else
            Fq = sum(F1*n(:)',2);
        end
        Fq = reshape(Fq, sizeqx);
        Iq = abs(Fq).^2.*sinth;
        avgF2q = trapz(Iq, 2)*dtheta*dphi/pi;
    end
end

