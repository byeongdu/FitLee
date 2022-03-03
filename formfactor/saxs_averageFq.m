function Fqmean = saxs_averageFq(q, funcname, parameters, nr, method, Rotmat)
% Fqmean = saxs_averageFq(q, funcname, parameters, nr, method, Rotmat)
% Orientational or distribution average for particle scattering F(q).
% It will return <F(q)>. See also saxs_average.m, which returns <|F(q)|^2>
% and |<F(q)>|^2.
% 
% Input
% q : the scattering vector, a column vector
% funcname: the name of a function that will be used. for instance,
% saxsoctahedron
% parameters: parameters for the form factor such as [R, H...] or a cell
% nr : number distribution function, when parameters contains multiple
% rows, numel(nr) == size(parameters, 1), or a cell.
% method: integration method
%   1: spiral average
%   2: 3D average; q should be |q|
%   3: 2D average (xyplane); q should be qxy. or (qxy, qz)
%   4: 1D average (x direction); q should be qx. or (qx, qy, qz)
% Rotmat: cell array of rotation matrices, which is required for 2D and 1D.
%
% Output
% <F(q)>
%
% The argument, parameters, should be a row vector for a particle
% If there are multiple particles are concerned,
% It would be stacked into multiple rows for instance,
% [R1, H1; R2, H2; R3, H3;...]
% Some form factor functions may allow an input of multiple particles, 
% some may not.

if iscell(parameters)
    nrp = size(parameters{1}, 1);
else
    nrp = size(parameters, 1);
end
Fqmean = [];

if nargin < 4
    sizeaverage = 0;
    nr = ones(nrp, 1);
else
    sizeaverage = 1;
end

if nargin<5
    method = 1;
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
N = 2^7;
MAXpntat1compute = 25600*2; % This number provides the fast calculation.

%if N^2*nrp*numq > MAXpntat1compute
if size(q,1)==1 & size(q,2)>1
    q = q(:);
end
QQ = q;
numq = numel(q(:,1));

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

if strfind(funcname, 'pdbFq');
    NumQseg = 1;
    qseg = numq;
end

tic
for seg=1:qseg
    if seg < qseg
        q = QQ((seg-1)*NumQseg+1:NumQseg*seg);
    else
        q = QQ((seg-1)*NumQseg+1:end);
    end
    Fq = runcal;
    Fqmean = [Fqmean;Fq];
    fprintf('Running seg %i/%i, numel of q = %i\n', seg, qseg, numq)
end
toc    
%end

function Fqmean = runcal
    avgF2q = [];
    %strcmd = ['F = @', funcname, ';'];
    %eval(strcmd);
switch method
    case 1
        N2 = N^2;
        [Qx, Qy, Qz] = spiral_on_sphere(1, N2);
        q = q(:);
        qx = q*Qx';
        qy = q*Qy';
        qz = q*Qz';
        dphi = 1;
        dtheta = pi/2;
        sizeqx = size(qx);
        sinth = 1;
    case 2
        %N = 20;
        dtheta = pi/(2*N);
        dphi = pi/(2*N);
        max_theta = pi/2;
        max_phi = pi/2;
        theta = 0:dtheta:max_theta;
        phi = 0:dphi:max_phi;
        % Integration area
        % When max_theta = pi and max_phi=2*pi, integration can be done for
        % a full circle, then the integration area becomes 4*pi.
        % When max_theta=pi/2 and max_phi=pi/2, integration can be done for
        % a 1/8 circle, then the integration area becomes pi/2.
        % Therefore, 
        % the integration area fraction(IAF) = (max_theta/pi * max_phi/pi / 2)
        % and the integration area = IAF*4*pi;
        IA = (max_theta/pi * max_phi/pi / 2) * 4*pi;
        %theta = 0:dtheta:pi;
        %phi = 0:dphi:2*pi;
        [q, theta, phi] = ndgrid(q, theta, phi);
        qx = -q.*sin(theta).*cos(phi);
        qy = q.*sin(theta).*sin(phi);
        qz = q.*cos(theta);
        sinth = sin(theta);
        sizeqx = size(qx);
    case 3 %2D average
        dtheta = pi/(2*N);
        max_theta = pi; % When max_theta = pi, it only integrate 1/2.
        theta = 0:dtheta:max_theta;
        IA = max_theta;
        %theta = 0:dtheta:pi;
        %phi = 0:dphi:2*pi;
        if size(q,2) == 1
            [r, theta] = ndgrid(q, theta);
            qx = r.*cos(theta);
            qy = r.*sin(theta);
            qz = zeros(size(qx));
        else
            [r, theta] = ndgrid(q(:,1), theta);
            qx = r.*cos(theta);
            qy = r.*sin(theta);
            qz = q(:,2);
        end    
        sizeqx = size(qx);
        sinth = 1;
        dphi = 1;
    case 4 %1D average
        %theta = 0:dtheta:pi;
        %phi = 0:dphi:2*pi;
        IA = 1;
        if size(q,2) == 1
            qx = q;
            qy = zeros(size(qx));
            qz = qy;
        else
            qx = q(:,1);
            qy = q(:,2);
            qz = q(:,3);
        end
        sizeqx = size(qx);
        sinth = 1;
        dphi = 1;
        dtheta = 1;
end

qx = qx(:);
qy = qy(:);
qz = qz(:);

if iscell(nr)
    n = [];
    F1 = zeros(numel(qx), numel(nr));
    for inr=1:numel(nr)
        n = nr{inr};
        if strcmp(funcname{inr}, 'sphere')
            Q = sqrt(qx.^2+qy.^2+qz.^2);
            [~, ~,~,~, ~, mF] = SchultzsphereFun(Q, parameters{inr}(1), parameters{inr}(2));
            Ftemp = mF;
        else
            Q = [qx, qy, qz]*Rotmat{inr};
            qx=Q(:,1);
            qy=Q(:,2);
            qz=Q(:,3);
            Ftemp = feval(funcname{inr}, qx, qy, qz, parameters{inr});
        end
        F1 = Ftemp;
        Fqmean(:,inr) = calaverage;
        n(inr) = sum(nr{inr});
    end
else
    if iscell(Rotmat)
        R = Rotmat{1};
    else
        R = Rotmat;
    end
    Q = [qx, qy, qz]*R;
    qx=Q(:,1);
    qy=Q(:,2);
    qz=Q(:,3);
    F1 = feval(funcname, qx, qy, qz, parameters);
    Fqmean = calaverage;
end
%[F1, F2] = feval(F, qx, qy, qz, parameters);

function Fq = calaverage
    % size and orientation average for |F|^2
    Fq = [];
    
     if sizeaverage == 0 % only orientational average
         Fq = zeros(numq, np);
         for k = 1:np
             Fqt = F1(:,k);
             Fqt = reshape(Fqt, sizeqx).*sinth;
             switch method
                 case 1
                     Fq(:,k) = mean(Fqt, 2);
                 case 2
                     Fq(:,k) = trapz(trapz(Fqt, 3), 2)*dtheta*dphi/IA;
                 case 3
                     Fq(:,k) = trapz(Fqt, 2)*dtheta*dphi/IA;
                 case 4
                     Fq = Fqt;
                     %y = Iq;
                     %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
             end
         end
     end

     if sizeaverage == 1 % both size and orientational average
         try
            Fqt = sum(F1*n(:),2); % size average
         catch
             Fqt = sum(F1*n,2);
         end
        Fqt = reshape(Fqt, sizeqx);
    %    Iq = abs(F1).^2.*sinth;
        Fqt = Fqt.*sinth;
        if method == 1
            Fq = mean(Fqt,2);
        else
            if (method==2)
                Fq = trapz(trapz(Fqt, 3), 2)*dtheta*dphi/IA;
                %Iq = sum(sum(Iq, 3), 2)*dtheta*dphi*2/pi;
            elseif (method==3)
                Fq = trapz(Fqt, 2)*dtheta*dphi/IA;
            elseif (method==4)
                %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
            end
        end
     end
end
end
end
    
