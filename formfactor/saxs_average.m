function [yall, avgFo, avgPqo, avgFq_os] = saxs_average(q, funcname, parameters, nr, method, Rotmat, edgelength, varargin)
% [y, y2, y3, y4] = saxs_average(q, funcname, parameters, nr, method, Rotmat)
% Orientational average for particle scattering
% Input
% q : the scattering vector, a column vector or (qxy, qz) or (qx, qy, qz)
% funcname: the name of a function that will be used. for instance,
% saxsoctahedron
% parameters: parameters for the form factor such as [R, H...] or a cell
% nr : number distribution function, when parameters contains multiple
% rows, numel(nr) == size(parameters, 1), or a cell.
% When you calculate polydisperse particles, put one paritcles at a time.
%
% method: integration method
%   1: spiral average
%   2: 3D average; q should be |q|
%   3: 2D average (xyplane); q should be qxy. or (qxy, qz)
%   4: 1D average (x direction); q should be qx. or (qx, qy, qz)
% Rotmat: cell array of rotation matrices, which is required for 2D and 1D.
%
% Output
% y : <|F|^2>_orient&size
% y2 : <F>_orient : If nr is empty or 1 (monodisperse), the orientational average
% will be done. 
% y3 : <|F|^2>_size: Size average will be done.
% y4 : <F>_orient&size
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
yall = [];
avgFo = [];
avgPqo = [];
avgFq_os = [];

if nargin < 4
    sizeaverage = 0;
    nr = ones(nrp, 1);
else
    sizeaverage = 1;
end
if nargin >= 4
    sizeaverage = numel(nr)-1;
end
% if iscell(nr)
%     sizeaverage = cellfun(@(x) numel(x)-1, nr);
% else
%     sizeaverage = numel(nr)-1;
% end
ispolydisperse = 0;

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
%    n = nr/sum(nr);
    n = nr;
    np = numel(nr);
end
if np > 1
    ispolydisperse = 1;
end

if nargout < 2
    averageoverF = 0;
elseif nargout == 2
    averageoverF = 1;
elseif nargout >= 3
    averageoverF = 2;
end

%strcmd = ['F = @', funcname, ';'];
%eval(strcmd);
N = 2^6;
%N = 2^7;
% if ~contains(funcname, 'diazc')
%     MAXpntat1compute = 25600*2; % This number provides the fast calculation.
% else
    MAXpntat1compute = 25600*2; % This number provides the fast calculation.
% end
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

for seg=1:qseg
    if seg < qseg
        q = QQ((seg-1)*NumQseg+1:NumQseg*seg, :);
    else
        q = QQ((seg-1)*NumQseg+1:end, :);
    end
    [y, avgF1, avgPq, avgF_so] = runcal;
    yall = [yall;y];
    avgFo = [avgFo; avgF1];
    avgPqo = [avgPqo; avgPq];
    avgFq_os = [avgFq_os; avgF_so];

    % Display progress every 10 computation.
    disp_progress = 1;
    if qseg>10
        disp_progress = 0;
        if mod(seg, 10) == 0
            disp_progress = 1;
        end
    end
    if disp_progress == 1
        fprintf('Running seg %i/%i, numel of q = %i\n', seg, qseg, numq)
    end
end

%end

function [Iqtemp, avgFq_o, avgPqo, avgF_so] = runcal
    avgFq_o = [];
    avgF_so = [];
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
                qz = repmat(q(:,2), 1, size(qx, 2));
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

    %if numel(Rm > 1)
    %    [R, qx] = meshgrid(Rm, qx);
    %    [H, qy] = meshgrid(Hm, qy);
    %    [~, qz] = meshgrid(Rm, qz);
    %else
    %    R = Rm;
    %    H = Hm;
    %end

    %F1 = feval(F, qx, qy, qz, parameters);
    if iscell(nr)
        n = [];
        F1 = zeros(numel(qx), numel(nr));nn=0;nnold = 1;
        %Iqtemp = zeros(sizeqx(1), np);
        %avgF2q = zeros(sizeqx(1), np);
        for inr=1:numel(nr)
            n = nr{inr};
            if strcmp(funcname{inr}, 'sphere')
                Q = sqrt(qx.^2+qy.^2+qz.^2);
                [~, ~,~,~, ~, mF] = SchultzsphereFun(Q, parameters{inr}(1), parameters{inr}(2));
                Ftemp = mF;
                %Ftemp = feval(funcname{inr}, qx, qy, qz, parameters{inr});
            elseif strcmp(funcname{inr}, 'diaz') | strcmp(funcname{inr}, 'diazc')
                Q = [qx, qy, qz]*Rotmat{inr};
                qx=Q(:,1);
                qy=Q(:,2);
                qz=Q(:,3);
                if ispolydisperse==0
                    Ftemp = feval(qx, qy, qz, funcname{inr}, parameters{inr});
                else
                    Ftemp = feval(qx, qy, qz, funcname{inr}, [parameters{inr}; edgelength]);
                end                
            else
                Q = [qx, qy, qz]*Rotmat{inr};
                qx=Q(:,1);
                qy=Q(:,2);
                qz=Q(:,3);
                if ispolydisperse==0
                    if numel(varargin)==0
                        Ftemp = feval(funcname{inr}, qx, qy, qz, parameters{inr});
                    else
                        if numel(varargin)==2
                            Ftemp = feval(funcname{inr}, qx, qy, qz, parameters{inr}, varargin{1}, varargin{2});
                        end
                    end
                else
                    if numel(varargin)==0
                        Ftemp = feval(funcname{inr}, qx, qy, qz, [parameters{inr}; edgelength]);
                    else
                        if numel(varargin)==2
                            Ftemp = feval(funcname{inr}, qx, qy, qz, [parameters{inr}; edgelength], varargin{1}, varargin{2});
                        end
                    end
                end
            end
            %nn = nn+numel(nr{inr});
            %F1(:,nnold:nn) = Ftemp;
            F1 = Ftemp;
            [Iqtemp(:,inr), avgFq_o(:,inr), avgPqo(:,inr), avgF_so(:,inr)] = calaverage;
            %F1(:,inr) = Ftemp;
            n(inr) = sum(nr{inr});
            nnold = nn+1;
            %n = [n, nr{i}];
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
        if strcmp(funcname, 'diaz') | strcmp(funcname, 'diazc')
            if numel(parameters) == 1
                F1 = feval(funcname, qx, qy, qz, parameters);
            else
                if ispolydisperse==0
%                    switch numel(parameters)
%                        case 2
                            F1 = feval(funcname, qx, qy, qz, parameters{1}, parameters{2});
%                        case 3
%                            F1 = feval(funcname, qx, qy, qz, parameters{1}, parameters{2}, parameters{3});
%                    end
                else
                    fprintf('No polydispersity for this function yet.\n');
                    error;
                end 
            end
        else
            if ispolydisperse == 0
                if numel(varargin) ==0
                    F1 = feval(funcname, qx, qy, qz, parameters);
                else
                    if numel(varargin) == 2
                        F1 = feval(funcname, qx, qy, qz, parameters, varargin{1}, varargin{2});
                    end
                end
            else
%                F1 = feval(funcname, qx, qy, qz, [parameters; edgelength]);
                for indpoly=1:size(parameters, 1)
                    if numel(varargin) ==0 
                        F1(:,indpoly) = feval(funcname, qx, qy, qz, parameters(indpoly, :));
                    else
                        if numel(varargin)==2
                            F1(:,indpoly) = feval(funcname, qx, qy, qz, parameters(indpoly, :), varargin{1}, varargin{2});
                        end
                    end
                end
                %F1(:,indpoly+1) = feval(funcname, qx, qy, qz, edgelength);
            end
        end
        [Iqtemp, avgFq_o, avgPqo, avgF_so] = calaverage;
    end
%[F1, F2] = feval(F, qx, qy, qz, parameters);

function [Iq, avgF, avgPq, avgF_so] = calaverage
    % Iq : <|F|^2_bar> - if n(r) is provided, both size and orientation
    % averaged.
    % avgF : only orientation averaged F(q) no matter n(r) is provided or
    % not.
    % avgPq: only orientation averaged P(q)- monodisperse
    % avgF_so : <F_bar> - if n(r) is provided, both size and orientation
    % averaged.
    Fq_mono = [];
    if ispolydisperse 
        if size(F1, 2) > numel(n)
            Fq_mono = F1(:, end);
            F1 = F1(:, 1:end-1);
        else
            Fq_mono = F1(:, end);
        end
    end
    avgF = [];
    avgF2 = [];
    avgF_so = [];
    Iq = [];
    
     if sizeaverage == 0 % only orientational average
         Iq = zeros(sizeqx(1), np);
         avgF = zeros(sizeqx(1), np);
         for k = 1:np
             FF = F1(:,k);
             Fq2 = abs(FF).^2;
             Iqt = reshape(Fq2, sizeqx).*sinth;
             Fqt = reshape(FF, sizeqx).*sinth;
             switch method
                 case 1
                     Iq(:,k) = mean(Iqt, 2);
                     if averageoverF >= 1
                        avgF(:,k) = mean(Fqt, 2);
                     end
                 case 2
                     Iq(:,k) = trapz(trapz(Iqt, 3), 2)*dtheta*dphi/IA;
                     if averageoverF >= 1
                        avgF(:,k) = trapz(trapz(Fqt, 3), 2)*dtheta*dphi/IA;
                     end
                 case 3
                     Iq(:,k) = trapz(Iqt, 2)*dtheta*dphi/IA;
                     if averageoverF >= 1
                        avgF(:,k) = trapz(Fqt, 2)*dtheta*dphi/IA;
                     end
                 case 4
                     Iq = Iqt;
                     if averageoverF >= 1
                        avgF = Fqt;
                     end
                     %y = Iq;
                     %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
             end
         end
         %
         avgF_so = avgF;
         avgPq = Iq;
         return;
     end

     if sizeaverage > 0 % both size and orientational average
         % the number of np should be 1......
         [F1x, F1y] = size(F1);
         [nx, ny] = size(n(:));
         if F1y ~= nx
            Iq = sum(abs(F1).^2*n(:)',2); % size average for I(q)
            Fq = sum(F1*n(:)',2); % size average for F(q)
         else
            Iq = sum(abs(F1).^2*n(:),2); % size average for I(q)
            Fq = sum(F1*n(:),2); % size average for F(q)
         end
        Iq = reshape(Iq, sizeqx);
        Fq = reshape(Fq, sizeqx); 
        Fqt = reshape(Fq_mono, sizeqx).*sinth;
        Iqt = abs(Fqt).^2;
    %    Iq = abs(F1).^2.*sinth;
        Iq = Iq.*sinth;
        Fq = Fq.*sinth;
        if method == 1
            Iq = mean(Iq,2);
            avgF = mean(Fqt, 2);
            avgF_so = mean(Fq, 2);
            avgPq = mean(Iqt, 2);
        else
            if (method==2)
                Iq = trapz(trapz(Iq, 3), 2)*dtheta*dphi/IA;
                avgF = trapz(trapz(Fqt, 3), 2)*dtheta*dphi/IA;
                avgPq = trapz(trapz(Iqt, 3), 2)*dtheta*dphi/IA;
                avgF_so = trapz(trapz(Fq, 3), 2)*dtheta*dphi/IA;
                %Iq = sum(sum(Iq, 3), 2)*dtheta*dphi*2/pi;
            elseif (method==3)
                Iq = trapz(Iq, 2)*dtheta*dphi/IA;
                avgF = trapz(Fqt, 2)*dtheta*dphi/IA;
                avgPq = trapz(Iqt, 2)*dtheta*dphi/IA;
                avgF_so = trapz(Fq, 2)*dtheta*dphi/IA;
            elseif (method==4)
                avgF = Fqt;
                avgPq = Iqt;
                avgF_so = Fq;
                %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
            end
        end
        %avgF = avgF2;

     end
%     if averageoverF == 2 % average over F(qx,qy,qz) for size first, then perform orientational average on its absolute square. 
%         if size(F1, 2)==numel(n)
%             Fq = sum(F1*n(:),2);
%         else
%             Fq = sum(F1*n(:)',2);
%         end
%         Fq = reshape(Fq, sizeqx);
%         Iqt = Fq.*sinth;
%         switch method
%             case 1
%                 avgF2 = mean(Iqt, 2);
%             case 2
%                 avgF2 = trapz(trapz(Iqt, 3), 2)*dtheta*dphi/IA;
%             case 3
%                 avgF2 = trapz(Iqt, 2)*dtheta*dphi/IA;
%             case 4
%                 avgF2 = Iqt;
%                 %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
%         end
%     end

% The following was original code....    
%     if averageoverF == 1 % average over F(qx,qy,qz) for size first, then perform orientational average on its absolute square. 
%         if size(F1, 2)==numel(n)
%             Fq = sum(F1*n(:),2);
%         else
%             Fq = sum(F1*n(:)',2);
%         end
%         Fq = reshape(Fq, sizeqx);
%         %Iqt = abs(Fq).^2.*sinth;
%         Iqt = abs(Fq).^2.*sinth;
%         switch method
%             case 1
%                 IqavgF2q = mean(Iqt, 2);
%             case 2
%                 IqavgF2q = trapz(trapz(Iqt, 3), 2)*dtheta*dphi/IA;
%             case 3
%                 IqavgF2q = trapz(Iqt, 2)*dtheta*dphi/IA;
%             case 4
%                 IqavgF2q = Iqt;
%                 %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
%         end
%     end
    
end
end
end
    


% if np == 1
%     F1 = reshape(F1, sizeqx);
%     y = abs(F1).^2.*sinth;
%     if method == 1
%         y = mean(y,2);
%     else
%         if (method==2)
%             y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
%         elseif (method==3)
%             y = trapz(y, 2)*dtheta*dphi/pi;
%         elseif (method==4)
%             %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
%         end
%     end
% else
%     if sizeaverage == 0 % only orientational average
%         y = zeros(numq, np);
%         for k = 1:np
%             Fq2 = abs(F1(:,k)).^2;
%             Iq = reshape(Fq2, sizeqx).*sinth;
%             switch method
%                 case 1
%                     y(:,k) = mean(Iq, 2);
%                 case 2
%                     y(:,k) = trapz(trapz(Iq, 3), 2)*dtheta*dphi*2/pi;
%                 case 3
%                     y(:,k) = trapz(Iq, 2)*dtheta*dphi/pi;
%                 case 4
%                     y = Iq;
%                     %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
%             end
%         end
%     end
%     if sizeaverage == 1 % average over |F(qx,qy,qz)|^2 for size and orientation.
%         if size(F1, 2)==numel(n)
%             Fq2 = sum(abs(F1).^2*n(:), 2);
%         else
%             Fq2 = sum(abs(F1).^2*n(:)', 2);
%         end
%         Iq = reshape(Fq2, sizeqx).*sinth;
%             switch method
%                 case 1
%                     y = mean(Iq, 2);
%                 case 2
%                     y = trapz(trapz(Iq, 3), 2)*dtheta*dphi*2/pi;
%                 case 3
%                     y = trapz(Iq, 2)*dtheta*dphi/pi;
%                 case 4
%                     y = Iq;
%                     %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
%             end
%     end
%     if averageoverF == 1 % average over F(qx,qy,qz) for size first, then perform orientational average on its absolute square. 
%         if size(F1, 2)==numel(n)
%             Fq = sum(F1*n(:),2);
%         else
%             Fq = sum(F1*n(:)',2);
%         end
%         Fq = reshape(Fq, sizeqx);
%         Iq = abs(Fq).^2.*sinth;
%         switch method
%             case 1
%                 avgF2q = mean(Iq, 2);
%             case 2
%                 avgF2q = trapz(trapz(Iq, 3), 2)*dtheta*dphi*2/pi;
%             case 3
%                 avgF2q = trapz(Iq, 2)*dtheta*dphi/pi;
%             case 4
%                 avgF2q = Iq;
%                 %y = trapz(trapz(y, 3), 2)*dtheta*dphi*2/pi;
%         end
%     end
% end

