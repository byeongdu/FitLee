function [Fqhkls, Pqhkls, Pqall, Fqmall, Pqmall, Fqall] = AnisoFormFactor(particles, hklsorq, cellinfo, dim)
% [Fqhkls, Pqhkls] = AnisoFormFactor(particles, hkls, cellinfo, method)
% Fqhkls = AnisoFormFactor(particles, hkls, cellinfo, method)
% [~, Pq] = AnisoFormFactor(particles, q, cellinfo, method)
% [~, ~, Pqall, Fqmall, Pqmall] = AnisoFormFactor(particles, q, cellinfo, method)
%
% This function is to calculate the form factors of particles in a unit
% cell. 
% Output
%   Fqhkls: F(qx,qy,qz, , N_particles)
%   Pqhkls : sum<|F|^2>
%   Pqall: Pq(q(:), N_particles)
%   Fqmean(q(:), N_particles) : <F(q)> - only when the second input (hkls) is scalar q.
%   Pqmall(q(:), N_particles) : <|F(q)|^2>_orientation : only the
%   orientation is averaged not the size distribution...
%   Fqall(q(:), N_particles) :  <F(q)_bar> (size and orientation averaged)
%           - only when the second input (hkls) is scalar q.
%
% The classical electron radius will be multiplied to the electron density.
%
% Input
%   particles: particles, struct format
%   hklsorq: it coulbe hkls(struct), [qx, qy, qz] (N*3 matrix), or q (a colum vector)
%   cellinfo: cellinfo (struct)
%   method for Pq integration (1~4, 1 and 2 for 3D, 3 for 2D and 4 for 1D)
%
% When inputmode == 3(q for input), then Pqhkls and Pqhkls2 will be the
% integrated intensity for all particles. Otherwise Fqhkls or Pqhkls will
% be that of individual particles. So, SOF is included in the mode 3 only.
%
% particles.shape
%
% for cube
%   edge length of cube = 2*R, where R is particles.radius
%
%
% when anisotropic particles needs to be rotated....
% rotate pyramid 45degree in xy plane.
%x = qx;
%y = qy;
%qx = cos(xyangle)*x + sin(xyangle)*y;
%qy = -sin(xyangle)*x + cos(xyangle)*y;
% rotate pyramid 45degree in yz plane.
%y = qy;
%z = qz;
%qy = cos(yzangle)*y + sin(yzangle)*z;
%qz = -sin(yzangle)*y + cos(yzangle)*z;
%switch cellinfo.dim;
%    case 1
%        dim = 4;
%    case 2
%        dim = 3;
%    case 3
%        dim = 1;
%end
% Unit of Fq: cm.
% Unit of Pq: cm^2.

r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;
qxn = [];qyn = []; qzn = [];

%r_e = r_e*Angstrom2Centimeter;

numparticle = numel(particles);
numP = zeros(numparticle, 1);
Np = 0;
isonlysphere = 1;
if ~iscell(particles)
    p = particles;
    particles = [];
    particles{1} = p;
end

m = 1;
ptinput = [];
for pindx=1:numparticle
    for k=1:numel(particles{pindx}.position)/3
        ptinput = [ptinput; m, pindx, particles{pindx}.position(k, :)];
        m = m + 1;
    end
    numP(pindx) = numel(particles{pindx}.position)/3*particles{pindx}.SOF;
    particles{pindx}.rho = particles{pindx}.rho*r_e;
    Np = Np + numP(pindx);   
    if ~strcmp(particles{pindx}.shape, 'sphere')
        isonlysphere = 0;
    end
end

[Ppos, Ptindex, Pindex] = unique_m(ptinput(:, 3:5)); % Unique by positions.

%numPfrac = numP/Np;
numPfrac = numP;
if nargin < 4
    dim = 1;
end


%typeofinput = 0;

if isstruct(hklsorq)
    hkls = hklsorq;
    typeofinput = 1;
    numq = numel(hkls);
else
    if size(hklsorq, 2) == 3 % (qx, qy, qz)
        hkls = [];
        qx = hklsorq(:,1);
        qy = hklsorq(:,2);
        qz = hklsorq(:,3);
        
        typeofinput = 2;
        numq = numel(qx);
        q = sqrt(qx.^2+qy.^2+qz.^2);
        switch dim
            case 4 % single crystal
                QforPq = [qx, qy, qz];
            case 3 % 2D powder
                QforPq = [sqrt(qx.^2+qy.^2), qz];
            case {1,2}
                QforPq = q;
        end
    elseif size(hklsorq, 2) ==2 % (qxy, qz)
        hkls = [];
        qxy = hklsorq(:,1);
        qz = hklsorq(:,2);
        numq = numel(qxy);
        q = sqrt(qxy.^2+qz.^2);
        switch dim
            case 3 % 2D powder
                QforPq = [qxy, qz];
            case {1,2}
                QforPq = q;
        end
        typeofinput = 4;
    else
        hkls = [];
        q = hklsorq(:);
        numq = numel(q);
        QforPq = q;
        typeofinput = 3;
    end
end

switch typeofinput
    case 1 % when the input is the hkls.
        Fqhkls = zeros(numq, max([hkls.multiplicity]), numel(particles));
        Pqhkls = zeros(numq, numel(particles));
    case 2  % when the input is the q vector or qx,qy,qz
        if isonlysphere
            Fqhkls = zeros(numq, numel(particles{1}.radius), numel(particles));
        else
            Fqhkls = zeros(numq, 1, numel(particles));
        end
    case 3 % when the input is the total q. In this case dim should be 3
        if dim > 2
            cprintf('red', 'Error in AnisoFromFactor.m.\n');
            error('When an array of q is input, orientation type should be 3D random powder.\n');
        end
        Fqhkls = [];
        Pqhkls = zeros(numq, numel(particles));
    case 4 % 2D powder.
        Fqhkls = [];
        Pqhkls = zeros(numq, numel(particles));
end


for pshape = 1:numel(particles)
    particleshape = particles{pshape}.shape;
    isfromparticlemaker = 0;
    if isfield(particles{pshape}, 'edgelength')
        edL = particles{pshape}.edgelength;
        isfromparticlemaker = 1;
    end
    if isfromparticlemaker
        if isfield(particles{pshape}, 'edgelength_sig')
            sig = particles{pshape}.edgelength_sig;
        else
            sig = 0;
        end
    else
        if isfield(particles{pshape}, 'radius_sig')
            sig = particles{pshape}.radius_sig;
        else
            sig = 0;
        end
    end
    
    switch lower(particleshape)
        case 'rho-spherical'
            fffunctionname = 'rho2Fq';
            parameter = particles{pshape}.rhomodel;
            parameter(:,2) = parameter(:,2);
        case 'boxmodel'
            fffunctionname = 'boxFq';
            parameter = particles{pshape}.BOXmodel;
        case 'pdb'
            fffunctionname = 'pdbFq';
            if isfield(particles{pshape}, 'position_atoms')
                parameter = particles{pshape}.position_atoms;
            else
                parameter = particles{pshape}.pdb;
            end
        case 'sphere'
            R = particles{pshape}.radius;
            if isfield(particles{pshape}, 'radius_sig')
                sig = particles{pshape}.radius_sig;
            else
                sig = 0;
            end
            parameter = [R, sig];
            fffunctionname = 'sphere';
            if typeofinput ~= 3
                if ~isempty(hkls)
                    d = [hkls(:).D];
                    qn = 2*pi./d;
                    qn = qn(:);
                else
                    qn = q;
                end
                pt{1} = particles{pshape};
                [Fqtemp, Pqtemp] = SphFormFactorforSq2(qn, pt);
                
%                 Fqtemp = Fqtemp/Angstrom2Centimeter;
%                 Pqtemp = Pqtemp/Angstrom2Centimeter^2;
                
                Fqhkls(1:size(Fqtemp,1), 1, pshape) = Fqtemp;
                Pqhkls(:, pshape) = Pqtemp;
                %continue
            end
        case {'rhombicdodecahedron', 'rhombic dodecahedron', 'rd'}
            if isfromparticlemaker
                R = edL;
                sig = sig;
            else
                R = particles{pshape}.radius;
            end
            parameter = R;
            fffunctionname = 'saxsrhombicdodecahedron';
        case 'octahedron'
            if isfromparticlemaker
                R = edL/2;
                sig = sig/2;
            else
                R = particles{pshape}.radius;
            end
            parameter = R;
            % amplitude
            fffunctionname = 'saxsoctahedron';

        case 'pyramid'
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) == 1
                    H = R*sqrt(2);
                else
                    H = edL(2);
                end
            else
                R = particles{pshape}.radius;
                H = particles{pshape}.height;
            end
            parameter = [R,H];
            % amplitude
            fffunctionname = 'saxspyramid';

        case 'icosahedron'
            if isfromparticlemaker
                R = edL/2;
                sig = sig/2;
            else
                R = particles{pshape}.radius;
            end
            parameter = R;
            % amplitude
            fffunctionname = 'saxsicosahedron';

        case 'truncatedoctahedron'
            tr = [];
            if isfromparticlemaker
                R = edL(1);
                if numel(edL)==2
                    t = edL(2);
                    tr = t*R;
                end
            else
                R = particles{pshape}.radius;
                if isfield(particles{pshape}, 'truncatedfraction')
                    t = particles{pshape}.truncatedfraction;
                    tr = t*R;
                end
            end
            if isempty(tr)
                parameter = R;
            else
                parameter = [R, tr];
            end

            fffunctionname = 'saxstruncatedoctahedron';
            
        case 'facefilledto'
            if isfromparticlemaker
                if numel(edL) ~= 3
                    error("Need three parameter for a facefilledTO")
                else
                    R = edL(1);
                end
            else
                error("Define facefilledTO from particlemaker_PDB.")
            end
            parameter = edL;

            fffunctionname = 'saxsfacefilledto';
            
        case 'pentagonalbipyramid'
            if isfromparticlemaker
                R = edL(1);
                H = edL(2);
            end
            parameter = [R, H];

            fffunctionname = 'saxspentagonalbipyramid';

        case 'longrd'
            H = [];
            if isfromparticlemaker
                R = edL(1);
                if numel(edL) > 1
                    H = edL(2);
                else
                    H = 0;
                end
            else
                R = particles{pshape}.radius(1);
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end
            if isempty(H)
                error('field height is missing')
            end

            parameter = [R, H];
            fffunctionname = 'saxslongRD';

        case 'truncatedcube'
            if isfromparticlemaker
                R = edL;
            else
                R = particles{pshape}.radius;
            end
            parameter = R;
            fffunctionname = 'saxstruncatedcube';
            
        case 'cubooctahedron'
            tr = [];
            if isfromparticlemaker
                R = edL;
            else
                R = particles{pshape}.radius;
                if isfield( particles{pshape}, 'truncatedfraction')
                    f = particles{pshape}.truncatedfraction;
                    tr = f*R;
                end
            end
            parameter = [R, tr];
            fffunctionname = 'saxscubooctahedron';
                        
        case 'cube'
            if isfromparticlemaker
                R = edL(1);
                parameter = [R];
                if numel(edL) == 2
                    H = edL(2);
                    parameter = [R, R, H];
                    sig = [sig(1), sig];
                end
                parameter = parameter/2;
                sig = sig/2; % since radius becomes half.
            else
                R = particles{pshape}.radius(1);
                %H = particles{1}.height;
                parameter = R;
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                    parameter = [R, R, H];
                end
            end
            fffunctionname = 'saxscube';
            
        case 'cylinder'
            H = [];
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) > 1
                    H = edL(2);
                end
                sig = sig/2;
            else
                R = particles{pshape}.radius(1);
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end
            if isempty(H)
                error('field height is missing')
            end

            parameter = [R, H];
            fffunctionname = 'saxscylinder';
        case 'concavecube'
            H = [];
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) > 1
                    H = edL(2);
                end
                sig = sig/2;
            else
                R = particles{pshape}.radius(1);
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end
            if isempty(H)
                error('field height is missing')
            end

            parameter = [edL(1), H];
            fffunctionname = 'saxsconcavecube';
        case 'convexcube'
            H = [];
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) > 1
                    H = edL(2);
                end
                sig = sig/2;
            else
                R = particles{pshape}.radius(1);
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end
            if isempty(H)
                error('field height is missing')
            end

            parameter = [edL(1), H];
            fffunctionname = 'saxsconvexcube';

        case 'thh'
            H = [];
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) > 1
                    H = edL(2)/2;
                end
                sig = sig/2;
            else
                if isfield(particles{pshape}, 'radius')
                    R = particles{pshape}.radius(1);
                end
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end                    
            if isempty(H)
                H = R;
            end

            parameter = [R, H];
            fffunctionname = 'saxsTHH';
    end
% size distribution effect.....

    if ~strcmp(fffunctionname, 'sphere')
        if sum(sig) ~= 0
            if numel(sig) == 1
                [nr, xr] = schultzdist99(parameter(1), sig, 10);nr = nr/sum(nr);
                xr = xr(:);
                nr = nr(:);
            else
                % this will compute all different combination of sig along
                % different dimensions of particles.. takes a long time to
                % compute.
                disp('You choose to compute polydispersity for all dimensions of particles');
                nrt = {};
                xrt = {};
                for sigi = 1:numel(sig)
                    if sig(sigi) == 0
                        nrt{sigi} = 1; xrt{sigi} = parameter(sigi);
                    else
                        [nrt{sigi}, xrt{sigi}] = schultzdist99(parameter(sigi), sig(sigi), 10);
                        nrt{sigi} = nrt{sigi}/sum(nrt{sigi});
                    end
                end
                NonezeroSig = cellfun(@numel, nrt) > 1;
                if sum(NonezeroSig) == 1
                    xr = repmat(parameter, 10, 1);
                    xr(:, NonezeroSig) = xrt{NonezeroSig};
                    %nr = 1/10*ones(10, numel(sig));
                    nr = nrt{NonezeroSig};

                elseif sum(NonezeroSig) > 1

                    [un_par, ~, bin] = unique(parameter, 'stable');
                    xr = repmat(parameter, 100, 1);
                    if numel(un_par) == 2

                        indd = find(NonezeroSig==1);
                        xr1 = xrt{indd(1)};
                        xr2 = xrt{indd(2)};
                        nr1 = nrt{indd(1)};
                        nr2 = nrt{indd(2)};
                        [xr1, xr2] = meshgrid(xr1, xr2);
                        [nr1, nr2] = meshgrid(nr1, nr2);
                        for coli = 1:numel(parameter)
                            if bin(coli) == 1
                                xr0 = xr1(:);
                            else
                                xr0 = xr2(:);
                            end
                            xr(:, indd(bin(coli))) = xr0;
                        end

                        nr = nr1(:).*nr2(:)/sum(nr1(:).*nr2(:));
                    else
                        error('Not implemented yet.')
                    end
                end
            end
        else
            xr = parameter; nr = 1;
        end
%         if numel(parameter)>1
%             xr = xr/parameter(1)*parameter;
%         else
%             xr = parameter;
%         end
        
        nr = nr(:);
    % =========================================

        if typeofinput < 3
            switch typeofinput
    %             case 1 % When hkls is input
    %                 QforPq = zeros(numel(hkls),1);
    %                 for qindx = 1:numel(hkls)
    %                     % Convert HKLs to q vectors
    %                     HKLv = hkls(qindx).HKLs;
    %                     zeroindx = find((HKLv(1,:)==0) & (HKLv(2,:)==0) & (HKLv(3,:)==0));
    %                     if ~isempty(zeroindx)
    %                         HKLv(:,zeroindx) = [];
    %                     end
    %                     h = HKLv(1,:);
    %                     k = HKLv(2,:);
    %                     l = HKLv(3,:);
    %                     qv = 2*pi*B*[h;k;l];
    %                     qx = qv(1, :)';
    %                     qy = qv(2, :)';
    %                     qz = qv(3, :)';
    % 
    %                     % Rotating particle or q.
    %                     q_transform; 
    %                     
    %                     % Compute Fq
    %                     Fq = calFq(fffunctionname, qxn, qyn, qzn, parameter, particles{pshape}, xr, nr);
    %                     if numel(particles)==1
    %                         Fqhkls(qindx, 1:numel(Fq)) = Fq';
    %                     else
    %                         Fqhkls(qindx, 1:size(Fq,1), pshape) = Fq';
    %                     end
    %                     
    %                     % Compute Pq if required.
    %                     %QforPq(qindx) = q(1);
    %                 end

                case 2 % when qx, qy, qz are the input.
                    q_transform
                    Fq = calFq(fffunctionname, qxn, qyn, qzn, parameter, particles{pshape}, xr, nr);
                    if numel(particles) ==1
                        Fqhkls = Fq;
                    else
                        Fqhkls(:, pshape) = Fq;
                    end
            end

        end
    end
    
    fname{pshape} = fffunctionname;
    if isfield(particles{pshape}, 'Rotmat')
        Rm{pshape} = particles{pshape}.Rotmat;
    else
        Rm{pshape} = eye(3);
    end

    ob{pshape} =  particles{pshape};
    switch fffunctionname
        case {'sphere', 'pdbFq', 'boxFq', 'rho2Fq'}
            n{pshape} = numPfrac(pshape);
            x{pshape} = parameter;
        otherwise
            n{pshape} = nr*numPfrac(pshape);
            x{pshape} = xr;
    end
    
end


%is_Fq_o_avgReq = 0;
Pqall_ = 0;
Fqmall_ = 0;
Pqmall_ = 0;
Fqall_ = 0;
Pq_ = 0;

if nargout >= 2
    % either compute <|F(q)|^2> or <F(q)>
    if nargout < 4
%        [Pq, Pqall] = calPqMultiparticle(fname, QforPq, x, particles, x, n, dim, Rm);
        [Pq_, Pqall_] = calPqMultiparticle(fname, QforPq, x, ob, x, n, dim, Rm);
        Fqmall_ = [];
%        Pqhkls = Pq_;
    elseif nargout >= 4  
        % <F(q)>_orient for only a q array for now.
%        [Pq, Pqall, Fqmall, Pqmall, Fqall] = calPqMultiparticle(fname, QforPq, x, particles, x, n, dim, Rm);
        if strcmp(fffunctionname, 'pdbFq')
            Pq_ = 1; Pqall_ = 1; Fqmall_ = 1; Pqmall_ = 1; Fqall_=1;
        else
            [Pq_, Pqall_, Fqmall_, Pqmall_, Fqall_] = calPqMultiparticle(fname, QforPq, x, ob, x, n, dim, Rm);
        end
%         Fqmall = zeros(numq, numel(fname));
%         for i=1:numel(fname)
%             Fqmall(:,i) = saxs_averageFq(QforPq, fname{i}, x{i}, n{i}, dim, Rm{i});
%         end
    end
end

isOverlapped = 0;
Pqall = Pqall_;
Fqmall = Fqmall_;
Pqmall = Pqmall_;
Fqall = Fqall_;
Pqhkls = Pq_;

for i=1:numel(Ptindex)  % this is to average particles occupying the same position.
    Noverlap = sum(Pindex==i);
    if Noverlap <=1
        continue
    end
    cprintf('red', 'There are %i particles occupying [%0.3f, %0.3f, %0.3f]\n', Noverlap, Ppos(i, :));
    cprintf('red', 'Their P(q) and F(q) are to be averaged.\n');
    isOverlapped = 1;
    ptInd2Avg = ptinput(Pindex==i, 2);
    Pqall0 = sum(Pqall_(:, ptInd2Avg), 2);
    Pqall(:,  ptInd2Avg) = repmat(Pqall0, 1, Noverlap);
    if nargout >=4
        Fqall0 = mean(Fqmall_(:, ptInd2Avg), 2);
        Fqmall(:,  ptInd2Avg) = repmat(Fqall0, 1, Noverlap);

        Pqall0 = mean(Pqmall_(:, ptInd2Avg), 2);
        Pqmall(:,  ptInd2Avg) = repmat(Pqall0, 1, Noverlap);

        Fqall0 = mean(Fqall_(:, ptInd2Avg), 2);
        Fqall(:,  ptInd2Avg) = repmat(Fqall0, 1, Noverlap);
    end
end
if isOverlapped
    Pqhkls = sum(Pqall, 2);
end

if exist('Pqhkls', 'var')
    Pqhkls = Pqhkls * Angstrom2Centimeter^2;
end
if exist('Pqall', 'var')
    Pqall = Pqall * Angstrom2Centimeter^2;
end
if exist('Pqmall', 'var')
    Pqmall = Pqmall * Angstrom2Centimeter^2;
end

if exist('Fqhkls', 'var')
    Fqhkls = Fqhkls * Angstrom2Centimeter;
end
if exist('Fqmall', 'var')
    Fqmall = Fqmall * Angstrom2Centimeter;
end
if exist('Fqall', 'var')
    Fqall = Fqall * Angstrom2Centimeter;
end

try
    assignin('base', 'Fqhkls', Fqhkls)
    assignin('base', 'Pqhkls', Pqhkls)
end
    function q_transform
            if isfield(particles{pshape}, 'Rotmat')
                Rmat = particles{pshape}.Rotmat;
                %Q' = inv(mat) * [qx; qy; qz];
                Q = [qx(:), qy(:), qz(:)] * Rmat;% for rotation matrix R^T = R^-1.
                qxn = Q(:,1);
                qyn = Q(:,2);
                qzn = Q(:,3);
            else
                qxn = qx;
                qyn = qy;
                qzn = qz;
            end
            q = sqrt(qxn.^2+qyn.^2+qzn.^2);
    end
end
    
function Fq = calFq(fffunctionname, qx, qy, qz, parameter, obj, xr, nr)
    % If obj has the field Fq, Fq will be interpolated.
    % This will be mostly used for sphere or randomly oriented particles.
    if strcmp(obj.shape, 'rho-spherical')
        Fq = rho2Fq(parameter(:,1), parameter(:,2), sqrt(qx.^2+qy.^2+qz.^2));
        return
    end
    rho = obj.rho;
    if isfield(obj, 'Fq')
        if ~isempty(obj.Fq)
            P = obj.Fq;
            Fq = interp1(P(:,1),P(:,2),sqrt(qx.^2+qy.^2+qz.^2));
            Fq = Fq*rho;
            return
        end
    end
        
    if size(xr, 1) == 1
        % Monodisperse ordered particles.
        if isempty(obj.Lshell)
            Fq = feval(fffunctionname, qx, qy, qz, parameter);
        else
            Fq = feval(fffunctionname, qx, qy, qz, parameter, 'shell thickness', obj.Lshell);
        end
    else
        % Polydisperse ordered particles.
        if obj.orientationfactor == 0
            Fq = zeros(numel(qx), 1);
            for isize=1:numel(nr)
                if isempty(obj.Lshell)
                    tempint = feval(fffunctionname, qx, qy, qz, parameter*xr(isize)/parameter(1));
                else
                    tempint = feval(fffunctionname, qx, qy, qz, parameter*xr(isize)/parameter(1), 'shell thickness', obj.Lshell);
                end
                tempint(isnan(tempint)) = 0;
                Fq = Fq + tempint*nr(isize);
            end
        else
            % Both size and orientation average.......
            if size(parameter, 2) > 1
                if size(parameter, 2) == size(xr, 2)
                    xr = xr';
                    nr = nr';
                end
            end
            Fq = saxs_averageFq(sqrt(qx.^2+qy.^2+qz.^2), fffunctionname, parameter*xr/parameter(1), nr, 1, obj.Rotmat);
        end
    end
    Fq = Fq*rho;
end

function [Pq, Fqm, Pqmono, Fq] = calPq(fffunctionname, q, parameter, obj, xr, nr, dim, Rm)
    % If obj has the field Pq, Pq will be interpolated.
    % This reduces computation time a lot.
    Angstrom2Centimeter = 1E-8;

    rho = obj.rho;
    Pq = [];
    Fqm = [];
    Pqmono = [];
    Fq = [];
    if isfield(obj, 'Pq')
        if ~isempty(obj.Pq)
            P = obj.Pq;
            Pq = interp1(P(:,1),P(:,2),q);
            Pq = Pq*rho^2;
            if nargout==1
                return
            end
        end
    end
    if isfield(obj, 'Fqmean')
        if ~isempty(obj.Fqmean)
            Fq = obj.Fqmean;
            Fqm = interp1(Fq(:,1),Fq(:,2),q);
            Fqm = Fqm*rho;
            return
        end
    end
    
%     if numel(xr) == 1
%         nr = 1;
%     end

    switch fffunctionname
        case 'rho2Fq'
            r_ = parameter(:,1);
            rho_ = parameter(:,2);
            Aq = rho2Fq(r_, rho_, q);
            Pq = abs(Aq).^2;
            Fqm = Aq;
            Fq = Aq;
            Pq = nr*Pq;
%            Fqm = sqrt(nr)*Fqm;
            %Fq = Fqm;
            %Pqmono = abs(Fq).^2;
        case 'boxFq'
            if (dim==1) | (dim==2)
                pdb_orient = 3;
            end
            if dim==3
                pdb_orient = 2;
            end
            if dim==4
                pdb_orient = 1;
            end
            if nargout == 1
                Pq = boxPq(q, parameter, pdb_orient);
            elseif nargout == 2
                [Pq, Fqm] = boxPq(q, parameter, pdb_orient);
            elseif nargout == 4
                [Pq, Fqm, ~, Fq] = boxPq(q, parameter, pdb_orient);
            end
            Pq = nr*Pq;
            %Fqm = sqrt(nr)*Fqm;
            %Pqmono = abs(Fq).^2;
        case 'pdbFq'
            if (dim==1) | (dim==2)
                pdb_orient = 3;
            end
            if dim==3
                pdb_orient = 2;
            end
            if dim==4
                pdb_orient = 1;
            end
            if nargout == 1
                Pq = pq_pdb(q, parameter, pdb_orient);
            else
                [Pq, Fqm] = pq_pdb(q, parameter, pdb_orient);
            end
            Pq = nr*Pq;
            %Fqm = sqrt(nr)*Fqm;
        case 'sphere'
            [Fq, Pq] = SphFormFactorforSq2(q, obj); 
            obj.radius_sig = 0;
            [Fqm, Pqmono] = SphFormFactorforSq2(q, obj); 
%             Fq = Fq / Angstrom2Centimeter;
%             Fqm = Fqm / Angstrom2Centimeter;
%             Pq = Pq / Angstrom2Centimeter^2;
%             Pqmono = Pqmono / Angstrom2Centimeter^2;
             return
        otherwise
            if nargout == 1
                if isempty(obj.Lshell)
                    Pq = saxs_average(q, fffunctionname, xr, nr, dim, Rm, obj.edgelength);
                else
                    Pq = saxs_average(q, fffunctionname, xr, nr, dim, Rm, obj.edgelength, 'shell thickness', obj.Lshell);
                end
            else
                if isempty(obj.Lshell)
                    [Pq, Fqm, Pqmono, Fq] = saxs_average(q, fffunctionname, xr, nr, dim, Rm, obj.edgelength);
                else
                    [Pq, Fqm, Pqmono, Fq] = saxs_average(q, fffunctionname, xr, nr, dim, Rm, obj.edgelength, 'shell thickness', obj.Lshell);
                end
                if numel(nr) == 1
                    Pq = Pq*nr;
                    Pqmono = Pqmono*nr;
                else
                %%%% New addition to avoid S(q --> inf) deviates from 1.
                %%%% Nov. 2017. BLee
%                    Pq = Pq*sum(nr);
%                   intead, n=nr/sum(nr) in saxs_average.m is captioned
%                   out.
                end
            end
    end
    Pq = Pq*rho^2;
    Pqmono = Pqmono*rho^2;
    Fqm = Fqm*rho;
    Fq = Fq*rho;
end

function [Pq, Pqall, Fqall, Pqm_all, Fq] = calPqMultiparticle(fffunctionname, q, parameter, obj, xr, nr, dim, Rm)

    Pq = zeros(numel(q(:,1)),1);
    Pqall = zeros(numel(q(:,1)),numel(fffunctionname));
    Fqall = zeros(numel(q(:,1)),numel(fffunctionname));
    Pqm_all = Fqall;
    Fq = Fqall;
    for i=1:numel(fffunctionname)
        Ptmp = []; Fqtmp = []; Ftmp = [];
        if nargout==2
            tmp = calPq(fffunctionname{i}, q, parameter{i}, obj{i}, xr{i}, nr{i}, dim, Rm{i});
        else
            [tmp, Ftmp, Ptmp, Fqtmp] = calPq(fffunctionname{i}, q, parameter{i}, obj{i}, xr{i}, nr{i}, dim, Rm{i});
        end
        
        if isempty(Ptmp)
            Ptmp = zeros(size(tmp));
        end
        if isempty(Fqtmp)
            Fqtmp = zeros(size(tmp));
        end
%         if ~strcmp(fffunctionname{i}, 'sphere')
%             tmp = tmp * Angstrom2Centimeter^2;
%             Ftmp = Ftmp * Angstrom2Centimeter;
%             Ptmp = Ptmp * Angstrom2Centimeter^2;
%             Fqtmp = Fqtmp * Angstrom2Centimeter;
%         end
        Pq = Pq + tmp;
        Pqall(:,i) = tmp(:);
        if nargout >= 3
            Fqall(:,i) = Ftmp(:);
            Pqm_all(:,i) = Ptmp(:);
            if ~isempty(Fqtmp)
                Fq(:,i) = Fqtmp(:);
            end
        end
    end
end
