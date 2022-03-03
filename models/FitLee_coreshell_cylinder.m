function [out, report] = FitLee_coreshell_cylinder(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse coreshell cylinder radial direction fit. ' ,...
    'Data needs to be plotted as q vs Iq*q',...
'I(q) = I0*P(q; r0, sig0)*S(q) + background',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif numel(varargin) == 1
    p = varargin{1};
    isini = 1;
    if ischar(p)
        out = FitLee_helpstr;
        return
    end
elseif numel(varargin) == 0
    p = 1;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    bestP = [];
    bestP.I0 = 4.6095e-13;
    bestP.coreR = 25;
    bestP.sigma = 0.1;
    bestP.shellthick = 10;
    bestP.core_rho = 0;
    bestP.shell_rho = 0.39;
    bestP.solvent_rho = 0.3344;
    bestP.H = 200;
    bestP.sigH = 0.1;

    bestP.I0_2 = 1;
    bestP.s_2 = 2;
    bestP.Rg_2 = 300;
    bestP.P_2 = 4;
    
    bestP.D = 100;
    bestP.vf = 0.01;
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q)
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

% radial direction
q = q(:);
radius = [p.coreR, p.coreR + p.shellthick];
edensity = [p.core_rho, p.shell_rho, p.solvent_rho];
Iq = saxscylinder_CS(q, radius, edensity, p.sigma);

% length direction
if p.H == 0
    ILq = 1;
else
    numpnt = 55;
    x = linspace((p.H/10).^(1/1.5), (p.H + p.H*p.sigH*9).^(1/1.5), numpnt);x = x(:)';
    x = x.^1.5;
    nr = schultzdist(x, p.H, p.H*p.sigH);
    nr = nr(:)';
    ILq = sinc(q*x/2).^2.*repmat((sqrt(nr).*x).^2, length(q), 1);
    ILq = sum(ILq, 2);
end

Pq_1 = guinierporodmodel(q, p.I0_2, p.P_2, p.Rg_2, p.s_2);
Sq1 = strfactor2(q, p.D, p.vf)+Pq_1;
%Sq = p.powI*q.^p.PorodExp + Sq1;

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = p.I0*Iq.*Sq1.*ILq + back;
if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    maxR = p.coreR+p.shellthick;
    x = linspace(maxR/10, maxR + maxR*sigma*9, 100);
    nr = schultzdist(x, maxR, maxR*p.sigma);

    
    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(R)')
    report = '';
end


