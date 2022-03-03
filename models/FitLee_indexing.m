function [out, report] = FitLee_indexing(varargin)
    % Fit peaks found from indexing.m
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Fitting intensities of peaks found from indexing.m',...
    'Use for Lorentz corrected data', ...
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
    if ishandle(p)
        guiadd(p);
        return
    end
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
    if exist('pseudovoigt', 'file') ~= 2
        error('Need "pseudovoigt.m". Download and keep in with FitLee_multipseudovoigt.m.')
    end

    Nf = p;
    bestP = [];
    bestP.scaleI0 = 1;
    bestP.dsize = 2000;
    bestP.mstrain = 0.01;
    bestP.gw = 0.0002;
    bestP.amorph_a = 1;
    bestP.amorph_x0 = 0.02;
    bestP.amorph_sigG = 0.01;
    bestP.amorph_sigL = 0.01;
    bestP.amorph2_a = 1;
    bestP.amorph2_x0 = 0.04;
    bestP.amorph2_sigG = 0.01;
    bestP.amorph2_sigL = 0.01;
    bestP.bscale = 1E-5;
    bestP.exponent = 4;
    bestP.back_a = 0.1;
    bestP.back_ax = 0.1;
    bestP.back_ax2 = 0.01;
    
    bestP.string = 'FitLee_indexing(figH)';
    peak = evalin('base', 'indexing_hkls');
    for i=1:numel(peak(:,1))
        bestP.(sprintf('no_display_amp%i', i)) = peak(i, 2)/10;
        bestP.(sprintf('no_display_back%i', i)) = peak(i, 2)/10000;
    end
    
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q)
    if numel(q) > 1
        error('FitLee_indexing.m is for fitting a set of data, for now')
    end
    q = q{1};
end


q = q(:);
peak = evalin('base', 'indexing_hkls');
try
    backdata = evalin('base', 'background_drawn');
catch
    b = findobj('color', 'r', 'tag', 'back');
    if ~isempty(b)
        xd = get(b, 'xdata');
        yd = get(b, 'ydata');
        backdata = [xd(:), yd(:)];
        assignin('base', 'background_drawn', backdata);
        delete(b);
    else
        backdata = [];
    end
end
pos = zeros(size(peak, 1), 1);
amp = pos;
backamp = pos;
%ud = get(gcbf, 'userdata');
for i=1:numel(peak(:,1))
    pos(i) = peak(i, 1);
    amp(i) = p.(sprintf('no_display_amp%i', i));
    backamp(i) = p.(sprintf('no_display_back%i', i));
    peak(i, 2) = amp(i);
end

%backsp = zeros(size(q));
posn = [0;pos;q(end)];
backamp = [backamp(1);backamp;backamp(end)];
backamp = smooth(backamp, 3);
%for i=1:numel(pos)
backsp = interp1(posn,backamp, q);
%    backsp = ppval(PP, q);
%end

Iq = powderpeaks(q, pos, amp, p);
Iq = p.scaleI0*Iq;
%back = p.bscale*(1-exp(-p.constback^2*q.^2))./q.^(p.exponent);
back = p.bscale./q.^(p.exponent) + p.back_a + p.back_ax*q+p.back_ax2*q.^2;
%backPeak = gaussa(q, [p.amorph_a, p.amorph_x0, p.amorph_sig, 0]);
backPeak = pseudovoigt(q, [p.amorph_a, p.amorph_x0, p.amorph_sigG, p.amorph_sigL]);
backPeak2 = pseudovoigt(q, [p.amorph2_a, p.amorph2_x0, p.amorph2_sigG, p.amorph2_sigL]);
if ~isempty(backdata)
    back = back+interp1(backdata(:,1), backdata(:,2), q);
end
back = back + backPeak+backPeak2 + backsp;
out = [Iq(:) + back, back];

if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    report = savehkls(peak);
end

function guiadd(figH)
    uph = findobj(figH, 'type', 'uipanel');
    pos = get(figH, 'position');
    hFigHeight = pos(end);
    %[1, hFigHeight-110,206,20]
    popup = uicontrol(...
          'parent', uph,...
          'horizontalalignment'       , 'Left', ...
          'Style', 'pushbutton',...
          'String', 'Draw Background',...
          'Position', [50, hFigHeight-230, 150, 30],...
          'Callback', {@drawback, figH});
      
function drawback(varargin)
    edfig = findobj(varargin{3}, 'tag', 'edit_Figure');
    figh = str2double(get(edfig, 'string'));
    drawContinousBack('start', figh)
    
function report = savehkls(varargin)
    peak = varargin{1};
    [filename, pathname] = uiputfile({'*.txt'; '*.dat'; '*.inflip'; '*.*'}, 'Save hkls As');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
        return
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    filename = fullfile(pathname, filesep, filename);
    [~, inputName, exttype] = fileparts(filename);
    fid = fopen(filename, 'w');
    if contains_replace(exttype, 'inflip')
        sg = evalin('base', 'indexing_sg');
        sginfo = sgroup(sg);
        hkls = evalin('base', 'indexing_hklstr');
        cellinfo = evalin('base', 'indexing_cellinfo');
        
        fprintf(fid, 'title  %s,  Spacegroup simulation\n', inputName);
        fprintf(fid, '\n');
        fprintf(fid, '# crystallographic data\n');
        fprintf(fid, '#----------------------\n');
        fprintf(fid, 'cell  %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n', cellinfo.A, cellinfo.B, cellinfo.C, cellinfo.alpha, cellinfo.beta, cellinfo.gamma);
        fprintf(fid, '\n');
        fprintf(fid, '# spacegroup %s\n', sg);
        fprintf(fid, 'centers\n');
        for i=1:sginfo.NoLatticeCenteringVector
            cv = sginfo.LatticeCenteringVector(:, i);
            fprintf(fid, '  %0.6f %0.6f %0.6f\n', cv(1), cv(2), cv(3));
        end
        fprintf(fid, 'endcenters\n');
        fprintf(fid, '\n');
        fprintf(fid, '\n');
        fprintf(fid, '# Grid definition for density maps\n');
        fprintf(fid, '#---------------------------------\n');
        fprintf(fid, 'dimension  3\n');
        fprintf(fid, 'voxel   100  100  100\n');
        fprintf(fid, '\n');
        fprintf(fid, '\n');
        fprintf(fid, '# symmetry\n');
        fprintf(fid, 'symmetry\n');
        for i=1:numel(sginfo.Symmetry)
            fprintf(fid, '%s\n', sginfo.Symmetry{i});
        end
        fprintf(fid, 'endsymmetry\n');
        fprintf(fid, '\n');
        fprintf(fid, '\n');
        fprintf(fid, 'outputfile %s\n', fullfile(inputName, '.xplor'));
        fprintf(fid, '# Keywords for charge flipping\n');
        fprintf(fid, '#-----------------------------\n');
        fprintf(fid, 'delta  AUTO\n');
        fprintf(fid, 'weakratio  0.0\n');
        fprintf(fid, 'Biso  0.0005\n');
        fprintf(fid, 'randomseed AUTO\n');
        fprintf(fid, 'searchsymmetry average\n');
        fprintf(fid, 'derivesymmetry yes\n');
        fprintf(fid, 'finevoxel 5 angstrom\n');
        fprintf(fid, 'coverage no\n');
        fprintf(fid, '# input data\n');
        fprintf(fid, '#-----------\n');
        fprintf(fid, '\n');
        fprintf(fid, 'dataformat intensity group\n');
        fprintf(fid, '\n');
        fprintf(fid, 'fbegin\n');
        fprintf(fid, '#   h    k    l      |Fobs|^2    group index\n');
        % Lorentz Correction.
        %peak(:,2) = peak(:,2).*peak(:,1).^2;
        % Multiplicity Correction..
        hklinfo.DMIN = 2*pi/max(peak(:,1));
        [~, HKLs] = indexref(cellinfo, sginfo, hklinfo);
        
        maxInt = max(peak(:,2));
        Npeak = numel(peak(:,2));
        cnt = 1;
        for i=1:numel(hkls)
            hkl0 = hkls(i).hkl;
            if i>Npeak
                break
            end
            for k=1:size(hkl0, 1)
                for m=1:size(HKLs(cnt).HKLs, 2)
                    hkl = HKLs(cnt).HKLs(:, m);
                    fprintf(fid, '%i %i %i %0.6f %i\n', hkl(1), hkl(2), hkl(3), peak(i, 2)/maxInt*100, i);
                end
                cnt = cnt+1;
            end
        end
        fprintf(fid, 'endf\n');
    else
        for i=1:size(peak, 1)
            fprintf(fid, '%0.8e      %0.8e\n', peak(i, 1), peak(i, 2));
        end
    end
    fclose(fid);
    report = 'Done';