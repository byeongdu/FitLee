# Overview
This is a matlab package for curve fitting. 

# Requirement
This is a program written in matlab, requiring a matlab license.
No additional matlab toolbox is needed.

# Installation
Download codes from the 'code' pulldown menu, for example as a zip. 

# How to start
1. Define the source, data, and download folders as your matlab path. If not familiar with this, take a look at the direction from Mathworks. [https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html#:~:text=Change%20Folders%20on%20Search%20Path%20Interactively,-Use%20the%20Set&text=On%20the%20Home%20tab%2C%20in,folders%20to%20MATLAB%20search%20path.]

2. To run FitLee, on matlab prompt >> startFitLee. 
3. Or, if you write your own fitting codes, you can load into FitLee, for example, >> FitLee('myfitting_code.m')

# Credits
fminsearchcon.m is written by John D'Errico. See the detail from the m file.

## How to write your own fitting code.
- The fitting code consists of four blocks. For example, have a look at models>FitLee_schultzsphere2.m.
1. Block I: the text block that will appear when "Reference" button is pressed.

function [out, report] = FitLee_schultzsphere2(varargin)
FitLee_helpstr = {'Schultz polydisperse sphere fit in absolute unit. ' ,...
'$I(q) = fn_0\cdot(Sq\cdotP(q; r_0, \sigma_0) + N_{ratio}\cdot\delta_{\rho1}^2*P(q; r_1, \sigma_1)) + Ib$',...
...
'Ref: ',...
'    1. Kwon et al. Nature Materials, 2015, 14, 215–223. ',...
'    2. Wang et al. J. Phys. Chem. C. 2013. 117(44), 22627. '};

2. Block II: the block contains fit parameters and intialize the fitting GUI.
- Do not change this block...
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
- Edit the block as needed.
if isini
    Nf = p;         # Do not change this line.
    bestP = [];     # Do not change this line.
...
    bestP.poly4 = 0;                    # List your parameters as a structure. The field name will appear on GUI.
    assignin('base', 'bestP', bestP);   # Do not change this line.
    out = bestP;                        # Do not change this line.
    return
end

%% fitting codes .......................
if iscell(q)
    if numel(q) > 1
        error('your own error message')
    end
    q = q{1};
end

3. Block III : fitting code
r_e = 2.818E-5; % Angstrom
....
out = pnumberfraction*r_e^2*Angstrom2Centimeter^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;
- out is the output. It can have multiple columns. 
if isnan(out)
    out = ones(size(out));
end

4. Block IV: this block will run when "Print aux result" button.
- Use this block for compute additional information, for example to draw the model or size distribution as below:
if nargout == 2 # Do not change this line
    x = 0:1:(max(p.r0, p.r1)+max(p.sig0, p.sig1)*10);
    ....
    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(r)')
    fprintf('Statistical information of the particle 0 ======================================\n');
    fprintf('Number-mean radius of a single particle : %0.3e %c.\n', p.r0, char(197));
    ...
    fprintf('Weight concentration (g/mL) can be obtained by multiplying your particles'' density (g/mL) to the fn0.\n');
    fprintf('==============================================================\n');
    
    report = '';    # do not change this line.
end
