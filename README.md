# Overview
This is a MATLAB package for curve fitting.

## Requirements
This program is written in MATLAB and requires a MATLAB license.
No additional MATLAB toolbox is needed.

## Installation
Download the code from the 'Code' pull-down menu, for example as a zip file.

## Credits
This package uses fminsearchcon.m, written by John D'Errico. For your convenience, it is downloaded and included.
See the details in the .m file.

# Usage

## How to start
1. Add the source, data, and download folders to your MATLAB path. If you are not familiar with this, take a look at the [directions from MathWorks](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html).

2. To run FitLee, type `>> startFitLee` at the MATLAB prompt.
3. Or, if you write your own fitting code, you can load it into FitLee, for example, `>> FitLee('myfitting_code.m')`.

## Fitting multiple curves
It is possible to fit multiple datasets with the same model.
See the example below and modify it accordingly.
```matlab
% Assuming you have already done a curve fitting using FitLee ...
% for example, the pseudo-Voigt function fitting with 4 peaks.
%   ex) FitLee('FitLee_multipseudovoigt', 4)
%   Choose the right q range and parameters. Then, fit.
%   You may have to adjust LB and UB.
%   Let's say your fit looks good.

% Now you want to load new data, fit, and save the result:

% load new data:
% below, you should replace mydata.txt with your filename.
a = load('mydata.txt');

% intensity value
ydata = a(:,2);

% put the intensity to the FitLee
set(fit.datahandle, 'ydata', ydata);
fit.yd = ydata;
fit.data{1} = [a(:,1), a(:,2)];

% if needed change the Lower boundary (LB) or Upper boundary (UB)
while 1
    istherered = false;
    for i=1:numel(fit.param)
        LB = ['edit_LB', num2str(i)];
        UB = ['edit_UB', num2str(i)];
        LB_val = eval(fit.handles.(LB).String);
        UB_val = eval(fit.handles.(UB).String);
        % if color of LB%i is red
        if all(fit.handles.(LB).ForegroundColor == [1,0,0]) 
            % decrease LB by 50%
            sn = sign(LB_val);
            if sn>0
                fit.handles.(LB).String = LB_val*0.5;
            else
                fit.handles.(LB).String = sn*(abs(LB_val)*1.5);
            end
            istherered = true;
        end
        % if color of UB%i is red
        if all(fit.handles.(UB).ForegroundColor == [1,0,0]) 
            % increase the UB by 50%.
            sn = sign(UB_val);
            if sn>0
                fit.handles.(UB).String = UB_val*1.5;
            else
                fit.handles.(UB).String = sn*(abs(UB_val)*0.5);
            end
            istherered = true;
        end
    end
    
    % execute fitting
    feval(fit.handles.pushbutton_Fit.Callback)

    % when fit value is within the limit, exit the loop
    if ~istherered
        break
    end
end

% save the result
% below, you should replace test.txt with your result filename.
fit.handles.edit_savefilename.String = "test.txt";
feval(fit.handles.pushbutton_Save.Callback)
```

## How to write your own fitting code
- The fitting code consists of four blocks. For example, have a look at models>FitLee_schultzsphere2.m.

1. Block I:
- The text block that will appear when the "Reference" button is pressed.

```matlab
function [out, report] = FitLee_schultzsphere2(varargin)
FitLee_helpstr = {'Schultz polydisperse sphere fit in absolute unit. ' ,...
'$I(q) = fn_0\cdot(Sq\cdotP(q; r_0, \sigma_0) + N_{ratio}\cdot\delta_{\rho1}^2*P(q; r_1, \sigma_1)) + Ib$',...
...
'Ref: ',...
'    1. Kwon et al. Nature Materials, 2015, 14, 215–223. ',...
'    2. Wang et al. J. Phys. Chem. C. 2013. 117(44), 22627. '};
```

2. Block II:
- This block contains the fit parameters and initializes the fitting GUI.
- Do not change this block.

```matlab
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

```

- The following block lists your fitting parameters.
- Edit it as needed.

```matlab
%% initialize fit parameter bestP
if isini
    Nf = p;         % Do not change this line.
    bestP = [];     % Do not change this line.
...
    bestP.poly4 = 0;                    % List your parameters as a structure. The field name will appear in the GUI.
    assignin('base', 'bestP', bestP);   % Do not change this line.
    out = bestP;                        % Do not change this line.
    return
end

```

3. Block III:
- Your own fitting code.

```matlab
%% fitting code .......................
if iscell(q)
    if numel(q) > 1
        error('your own error message')
    end
    q = q{1};
end

r_e = 2.818E-5; % Angstrom
....
% out is the output. It can have multiple columns.
out = pnumberfraction*r_e^2*Angstrom2Centimeter^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;
if isnan(out)
    out = ones(size(out));
end
```

4. Block IV:
- This block runs when the "Print aux result" button is pressed.
- Use this block to compute additional information, for example to draw the model or size distribution as below:

```matlab
if nargout == 2 % Do not change this line.
    x = 0:1:(max(p.r0, p.r1)+max(p.sig0, p.sig1)*10);
    ....
    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(r)')
    fprintf('Statistical information of particle 0 ======================================\n');
    fprintf('Number-mean radius of a single particle : %0.3e %c.\n', p.r0, char(197));
    ...
    fprintf('Weight concentration (g/mL) can be obtained by multiplying your particles'' density (g/mL) by fn0.\n');
    fprintf('==============================================================\n');
    
    report = '';    % Do not change this line.
end
```
