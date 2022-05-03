function figH = FitLee(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.         %
%                                                                         %
% Copyright 2016. UChicago Argonne, LLC. This software was produced       %
% under U.S. Government contract DE-AC02-06CH11357 for Argonne National   %
% Laboratory (ANL), which is operated by UChicago Argonne, LLC for the    %
% U.S. Department of Energy. The U.S. Government has rights to use,       %
% reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR    %
% UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR        %
% ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is     %
% modified to produce derivative works, such modified software should     %
% be clearly marked, so as not to confuse it with the version available   %
% from ANL.                                                               %
%                                                                         %
% Additionally, redistribution and use in source and binary forms, with   %
% or without modification, are permitted provided that the following      %
% conditions are met:                                                     %
%                                                                         %
%     * Redistributions of source code must retain the above copyright    %
%       notice, this list of conditions and the following disclaimer.     %
%                                                                         %
%     * Redistributions in binary form must reproduce the above copyright %
%       notice, this list of conditions and the following disclaimer in   %
%       the documentation and/or other materials provided with the        %
%       distribution.                                                     %
%                                                                         %
%     * Neither the name of UChicago Argonne, LLC, Argonne National       %
%       Laboratory, ANL, the U.S. Government, nor the names of its        %
%       contributors may be used to endorse or promote products derived   %
%       from this software without specific prior written permission.     %
%                                                                         %
% THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS     %
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT       %
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS       %
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago     %
% Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,        %
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,    %
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;        %
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER        %
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      %
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN       %
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE         %
% POSSIBILITY OF SUCH DAMAGE.                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Usage:
% FitLee ('function name', parameter)
% FitLee ('function name', action)
%   for actions, type the callback function name.
%   for instance, FitLee('multivoight', 'pushbutton_Fit_callback')
% User can define 'pre_userfit' function handle.
%   This will be run before fit and simulation.
% Example:
% FitLee('FitLee_schultzsphere', 'pushbutton_Fit_Callback')
% FitLee('FitLee_polyhedra', 'pushbutton_Fit_Callback')
%
% If you like to add uicontrol to FitLee, you can pass FitLee's handle to
% particular macro you will use and then let the macro add uicontrols to
% FitLee. see FitLee_polyhedra for this. It adds uicontrol to pick a
% particle shape.

if exist('fminsearchcon', 'file') ~= 2
    error('Need "fminsearchcon.m". Download and keep it in FitLee folder.')
end

fitfunctionname = 'FitLee_multivoight';% default function;

if numel(varargin) >= 1
    fitfunctionname = varargin{1};
end
fitFunHandle =  str2func(fitfunctionname);
qNewFitLee = 0;
if numel(varargin) < 2
    qNewFitLee = 1;
    param = [];
else
    if isnumeric(varargin{2})
        qNewFitLee = 1;
        param = varargin{2};
    end
end
if qNewFitLee
    if isempty(param)
        p = fitFunHandle();
    else
        p = fitFunHandle(param);
    end

    figH = FitLee(fitfunctionname, p);
    % new addition =========================
    handles = guihandles(figH);
    %fit = evalin('base', 'fit');
    fit = [];
    fit.dattag = '__data';
    fit.fittag = '__fit';
    fit.simtag = '__sim';
    fit.handles = handles;
    setappdata(figH, 'Fit', fit)
    assignin('base', 'fit', fit)
    % ======================================
    return
end
call_obj = gcbf;
if ~isempty(call_obj)
    if strcmp(call_obj.Name, 'indexing')
        call_obj = [];
    end
end

% if numel(varargin) == 2 now proceed to make fit GUI.
% if isstruct(varargin{2})
%     bestP = varargin{2};
%     figH = initFigure;
%     get(figH, 'tag')
% else
%     figH = findall(0, 'type', 'figure', 'name', ['FitLee for ', fitfunctionname]);
%     if ~isempty(varargin{2})
%         eval(varargin{2});
%     end
% end
if ~isstruct(varargin{2})
%    figH = findall(0, 'type', 'figure', 'name', ['FitLee for ', fitfunctionname]);    
    figH = findall(0, 'type', 'figure', 'name', fitfunctionname);

    if ischar(varargin{2}) % when varargin{2} is a command such as 'xx_callback'
        eval(varargin{2});
    end
    return
end

bestP = varargin{2};
figH = initFigure(bestP, fitfunctionname);

% This code to pass the handle, figH, to Fit macro.
% bestP.string should be something like: 'FitLee_polyhedra(figH)'
if isfield(bestP, 'string')
    eval(bestP.string);
    bestP = rmfield(bestP, 'string');
end

handles = guihandles(figH);
try
    fit = evalin('base', 'fit');
catch
    fit = [];
    fit.dattag = '__data';
    fit.fittag = '__fit';
    fit.simtag = '__sim';
    fit.handles = handles;    
    setappdata(figH, 'Fit', fit)
    assignin('base', 'fit', fit);
end

% Define pre_fitFunHandle...
% This function will be run before simulation and fitting...
% This is added for the case requiring pre processing such as calculate
% Rmat for size distribution function fitting..
%if numel(varargin) == 3;
%    pre_fitfunctionname = varargin{3};
%    pre_fitFunHandle =  str2func(pre_fitfunctionname);
%else
%    pre_fitFunHandle = [];
%end
  
    function figH = initFigure(bestP, fitfunctionname)
        hFigHeight_default = 5*20 + 100;
        pnames = fieldnames(bestP);
        is_nodisplay = contains_replace(pnames, 'no_display');
        Npnames = numel(pnames)-sum(is_nodisplay);
        posScreen   = get(0,'screenSize');
        hFigWidth   = 650;
        hFigHeight  = Npnames*20 + 100;
        if hFigHeight<hFigHeight_default
            hFigHeight = hFigHeight_default;
        end
        hFigPos     = [...
            posScreen(3)/2-hFigWidth/2,...
            posScreen(4)/2-hFigHeight/2,...
            hFigWidth,hFigHeight];

%          'name'                        , ['FitLee for ', fitfunctionname], ...
        figH = figure(...
          'position'                    , hFigPos,...
          'visible'                       , 'on', ...
          'units'                         , 'pixel', ...
          'busyaction'                    , 'queue', ...
          'doublebuffer'                  , 'on', ...
          'handlevisibility'              , 'callback', ...
          'IntegerHandle'                 , 'off',...
          'interruptible'                 , 'on', ...
          'menubar'                       , 'none', ...
          'numbertitle'                   , 'off', ...
          'resize'                        , 'off', ...
          'name'                        , fitfunctionname, ...
          'tag'                           , 'FitLee', ...
          'toolbar'                       , 'none', ...
          'defaultaxesunits'              , 'pixels', ...
          'defaulttextfontunits'          , 'pixels', ...
          'defaulttextfontname'           , 'Verdana', ...
          'defaulttextfontsize'           , 12, ...
          'defaultuicontrolunits'         , 'pixels', ...
          'defaultuicontrolfontunits'     , 'pixels', ...
          'defaultuicontrolfontsize'      , 10, ...
          'defaultuicontrolfontname'      , 'Verdana', ...
          'defaultuicontrolinterruptible' , 'off');

        % Data loading and selection control
        uph = uipanel(...
          'units'                     , 'pixels', ...
          'parent'                    , figH, ...
          'BackgroundColor'             , [0, 0.5, 0.5],...
          'ForegroundColor'             ,[0, 0, 1],...
          'Borderwidth'                 ,2,...
          'Units'                       , 'pixels', ...
          'position'                    ,   [0, 0, 220, hFigHeight],...
          'bordertype'                , 'beveledin', ...
          'tag'                       , 'versionPanel');

        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [0,0.5,0.5], ...
          'foregroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'Left', ...
          'parent'                    , uph, ...
          'string'                  ,   'Load data either from disk or Figure',...
          'position'                  ,[1, hFigHeight-25,206,20], ...
          'tag'                       , 'text_dir1');
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [0,0.5,0.5], ...
          'foregroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'Left', ...
          'parent'                    , uph, ...
          'string'                  ,   'ROI: for each curve or a set for all',...
          'position'                  ,[1, hFigHeight-110,206,20], ...
          'tag'                       , 'text_dir2');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   'Load data',...
          'callback'                , @pb_loaddata_Callback,...
          'position'                  ,[10, hFigHeight-50,100,25], ...
          'tag'                       , 'pb_loaddata');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   'Plot options',...
          'callback'                , @pb_plotstyle_Callback,...
          'position'                  ,[110, hFigHeight-50,90,25], ...
          'tag'                       , 'pb_plotstyle');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   'Fit region read from cursors',...
          'callback'                , @pushbutton_qselect_Callback,...
          'position'                  ,[10, hFigHeight-130,200,20], ...
          'enable'                    , 'off',...
          'tag'                       , 'pushbutton_qselect');

        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '<',...
          'callback'                , @pb_LL_Callback,...
          'position'                  ,[180, hFigHeight-50-110,15,20], ...
          'tag'                       , 'pb_LL');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '>',...
          'callback'                , @pb_LU_Callback,...
          'position'                  ,[195, hFigHeight-50-110,15,20], ...
          'tag'                       , 'pb_LU');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '<',...
          'callback'                , @pb_RL_Callback,...
          'position'                  ,[180, hFigHeight-50-130,15,20], ...
          'tag'                       , 'pb_RL');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '>',...
          'callback'                , @pb_RU_Callback,...
          'position'                  ,[195, hFigHeight-50-130,15,20], ...
          'tag'                       , 'pb_RU');

        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '',...
          'CreateFcn'                , @edit_Figure_CreateFcn,...
          'callback'                , @edit_Figure_Callback,...
          'position'                  ,[100, hFigHeight-50-30,100,25], ...
          'tag'                       , 'edit_Figure');
        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '',...
          'callback'                , @edit_Lq_Callback,...
          'position'                  ,[65, hFigHeight-50-110,60,20], ...
          'tag'                       , 'edit_Lq');
        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '',...
          'callback'                , @edit_Rq_Callback,...
          'position'                  ,[65, hFigHeight-50-130,60,20], ...
          'tag'                       , 'edit_Rq');
        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '',...
          'callback'                , @edit_Lindex_Callback,...
          'position'                  ,[135, hFigHeight-50-110,40,20], ...
          'tag'                       , 'edit_Lindex');
        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '',...
          'callback'                , @edit_Rindex_Callback,...
          'position'                  ,[135, hFigHeight-50-130,40,20], ...
          'tag'                       , 'edit_Rindex');

        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [0,0.5,0.5], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   'Pick Figure',...
          'position'                  ,[10, hFigHeight-50-30,80,20], ...
          'tag'                       , 'text0');
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'right', ...
          'parent'                    , uph, ...
          'string'                  ,   'q min',...
          'position'                  ,[20, hFigHeight-50-110,40,20], ...
          'tag'                       , 'text1');
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'right', ...
          'parent'                    , uph, ...
          'string'                  ,   'q max',...
          'position'                  ,[20, hFigHeight-50-130,40,20], ...
          'tag'                       , 'text2');
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '@',...
          'position'                  ,[125, hFigHeight-50-110,10,20], ...
          'tag'                       , 'text3');
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , uph, ...
          'string'                  ,   '@',...
          'position'                  ,[125, hFigHeight-50-130,10,20], ...
          'tag'                       , 'text4');



        % Fitting control
        baseHpos = 200; % horizontal position where the controls for parameters start.
        baseVpos = 70; % vertical position where the controls for parameters start.
        % baseVpos is essentially the space that the following 5 controls take.
        Fitbuttonwidth = 120;
        Hbuttongap = 20;
        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'right', ...
          'parent'                    , figH, ...
          'position'                  ,[baseHpos+20, 20,Fitbuttonwidth,25], ...
          'tag'                       , 'edit_savefilename');

        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , figH, ...
          'string'                  ,   'Load Saved Result',...
          'callback'                , @pushbutton_Load_Callback,...
          'position'                  ,[baseHpos+Fitbuttonwidth+Hbuttongap,10,Fitbuttonwidth,25], ...
          'tag'                       , 'pushbutton_Load');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , figH, ...
          'string'                  ,   'Save Results',...
          'callback'                , @pushbutton_Save_Callback, ...
          'position'                  ,[baseHpos+Fitbuttonwidth+Hbuttongap, 35,Fitbuttonwidth,25], ...
          'tag'                       , 'pushbutton_Save');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , figH, ...
          'string'                  ,   'Fit',...
          'callback'                , @pushbutton_Fit_Callback,...
          'position'                  ,[baseHpos+Fitbuttonwidth+Hbuttongap, 60,Fitbuttonwidth,25], ...
          'tag'                       , 'pushbutton_Fit');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , figH, ...
          'string'                  ,   'Fit Options',...
          'callback'                , @pushbutton_setfitoptions,...
          'position'                  ,[baseHpos+Fitbuttonwidth*2+Hbuttongap, 60,Fitbuttonwidth,25], ...
          'tag'                       , 'pushbutton_setfitoptions');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , figH, ...
          'string'                  ,   'Print aux result',...
          'callback'                , @pushbutton_showAuxresult,...
          'position'                  ,[baseHpos+Fitbuttonwidth*2+Hbuttongap, 10,Fitbuttonwidth,25], ...
          'tag'                       , 'pushbutton_showAuxresult');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , figH, ...
          'string'                  ,   'Simulate',...
          'callback'                , @pushbutton_simulate_Callback,...
          'position'                  ,[baseHpos+20, 60,Fitbuttonwidth,25], ...
          'tag'                       , 'pushbutton_simulate');
        uicontrol(...
          'style'                     , 'pushbutton', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'center', ...
          'parent'                    , figH, ...
          'string'                  ,   'Reference',...
          'callback'                , @pushbutton_reference_Callback,...
          'position'                  ,[baseHpos+Fitbuttonwidth*2+Hbuttongap, 35,Fitbuttonwidth,25], ...
          'tag'                       , 'pushbutton_reference');

        if sum(is_nodisplay) > 0
            c = pnames(is_nodisplay);
            f = struct2cell(bestP);
            f = f(is_nodisplay);
            ud = cell2struct(f, c);
            set(gcf, 'userdata', ud)
        end
        for i=1:numel(pnames)
            position(1) = baseHpos + 20;
            position(2) = 20*(i-1) + 25 + baseVpos;
            label = pnames{i};
            if ~contains_replace(label, 'option')
                if ~contains_replace(pnames{i}, 'no_display') && ~strcmp(pnames{i}, 'string')
                    generate_uifit(position, label, bestP.(pnames{i}), i, figH);
                end
            end
        end
    end

    function generate_uifit(position, label, value, indx, figH)
        pLabel.h = 20;
        pLabN.h = pLabel.h;
        pVal.h = pLabel.h;
        pFit.h = pLabel.h;
        pLB.h = pLabel.h;
        pUB.h = pLabel.h;
        wgap = 5;
        pLabN.w = 20;
        pLabel.w = 100;
        pVal.w = 80;
        pFit.w = 40;
        pLB.w = 60;
        pUB.w = 60;
        if abs(value(1)) < 0.001
            editvalue = sprintf('%0.3e', value(1));
        else
            editvalue = sprintf('%0.3f', value(1));
        end
        
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'left', ...
          'string'                  , num2str(indx), ...
          'parent'                    , figH, ...
          'position'                  ,[position(1), position(2),pLabN.w,pLabN.h], ...
          'tag'                       , ['text_pN', num2str(indx)]);

        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'left', ...
          'string'                  , label, ...
          'parent'                    , figH, ...
          'position'                  ,[position(1)+pLabN.w, position(2), pLabel.w,pLabel.h], ...
          'tag'                       , ['text_p', num2str(indx)]);

        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'left', ...
          'string'                  , editvalue, ...
          'parent'                    , figH, ...
          'position'                  ,[position(1)+pLabN.w+pLabel.w+wgap, position(2),pVal.w,pVal.h], ...
          'tag'                       , ['edit_p', num2str(indx)]);

        uicontrol(...
          'style'                     , 'checkbox', ...
          'backgroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'left', ...
          'string'                  , 'Fit', ...
          'parent'                    , figH, ...
          'position'                  ,[position(1)+pLabN.w+pLabel.w+pVal.w + wgap*2, position(2),pFit.w,pFit.h], ...
          'tag'                       , ['checkbox_Fit', num2str(indx)]);


        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'callback'                  , @edit_LB_callback,...
          'horizontalalignment'       , 'left', ...
          'parent'                    , figH, ...
          'position'                  ,[position(1)+pLabN.w+pLabel.w+pVal.w + pFit.w + wgap*3, position(2),pLB.w,pLB.h], ...
          'tag'                       , ['edit_LB', num2str(indx)]);

        uicontrol(...
          'style'                     , 'edit', ...
          'backgroundcolor'           , [1,1,1], ...
          'callback'                  , @edit_UB_callback,...
          'horizontalalignment'       , 'left', ...
          'parent'                    , figH, ...
          'position'                  ,[position(1)+pLabN.w+pLabel.w+pVal.w + pFit.w + pLB.w + wgap*4, position(2),pUB.w,pUB.h], ...
          'tag'                       , ['edit_UB', num2str(indx)]);
    end
    function setplot(LineHandle, data)
        for i=1:numel(LineHandle)
            if iscell(data)
                dt = data{i};
            else
                dt = data;
            end
            set(LineHandle(i), 'ydata', real(dt));
        end
    end
    function pushbutton_setfitoptions(varargin)
        fit = evalin('base', 'fit');
        if isfield(fit, 'options')
            options = fit.options;
        else
            options = optimset('fminsearch');
            options = optimset(options, 'TolX',0.1E-6);
    %        options = optimset(options, 'PlotFcns',@optimplotx);
            options = optimset(options, 'OutputFcn',@BLFit_outfun);
            options = optimset(options, 'MaxIter',500);
            options = optimset(options, 'MaxFunEvals', 1000);
            fit.options = options;
        end
        prompt = {'Tolerance:','Maximum Iteration:', 'MaxFunEvals'};
        title = 'Change Your Fitting Options';
        dims = [1 55];
        definput = {sprintf('%0.2e', options.TolX),... 
            sprintf('%i', options.MaxIter),... 
            sprintf('%i', options.MaxFunEvals)};
        answer = inputdlg(prompt,title,dims,definput);
        if ~isempty(answer)
            options.TolX = str2double(answer{1});
            options.MaxIter = str2double(answer{2});
            options.MaxFunEvals = str2double(answer{3});
            fit.options = options;
            assignin('base', 'fit', fit);
        end
    end
    function pushbutton_reference_Callback(varargin)
        FitLee_helpstr = fitFunHandle('help');
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        if ~isempty(FitLee_helpstr)
            h = mymsgbox(FitLee_helpstr, 'Help', CreateStruct);
        else
            h = msgbox('No help available', 'Help');
        end
    end
    function edit_UB_callback(varargin)
        %fit = evalin('base', 'fit');
        fit = getappdata(gcbf, 'Fit');
        BLFit_checkinlimit(fit.handles);
    end
    function edit_LB_callback(varargin)
        %fit = evalin('base', 'fit');
        fit = getappdata(gcbf, 'Fit');
        BLFit_checkinlimit(fit.handles);
    end
    function pushbutton_simulate_Callback(varargin)

    % fit parameter reading
        %fit = evalin('base', 'fit');
        fit = evalin('base', 'fit');
        if ~isempty(fit)
            handles = fit.handles;
            FitLeehandle = fit.handles.FitLee;
        else
            if ~isempty(gcbf)
                FitLeehandle = gcbf;
            else
                FitLeehandle = handles.FitLee;
            end
            fit = getappdata(FitLeehandle, 'Fit');
        end
        
        handles = guihandles(FitLeehandle);
        [p, ~, ~, ~, ~, p_str] = BLFit_readparameters(handles);
        
        loadfitparam = 0 ;
        % fit parameter reading
        if sum(p) == 0
            loadfitparam = 1;
        end
        if sum(isnan(p))
            loadfitparam = 1;
        end
        if loadfitparam
            if isfield(fit, 'param')
                p = fit.param;
            else
                BLFit_setparameters(handles, p);
            end
        end

        try 
            pre_fitFunHandle = evallin('base', 'pre_userfit');
        catch
            pre_fitFunHandle = [];
        end
        if ~isempty(pre_fitFunHandle)
            pre_fitFunHandle(p_str);
        end
        
        if isfield(fit, 'tag')
            datatag = fit.tag;
        else
            datatag = [];
        end
        
        NdataSet = fit.NdataSet;
        % if NdataSet is not 0, multiple simultanous fit...
        if isfield(fit, 'simlineh')
            try
                delete(fit.simlineh);
            end
        end
        ax = findobj(fit.figh, 'type', 'axes');
        delete(findobj(ax, 'tag', 'temporary'));
        for i=1:NdataSet
            data2fit = fit.data{i};
            q{i} = data2fit(:,1);
            fit.simlineh(i) = line('xdata', data2fit(:,1), 'ydata', ...
                data2fit(:,2), 'color','m', 'tag', [fit.simtag, num2str(i)],...
                'parent', ax);
        end
        if numel(q) == 1
            q = q{1};
        end
        %end
        Iq = fitFunHandle(p_str, q, datatag);
        if iscell(Iq)
            for i=1:numel(Iq)
                fit.fit{i} = [fit.xd{i}, Iq{i}, fit.yd{i}];
            end
        else
            if numel(Iq) == length(size(Iq))
                fit.fit = [fit.xd(:), Iq(:), fit.yd(:)];
                setplot(fit.simlineh, Iq)
            else
                fit.fit = [fit.xd(:), Iq, fit.yd(:)];
                setplot(fit.simlineh, Iq(:,1))
                colr = jet(size(Iq, 2));
                for i=1:size(Iq, 2)-1
                    line('xdata', q, 'ydata', Iq(:,i+1), ...
                        'parent', ax, 'color', colr(i, :), ...
                        'linestyle', '-.', 'tag', 'temporary');
%                    line('xdata', q, 'ydata', Iq(:,i+1),'parent', ax, 'color', colr{i}, 'linestyle', '-.', 'tag', 'temporary');
                end
            end
        end
        fit.param = p;
        setappdata(FitLeehandle, 'Fit', fit);
        assignin('base', 'fit', fit);
    end
    function pushbutton_showAuxresult(varargin)
        fit = evalin('base', 'fit');
        FitLeehandle = gcbf;
        %fit = getappdata(FitLeehandle, 'Fit');
        handles = guihandles(FitLeehandle);

        [p, ~, ~, ~, ~, p_str] = BLFit_readparameters(handles);
        
        loadfitparam = 0 ;
        % fit parameter reading
        if sum(p) == 0
            loadfitparam = 1;
        end
        if sum(isnan(p))
            loadfitparam = 1;
        end
        if loadfitparam
            p = fit.param;
            BLFit_setparameters(handles, p);
        end

        try 
            pre_fitFunHandle = evallin('base', 'pre_userfit');
        catch
            pre_fitFunHandle = [];
        end
        if ~isempty(pre_fitFunHandle)
            pre_fitFunHandle(p_str);
        end
        
        if isfield(fit, 'tag')
            datatag = fit.tag;
        else
            datatag = [];
        end

        % define q
%        q = fit.xd(fit.roi(1));
        q = fit.xd(fit.roi(1):fit.roi(2));
        % Execute the function to run the auxilary function in each
        % fitting function if any.
        try
            q = q(:);
            [Iq, res] = fitFunHandle(p_str, q, datatag);
        catch
            fprintf('This fitting code does not have any auxilary function.\n');
        end
           
    end
    function pushbutton_Fit_Callback(varargin)
        % This function will fit fit.data to any model.
        % Setup data to fit.
        % load data and roi
        fit = evalin('base', 'fit');
        if ~isempty(fit)
%        if ~exist('handles')
%            fit = evalin('base', 'fit');
            handles = fit.handles;
            FitLeehandle = fit.handles.FitLee;
        else
            if ~isempty(gcbf)
                FitLeehandle = gcbf;
            else
                FitLeehandle = handles.FitLee;
            end
            
            fit = getappdata(FitLeehandle, 'Fit');
        end
        handles = guihandles(FitLeehandle);
        delete(findobj(fit.figh, '-regexp', 'tag', '__fit?'))
        delete(findobj(fit.figh, '-regexp', 'tag', '__sim?'))
        
        ax = findobj(fit.figh, 'type', 'axes');
        delete(findobj(ax, 'tag', 'temporary'));

        if isfield(fit, 'roi')
            roi = fit.roi;
        else
            roi = [];
        end
        
%        brush = logical(get(ax, 'BrushData'));
%         xd = get(ax1, 'XData');
%         yd = get(ax1, 'YData');
%         brushed_x = xd(brush);
%         brushed_y = yd(brush);
%        ind = linspace(1, numel(
        NdataSet = fit.NdataSet;
        % if NdataSet is not 0, multiple simultanous fit...
        for i=1:NdataSet
            data2fit = fit.data{i};
            brush = get(fit.datahandle(i), 'brushdata');
            if isempty(brush)
                brush = ones(1, numel(data2fit(:,1)));
            end
            ind = ones(size(brush));
            if ~isempty(roi)
                roi = sort(roi(i, 1:2));
                ind(1:(roi(1)-1)) = 0;
                ind((roi(2)+1):end) = 0;
            end
            ind = ind(:) & brush(:);
            q{i} = data2fit(ind(:),1);
            y{i} = data2fit(ind(:),2);

            fit.fitlineh(i) = line('xdata', q{i}, 'ydata', ...
                y{i}, 'color','r', 'tag', [fit.fittag, num2str(i)],...
                'parent', findobj(fit.figh, 'type', 'axes'));
        end
        if isfield(fit, 'tag')
            datatag = fit.tag;
        else
            datatag = [];
        end
        assignin('base', 'fit', fit)

%         if nargin < 3
%             handles = fit.handles;
%         end
        
        % If there is only a set of data of interest.
        if numel(q) == 1
            q = q{1}(:);
            y = y{1}(:);
        end
        
        % Conditioning data.
        k = isnan(y(:).*q(:));
        y(k) = [];
        q(k) = [];
        k = isinf(y(:).*q(:));
        y(k) = [];
        q(k) = [];
        

        %% Setup fitting parameters.
        % fit parameter reading ==============================
        [p, qfit, LB, UB, pnames, p_str] = BLFit_readparameters(handles);
        pnames = pnames(:);% for fitting p should be column array, so be pnames.

        %% Run pre-fit function

        try 
            pre_fitFunHandle = evallin('base', 'pre_userfit');
        catch
            pre_fitFunHandle = [];
        end
        
        if ~isempty(pre_fitFunHandle)
            pre_fitFunHandle(p_str);
        end
        
        %% Setup Constraints for fitting.
        % Load matrices for contraints from 'base'
        NONLCON = [];
        isLinIneq = 1;
        try
            A = evalin('base', 'A');
            B = evalin('base', 'B');
            if size(A, 2) ~= numel(p)
                isLinIneq = 0;
            end
            if size(A, 1) ~= numel(B)
                isLinIneq = 0;
            end
        catch
            isLinIneq = 0;
        end
        if ~isLinIneq
            A = [];B = [];
        end
        isLinEq = 1;
        try
            Aeq = evalin('base', 'Aeq');
            Beq = evalin('base', 'Beq');
            if size(Aeq, 2) ~= numel(p)
                isLinEq = 0;
            end
            if size(Aeq, 1) ~= numel(Beq)
                isLinEq = 0;
            end
        catch
            isLinEq = 0;
        end
        if ~isLinEq
            Aeq = [];Beq = [];
        end
        if isLinEq % if there is linear equality contraint, release other constraint.
            A = [];B = [];
        end
        
        % Vary or not ===========
         % ====================================================
         
        for i=1:numel(p)
            if ~qfit(i) % not checked for fitting then LB = UB = p
                LB(i) = p(i);
                UB(i) = p(i);
            else 
                vchanged = 0;
                if isnan(LB(i))
                    LB(i) = p(i) - abs(p(i))*0.5;
                    vchanged = 1;
                end
                if isnan(UB(i))
                    UB(i) = p(i) + abs(p(i))*0.5;
                    vchanged = 1;
                end

                if vchanged == 1
                    set(handles.(sprintf('edit_LB%i', i)), 'string', LB(i));
                    set(handles.(sprintf('edit_UB%i', i)), 'string', UB(i));
                end
            end
        end

        % ==========================        
        % Constraint done .................................................
        %% Setup options
        if isfield(fit, 'options')
            options = fit.options;
        else
            options = optimset('fminsearch');
            options = optimset(options, 'TolX',0.1E-6);
    %        options = optimset(options, 'PlotFcns',@optimplotx);
            options = optimset(options, 'OutputFcn',@BLFit_outfun);
            options = optimset(options, 'MaxIter',500);
            options = optimset(options, 'MaxFunEvals', 1000);
            fit.options = options;
        end
        
        %% Data fitting
        qfit = logical(qfit);
        p2fit = p(qfit);
        LB2fit = LB(qfit);
        UB2fit = UB(qfit);
        
        NLPstart = p2fit';
        LB = LB2fit;
        UB = UB2fit;
        
        pN{:,1} = pnames;
        pN{:,2} = p';
        pN{:,3} = qfit';
%         if exist('fmincon', 'file')
%             %options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');
%            INLP = fmincon(@(x) fitfunction(x, y, q, pnames, datatag),NLPstart,A,B,Aeq,Beq,LB,UB,NONLCON,options);
%         else
%            INLP = fminsearchcon(@(x) fitfunction(x, y, q, pnames, datatag),NLPstart,LB,UB, A, B, [], options);
%         end
        INLP = fminsearchcon(@(x) fitfunction(x, y, q, pN, datatag),NLPstart,LB,UB, A, B, [], options);
        INLP = restoreParameters(INLP, pN);
        fit.param = INLP;
        BLFit_setparameters(handles, INLP)
        setappdata(FitLeehandle, 'Fit', fit)
        assignin('base', 'fit', fit);
        
        %% Aftermath
        pushbutton_simulate_Callback

        % Return results ===========
    end
    function [p, pnames] = restoreParameters(p, pN)
        pv = pN{:,2};
        qfit = pN{:,3};
        pnames = pN{:,1};
        pv(qfit) = p;
        p = pv;
    end
        
    function cv = fitfunction(p, y, q, pnames, tag)

        % Calculation function
        %[predY, fit] = multivoight(parameters, q, datatag);
        % q : it can be either array or a cell for multiple data set
        % y : it can be either array or a cell for multiple data set
        % q and y should have the same length.
        % datatag: strings that are tags of plotted curves. 
        %       this can be either array or a cell for multiple data set.
        [p, pnames] = restoreParameters(p, pnames);
        p_str = cell2struct(num2cell(p), pnames, 1);
        Iq = fitFunHandle(p_str, q, tag);
        assignin('base', 'Iq', Iq);
        if iscell(Iq)
            Iqarray = cell2mat(Iq);
            Iqarray = Iqarray(:);
            yarray = cell2mat(y);
            yarray = yarray(:);
        else
            if (size(Iq, 2) > 1) && (size(Iq, 1)>size(Iq, 2))
                Iqarray = Iq(:,1);
            else 
                Iqarray = Iq(:);
            end
            yarray = y;
        end
        cv = chi_squared(yarray, Iqarray, 5);
        %fit.fit = [q(:), predY(:), y(:)];
    end

    function pushbutton_Save_Callback(varargin)
        fn = get(handles.edit_savefilename, 'string');
        BLFit_save(handles, fn);
    end
    function pb_loaddata_Callback(varargin)
        BLFit_Loaddata(handles);
        BLFit_pickFigure(str2double(get(handles.edit_Figure, 'string')), handles)
    end
    function pb_plotstyle_Callback(varargin)
        BLFit_plotstyle(handles);
    end
    function pb_LL_Callback(varargin)
        hObject = varargin{1};
        BLFit_Indexcallback(hObject, handles);
    end
    function pb_LU_Callback(varargin)
        hObject = varargin{1};
        BLFit_Indexcallback(hObject, handles);
    end
    function pb_RL_Callback(varargin)
        hObject = varargin{1};
        BLFit_Indexcallback(hObject, handles);
    end
    function pb_RU_Callback(varargin)
        hObject = varargin{1};
        BLFit_Indexcallback(hObject, handles);
    end
    function pushbutton_qselect_Callback(varargin)
        BLFit_qselect(gcbf);
        %setappdata(gcbf, 'Fit', fit0);
    end
    function edit_Lindex_Callback(varargin)
        BLFit_index_callback(gcbf);
    end
    function edit_Rindex_Callback(varargin)
        BLFit_index_callback(gcbf);
        %setappdata(gcbf, 'Fit', fit0);
    end 
    function edit_Lq_Callback(varargin)
        BLFit_qindex_callback(gcbf);
        %setappdata(gcbf, 'Fit', fit0);
    end
    function edit_Rq_Callback(varargin)
        BLFit_qindex_callback(gcbf);
        %setappdata(gcbf, 'Fit', fit0);
    end

    function edit_Figure_CreateFcn(varargin)
        %hObject = handles.edit_Figure;
        %hObject = varargin{1};
        if ~isempty(call_obj)
            fprintf('It is called from %s\n', get(call_obj, 'Tag'))
            switch lower(get(call_obj, 'Tag'))
                case 'saxslee_fig'
                    set(varargin{1}, 'string', 'SAXSLee');
                    set(varargin{1}, 'style', 'pushbutton');
                case 'load_fitlee'
                    
            end
%            set(varargin{1}, 'Callback', 'edit_Figure_update');
        end
        %Figh = str2double(get(hObject, 'string'));
        %BLFit_pickFigure(Figh, handles)
    end

    function edit_Figure_Callback(varargin)
        if ~exist('handles')
            fit = evalin('base', 'fit');
            %fit = getappdata(FitLeehandle, 'Fit');
            handles = fit.handles;
        end
        hObject = handles.edit_Figure;
        %hObject = varargin{1};
        if strcmp(get(hObject, 'style'), 'pushbutton')
            Figh = call_obj;
        else
            Figh = str2double(get(hObject, 'string'));
        end
        BLFit_pickFigure(Figh, handles);
        Fit = evalin('base', 'fit');
        try
            setappdata(gcbf, 'Fit', Fit);
        end
    end

    function pushbutton_Load_Callback(varargin)
        BLFit_LoadSaved(handles);
        fit = evalin('base', 'fit');
        set(handles.edit_Figure, 'string', get(fit.figh, 'Number'));
        BLFit_pickFigure(str2double(get(handles.edit_Figure, 'string')), handles)
    end

end

