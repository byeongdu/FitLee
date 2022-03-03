function startFitLee(varargin)
warning off
verNumber = '1.0';
fh = init;
%resize(fh);

%% Menu
hMenuRun = uimenu(fh,...
    'Label','Run',...
    'Position',1,...
    'HandleVisibility','callback',...
    'Tag','SSL_MenuDet',...
    'callback', @run);

hMenuHelp = uimenu(fh,...
    'Label','Help',...
    'Position',2,...
    'HandleVisibility','callback',...
    'Tag','SSL_MenuHelp',...
    'callback', @myhelp);

    function run(varargin)
        hn = findobj(gcbf, 'tag', 'popFuncNames');
        ho = findobj(gcbf, 'tag', 'edOptions');
        funslist = get(hn, 'string');
        [~, Fname] = fileparts(funslist{get(hn, 'value')});
        options = get(ho, 'string');
        if ~isempty(options)
            op = eval(options);
        else
            op = [];
        end
        
        if ~isempty(op)
            FitLee(Fname, op);
        else
            FitLee(Fname);
        end
    end
    function myhelp(varargin)
        hn = findobj(gcbf, 'tag', 'popFuncNames');
        funslist = get(hn, 'string');
        Fname = funslist{get(hn, 'value')};
        eval(sprintf('doc %s', Fname)); 
    end

function f = init(varargin)
pth = fileparts(which('startFitLee'));
FF = dir([pth, filesep, 'FitLee_*']);
Funcs = {};
for i=1:numel(FF)
    [~, Funcs{i}] = fileparts(FF(i).name);
end

f = figure;
set(f, 'NumberTitle', 'off')
set(f, 'Tag', 'Load_FitLee')
set(f, 'Name', 'Load FitLees') 
fpos = get(f, 'position');
fpos(3) = 500;
fpos(4) = 60;
set(f, 'position', fpos);
set(f, 'MenuBar', 'none');
%tbh = uitoolbar(f);
%set(f, 'uitoggletool', 'on');
%set(f, 'resizeFcn', @resize)

hMenubar = findall(gcf,'tag','figMenuFile');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(gcf,'tag','figMenuView');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(gcf,'tag','figMenuEdit');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(gcf,'tag','figMenuHelp');
set(findall(hMenubar), 'visible', 'off')

%get(findall(hMenubar),'tag')
hToolbar = findall(gcf,'tag','FigureToolBar');
hStandardTools = findall(hToolbar,'-regexp', 'tag', 'Standard');
set(hStandardTools, 'visible', 'off');
%get(findall(hToolbar),'tag')

f_txtnames = uicontrol('Parent', f, ...
    'tag' , 'txtFuncNames', ...
    'style', 'text', ...
    'position', [10, 30, 50, 20], 'string', 'Names: ');
f_ednames = uicontrol('Parent', f, ...
    'tag' , 'popFuncNames', ...
    'style', 'popupmenu', ...
    'position', [60, 25, 220, 25], 'string', Funcs);
f_txtoptions = uicontrol('Parent', f, ...
    'tag' , 'txtOptions', ...
    'style', 'text', ...
    'position', [10, 10, 50, 20], 'string', 'Options: ');
f_edOptions = uicontrol('Parent', f, ...
    'tag' , 'edOptions', ...
    'style', 'edit', ...
    'position', [60, 10, 100, 20], 'string', '');
end
end
