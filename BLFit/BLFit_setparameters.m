function BLFit_setparameters(handles, p, qfit, LB, UB)
%fit = evalin('base', 'fit');
ud = get(gcbf, 'userdata');
if ~isempty(ud)
    fnames = fieldnames(ud);
    pN = numel(p);
    nF = numel(fnames);
    k = 1;
    for i=1:numel(fnames)
        if contains_replace(fnames{i}, 'no_display')
            ud.(fnames{i}) = p(pN-nF+i);
        end
    end
    set(gcbf, 'userdata', ud);
end

if ishandle(handles)
    nump = findobj(handles, 'style', 'checkbox');
    for i=1:numel(nump)
        set(findobj(handles, 'tag', sprintf('edit_p%i', i)), 'string', p(i));
%        set(findobj(handles, 'tag', sprintf('checkbox_Fit%i', i)), 'value', qfit(i));
%        set(findobj(handles, 'tag', sprintf('edit_LB%i', i)), 'string', LB(i));
%        set(findobj(handles, 'tag', sprintf('edit_UB%i', i)), 'string', UB(i));
    end
    return
end

hnames = fieldnames(handles);
nump = numel(cell2mat(strfind(hnames, 'checkbox_Fit')));
%nump = numel(findobj(gcbf, 'Tag', 'checkbox_Fit'));
for i=1:nump
    if abs(p(i)) < 1E-2
        pv = sprintf('%0.4e', p(i));
    else
        pv = sprintf('%0.6f', p(i));
    end
    
    if nargin == 2
        set(handles.(sprintf('edit_p%i', i)), 'string', pv);
    elseif nargin ==5
        set(handles.(sprintf('edit_p%i', i)), 'string', pv);
        set(handles.(sprintf('checkbox_Fit%i', i)), 'value', qfit(i));
        set(handles.(sprintf('edit_LB%i', i)), 'string', LB(i));
        set(handles.(sprintf('edit_UB%i', i)), 'string', UB(i));
    end
end

BLFit_checkinlimit(handles)