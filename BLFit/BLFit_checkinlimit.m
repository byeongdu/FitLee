function BLFit_checkinlimit(handles)
tol = 0.001; % tolerance
%if ishandle(handles)
hnames = fieldnames(handles);
nump = numel(cell2mat(strfind(hnames, 'checkbox_Fit')));
for i=1:nump
    p = str2double(get(handles.(sprintf('edit_p%i', i)), 'string'));
%        set(findobj(handles, 'tag', sprintf('checkbox_Fit%i', i)), 'value', qfit(i));
    LB = str2double(get(handles.(sprintf('edit_LB%i', i)), 'string'));
    UB = str2double(get(handles.(sprintf('edit_UB%i', i)), 'string'));
    if p <=(LB+tol*abs(LB))
        set(handles.(sprintf('edit_LB%i', i)), 'Foreground', 'r');
    else
        set(handles.(sprintf('edit_LB%i', i)), 'Foreground', 'k');
    end
    if p>=(UB-tol*abs(UB))
        set(handles.(sprintf('edit_UB%i', i)), 'Foreground', 'r');
    else
        set(handles.(sprintf('edit_UB%i', i)), 'Foreground', 'k');
    end
end
%end