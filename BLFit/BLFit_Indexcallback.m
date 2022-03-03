function BLFit_Indexcallback(h, handles)
Lind = str2double(get(handles.edit_Lindex, 'string'));
Rind = str2double(get(handles.edit_Rindex, 'string'));
FitLeehandle = handles.FitLee;
fit = getappdata(FitLeehandle, 'Fit');
%fit = evalin('base', 'fit');
Numpnt = numel(fit.data{1}(:,1));
if isempty(Lind)
    Lind = 0;
end
if isempty(Rind)
    Rind = Numpnt;
end

switch get(h, 'tag')
    case 'pb_LL'
        Lind = Lind - 1;
        if Lind < 1
            Lind = 1;
        end
        set(handles.edit_Lindex, 'string', num2str(Lind));
    case 'pb_LU'
        Lind = Lind + 1;
        set(handles.edit_Lindex, 'string', num2str(Lind));
    case 'pb_RL'
        Rind = Rind - 1;
        set(handles.edit_Rindex, 'string', num2str(Rind));
    case 'pb_RU'
        Rind = Rind + 1;
        if Rind > Numpnt
            Rind = Numpnt;
        end
        set(handles.edit_Rindex, 'string', num2str(Rind));
end
BLFit_index_callback(FitLeehandle)
