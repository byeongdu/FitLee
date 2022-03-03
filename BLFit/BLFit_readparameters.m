function [p, qfit, LB, UB, pnames, p_str] = BLFit_readparameters(handles)

if ishandle(handles)
    nump = findobj(handles, 'style', 'checkbox');
    for i=1:numel(nump)
        p(i) = str2double(get(findobj(handles, 'tag', sprintf('edit_p%i', i)), 'string'));
        qfit(i) = get(findobj(handles, 'tag', sprintf('checkbox_Fit%i', i)), 'value');
        LB(i) = str2double(get(findobj(handles, 'tag', sprintf('edit_LB%i', i)), 'string'));
        UB(i) = str2double(get(findobj(handles, 'tag', sprintf('edit_UB%i', i)), 'string'));
        pnames{i} = get(findobj(handles, 'tag', sprintf('text_p%i', i)), 'string');
    end
    
    ud = get(gcbf, 'userdata');
    if ~isempty(ud)
        k = 1;
        fnames = fieldnames(ud);
        for i=1:numel(fnames)
            if contains_replace(fnames{i}, 'no_display')
                p((nump)+k) = ud.(fnames{i});
                qfit((nump)+k) = 1;
                LB((nump)+k) = p((nump)+k)/10000;
                UB((nump)+k) = p((nump)+k)*10000;
                pnames{(nump)+k} = fnames{i};
                k = k + 1;
            end
        end
    end
    
    p_str = cell2struct(num2cell(p), pnames, 2);
    return
end

hnames = fieldnames(handles);
nump = numel(cell2mat(strfind(hnames, 'checkbox_Fit')));
for i=1:nump
    p(i) = str2double(get(handles.(sprintf('edit_p%i', i)), 'string'));
    qfit(i) = get(handles.(sprintf('checkbox_Fit%i', i)), 'value');
    LB(i) = str2double(get(handles.(sprintf('edit_LB%i', i)), 'string'));
    UB(i) = str2double(get(handles.(sprintf('edit_UB%i', i)), 'string'));
    pnames{i} = get(handles.(sprintf('text_p%i', i)), 'string');
end

ud = get(gcbf, 'userdata');
if ~isempty(ud)
    k = 1;
    fnames = fieldnames(ud);
    for i=1:numel(fnames)
        if contains_replace(fnames{i}, 'no_display')
            p((nump)+k) = ud.(fnames{i});
            qfit((nump)+k) = 1;
            LB((nump)+k) = p((nump)+k)/10000;
            UB((nump)+k) = p((nump)+k)*10000;
            pnames{(nump)+k} = fnames{i};
            k = k + 1;
        end
    end
end

p_str = cell2struct(num2cell(p), pnames, 2);
