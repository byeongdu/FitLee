function FitLee_simultaneousfit(FitLeehandles, EqParam, ParamRatio)
% for simultaneous fitting of multiple data with the same model
% FitLee_simultaneousfit(FitLeehandles, Indices of Paramters to be the same, ParamRatio)
% Examples
% obj1 = FitLee('FitLee_coreshell_schultzsphere')
% obj2 = FitLee('FitLee_coreshell_schultzsphere')
% Load data and arrange q ranges to fit. Also put parameters and adjust
% parameters manually until satisfactory. Now assuming you want to keep 
% the parameter 1, 3, 4 across all objects the same.
% FitLee_simultaneousfit([obj1, obj2], [1,3,4]);
%
% When you would like to keep the ratio of parameter1 of obj1 over
% parameter1 of obj2 to 1.2,
% FitLee_simultaneousfit([obj1, obj2], [1,3,4], [1, 1.2]);
% When you would also like to keep the ratio of parameter5 of obj1 over
% parameter5 of obj2 to 1.3,
% FitLee_simultaneousfit([obj1, obj2], [1,3,4], [1, 1.2;5 1.3]);
% 
% When there are more than 2 objects, the ParamRatio matrix should be
% [Parameter Index, the Parameter of obj1 / the parameter of obj2, 
%    the Parameter of obj1 / the parameter of obj3, ...];
%
% ParamRatio = [N_parameter to have fixed ratio, N of FitLeehandles]
% for example, [1, ratio1, ratio2;
%               3, ratio1, ratio2];
% This is to use equality constraint.
% These parameters should be qfit and not be ParamEqual.

Nhandle = numel(FitLeehandles);
%fit = evalin('base', 'fit');
fitall.NdataSet = Nhandle;

NLPstart = [];
datatagf = [];
LBf = [];
UBf = [];
%ParamEqual = logical(ParamEqual);

for ind=1:Nhandle
    fit = getappdata(FitLeehandles(ind), 'Fit');
    delete(findobj(fit.figh, '-regexp', 'tag', '__fit?'));
    delete(findobj(fit.figh, '-regexp', 'tag', '__sim?'));
        
    ax = findobj(fit.figh, 'type', 'axes');
    delete(findobj(ax, 'tag', 'temporary'));

        if isfield(fit, 'roi')
            roi = fit.roi;
        else
            roi = [];
        end
        
        NdataSet = fit.NdataSet;
        if NdataSet > 1
            error('Works only when the number of data set in a Figure is 1.\n');
        end
        % if NdataSet is not 0, multiple simultanous fit...
        data2fit = fit.data{1};
        if isempty(roi)
            q{ind} = data2fit(:,1);
            y{ind} = data2fit(:,2);
        elseif numel(roi) == 2
            roifit = sort(roi(1, 1:2));
            q{ind} = data2fit(roifit(1):roifit(2),1);
            y{ind} = data2fit(roifit(1):roifit(2),2);
        else
            roifit = sort(roi(1, 1:2));
            q{ind} = data2fit(roifit(1):roifit(2),1);
            y{ind} = data2fit(roifit(1):roifit(2),2);
        end
        fit.fitlineh = line('xdata', q{ind}, 'ydata', ...
            y{ind}, 'color','r', 'tag', [fit.fittag, num2str(1)],...
            'parent', findobj(fit.figh, 'type', 'axes'));
        fitall.fitlineh(ind) = fit.fitlineh;
        
        if isfield(fit, 'tag')
            datatag = fit.tag;
        else
            datatag = [];
        end
        assignin('base', 'fit', fitall)

        handles = guihandles(FitLeehandles(ind));
        % If there is only a set of data of interest.
        
       

        %% Setup fitting parameters.
        % fit parameter reading ==============================
        [p, qfit, LB, UB, pnames, p_str] = BLFit_readparameters(handles);
        if ind == 1
            ParamEqual = zeros(1, numel(p));
            ParamEqual(EqParam) = 1;
            ParamEqual = logical(ParamEqual);
        end
        pnames = pnames(:);% for fitting p should be column array, so be pnames.

        %% Run pre-fit function

        try 
            pre_fitFunHandle = evalin('base', 'pre_userfit');
        catch
            pre_fitFunHandle = [];
        end
        
        if ~isempty(pre_fitFunHandle)
            pre_fitFunHandle(p_str);
        end
        
        %% Setup Constraints for fitting.
        % Load matrices for contraints from 'base'
        NONLCON = [];
        
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
        
        if ind==1
            pN{1} = pnames;
            pN{2} = p';
            pN{3} = qfit';
            pN{4} = ParamEqual(:);
            pN{5} = Nhandle;
        else
            % if parameters are not the same make equal.
            qfit = qfit & ~ParamEqual;
%            p(ParamEqual) = pN{2}(ParamEqual);
%            qfit(ParamEqual) = pN{3}(ParamEqual);
        end
        qfit = logical(qfit);
        p2fit = p(qfit)';
        LB2fit = LB(qfit)';
        UB2fit = UB(qfit)';
        
        
        NLPstart = [NLPstart; p2fit];
        datatagf{ind} = datatag;
        LBf = [LBf; LB2fit];
        UBf = [UBf; UB2fit];
        fitFunHandle{ind} =  str2func(handles.FitLee.Name);

end

% % Producing Equality constraint for keeping ratio.
% ParamEqual = ParamEqual(qfit);
% toFix = ParamEqual > 0;
% Nfitparam = sum(qfit);
Aeq = [];
Beq = [];
if nargin>2
    NtoFix = size(ParamRatio, 1);
    Nfitparam = numel(NLPstart);
    Aeq = zeros((Nhandle-1)*NtoFix, Nfitparam);
    Beq = zeros((Nhandle-1)*NtoFix, 1);
    %nInd = findindexParam(ind, pN);
    %Nratio = size(ParamRatio, 1);
    for ind = 1:NtoFix
        %nInd = zeros(1, Nhandle);
        for i=1:Nhandle
            nInd = findindexParam([i, ParamRatio(ind, 1)], pN);
            if i==1
                Aeq((Nhandle-1)*(ind-1)+[1:(Nhandle-1)], nInd) = 1;
            else
                Aeq((Nhandle-1)*(ind-1)+(i-1), nInd) = -1*ParamRatio(ind, i);
            end
        end
    end
end

        % ==========================        
        % Constraint done .................................................
        %% Setup options
        options = optimset('algorithm', 'levenberg-marquardt');
        options = optimset(options, 'TolX',0.1E-6);
        options = optimset(options, 'OutputFcn',@BLFit_outfun);
        options = optimset(options, 'MaxIter',500);
        options1 = optimset(options, 'MaxFunEvals', 1000);

        options2 = setoptimoptions('algorithm', 'levenberg-marquardt',...       
    'hessupdate', 'bfgs',...
    'Largescale', 'on',...
    'storeN', 20,...
        'TolX' , 1e-7,...
        'OutputFcn',@BLFit_outfun,...
        'TolFun' ,1e-13,...
        'Gradobj', 'on',...
        'GradConstr', 'on',...
        'MaxIter', 1000,...
        'MaxFunEvals', 5e4,...
        'popsize', 10);
        
        %% Data fitting
        rosen = @(x) fitfunction(x, y, q, pN, datatagf, fitFunHandle);
%        INLP = fminsearchcon(@(x) fitfunction(x, y, q, pnamesf, datatagf, fitFunHandle),NLPstart,LBf,UBf, Aeq, Beq, [], options1);
    if isempty(Aeq)
        INLP = fminsearchcon(rosen,NLPstart,LBf,UBf, [], [], [], options1);
        %var = {y, q, pnamesf, datatagf, fitFunHandle};
    else
        INLP = minimize(rosen,NLPstart,[],[], Aeq, Beq, LBf, UBf, [], options2);
    end

fit.param = INLP;

Nhandle = numel(q);
pv = restoreParameters(INLP, pN);
for ind=1:Nhandle
    %rng = ((ind-1)*Nparam+1:ind*Nparam);
    %pv = INLP(rng);
    handles = guihandles(FitLeehandles(ind));
    BLFit_setparameters(handles, pv{ind})
end
        %% Aftermath
%        pushbutton_simulate_Callback
end
function [p, pnames] = restoreParameters(p, pN)
    pvoriginal = pN{2};
    qfit = pN{3};
    pnames = pN{1};
    pv = {};
    toFix = pN{4};
    Nhandle = pN{5};
    Fitted = qfit & ~toFix;
    for i=1:Nhandle
        if i==1
            pv{i} = pvoriginal;
            idx = qfit;
        else
            pv{i} = pv{1};
            idx = Fitted;
        end
        idx = logical(idx);
        rng = 1:sum(idx);
        pv{i}(idx) = p(rng);
        p(rng) = [];
    end
    p = pv;
end

function nInd = findindexParam(ind, pN)
% ind is a coordinate of [FitLee object, Parameter Number]
% nInd is the index of the parameter ind in the new parameter for fitting.
    %pvoriginal = pN{2};
    qfit = pN{3};
    %pnames = pN{1};
    toFix = pN{4};
    Nhandle = pN{5};
    Fitted = qfit & ~toFix;
    if ind(1) > 1
        nInd = sum(qfit);
        nInd = nInd + sum(Fitted)*(ind(1)-2) + ind(2);
    else
        nInd = ind(2);
    end
end

function cv = fitfunction(p, varargin)

    y = varargin{1};
    q = varargin{2};
    pnames = varargin{3};
    tag = varargin{4};
    fitFunHandle = varargin{5};
    [p, pnames] = restoreParameters(p, pnames);
    % Calculation function
    %[predY, fit] = multivoight(parameters, q, datatag);
    % q : it can be either array or a cell for multiple data set
    % y : it can be either array or a cell for multiple data set
    % q and y should have the same length.
    % datatag: strings that are tags of plotted curves. 
    %       this can be either array or a cell for multiple data set.
    yv = [];
    Iv = [];
    Iqout = [];
    NDataSet = numel(q);
%    Nparam = numel(p)/NDataSet;
    for i=1:NDataSet
%        rng = ((i-1)*Nparam+1:i*Nparam);
%        p_str = cell2struct(num2cell(p(rng)), pnames(rng), 1);
        p_str = cell2struct(num2cell(p{i}), pnames, 1);
        fitfunc = fitFunHandle{i};
        Iq = fitfunc(p_str, q{i}, []);
        Iqout{i} = Iq;
        if iscell(Iq)
            Iqarray = cell2mat(Iq);
            Iqarray = Iqarray(:);
            yarray = cell2mat(y{i});
            yarray = yarray(:);
        else
            if (size(Iq, 2) > 1) && (size(Iq, 1)>size(Iq, 2))
                Iqarray = Iq(:,1);
            else 
                Iqarray = Iq(:);
            end
            yarray = y{i};
        end
        yv = [yv(:); yarray(:)];
        Iv = [Iv(:); Iqarray(:)];
    end
    cv = chi_squared(yv, Iv, 5);
    assignin('base', 'Iq', Iqout);

    %fit.fit = [q(:), predY(:), y(:)];
end
