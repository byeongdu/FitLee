function [y, Model, ox] = maxent_pr(qvec, MeasuredData, Errors, x_radius, MaxIterations,MaxEntStabilityParam,FormMtrx, Background)
% function [y, Model, ox] = maxent(qvec, MeasuredData, Errors, x_radius,
% MaxIterations,MaxEntStabilityParam,FormMtrx)
% Example
% a = load('C:\Users\blee\Google Drive\Masaki\SE2PC_NaP50.bsub_avg');
% a = a(1:300, :);
% [y, md,ox] = maxent_pr(a(:,1), a(:,2)/0.01, sqrt(a(:,2)), 1:1:600, 10000, 0.00001,[],[]);
global Chisquare;
global chtarg;
global chizer;
global fSum;
global blank;
global betaMX;
global flWv;
global blWv;
global c1;
global c2;
global s1;
global s2;
if nargin<8
    Background = zeros(size(MeasuredData));
end
if isempty(Errors)
    Errors = 0.005*sqrt(MeasuredData);
end
if isempty(FormMtrx)
    %FormMtrx = sinc(qvec*x_radius);
    for kk=1:length(x_radius)
        FormMtrx(:,kk) = sinc(qvec*x_radius(kk)); 
        % this means determined distribution will be volume distribution..
    end
    FormMtrx = FormMtrx*4*pi; %% from Jan's software
end

InitialModelBckg = 1E-6;
Model = InitialModelBckg*ones(numel(x_radius), 1);
if (numel(InitialModelBckg)==1)
    InitialModelBckg = InitialModelBckg*ones(numel(Model), 1);
end
blank = InitialModelBckg(1);
[m, n] = size(MeasuredData);
if m<n
    MeasuredData=MeasuredData';
    Errors = Errors';
end
% maximum entropy method
% measure data : wave
% error : wave
% initialmodelBckg : distribution model (R), wave
% Model : distribution model (R)
% variable : maxiterations, maxEntstabilitypram

%   FuncRef IR1R_ModelOpus OpusFnct		//converts the Model to MeasuredData
%	FuncRef IR1R_ModelOpus TropusFnct	//converts the MeasuredData into Model
%	FuncRef IR1R_UpdateDataForGrph UpdateGraph	//converts the MeasuredData into Model
% OpusFnct(measured data, model)

% constant =============
    e = 2.71828;
%==============
    npnts=numel(MeasuredData);			%number of measure points
    nbins=numel(Model);					%number of bins in model
	
%	WAVE flWv=root:Packages:Sizes:flWv
%	WAVE blWv=root:Packages:Sizes:blWv
	flWv=zeros(3);
	blWv=zeros(1,3);

%	chtarg, chizer, fSum, blank, CurrentEntropy, CurrentChiSq, CurChiSqMinusAlphaEntropy, Chisquare
%	fSum=0;
	tolerance=MaxEntStabilityParam*sqrt(2*npnts); %for convergnence in Chisquare
	tstLim=0.05;		%for convergence for entropy terms

%	Chisquare=0;
	chizer = npnts;		%setup some starting conditions
	chtarg = chizer; 
    ctarg = 0;
%	iter=0; snorm=0; cnorm=0; tnorm=0; a=0; b=0; test=0; i=0; j=0; l=0; fchange=0; df=0; 
    sEntropy=0; 
%    k=0;

%    ox = MeasuredData; 
%    ascratch = ox;bscratch=ox; etaScratch=ox;zscratch=ox;zscratch2=ox; % create work waves with measured Points length
%    cgrad = Model; sgrad = Model; ModelScratch = Model; ModelScratch2 = Model; xiScratch=Model; %create work waves with bins length
    xi = zeros(numel(Model), 3);
    eta = zeros(numel(MeasuredData), 3);
    c1 = zeros(1,3);s1 = c1; 
    betaMX = zeros(1,3);
    c2 = zeros(3); s2 = zeros(3);
	
	for iter=0:(MaxIterations-1)		% this is the main loop which does the searching for solution
		ox = opus(Model, FormMtrx);     % calculate theoritical curve(ox) from Model
		Chisquare=0;
		ascratch = (ox - MeasuredData) ./ Errors;
		ox = 2 * ascratch ./ Errors;
		ascratch=ascratch.^2;
		Chisquare = sum(ascratch);
		cgrad = tropus(ox,FormMtrx);
        test=0;
		fSum = sum(Model);
		sgrad = -log(Model./InitialModelBckg) ./ (InitialModelBckg * e);
		ModelScratch = Model .* sgrad.^2;
		snorm = sum(ModelScratch);
        ModelScratch = Model .* cgrad.^2;
		cnorm = sum(ModelScratch);
		ModelScratch = Model .* sgrad .* cgrad;
		tnorm = sum(ModelScratch);
		
		snorm = sqrt(snorm);
		cnorm = sqrt(cnorm);
		a = 1;
		b = 1/cnorm;
		if (iter>0)
			test = sqrt(0.5*(1-tnorm/(snorm*cnorm)));
			a = 0.5 / (snorm * test);
			b = 0.5 * b / test;
        end
		xi(:,1) = Model .* cgrad / cnorm;
		xi(:,2) = Model .* (a * sgrad - b * cgrad);

		xiscratch=xi(:,1);
		eta(:,1) = opus(xiscratch, FormMtrx);
	
		xiscratch=xi(:,2);
		eta(:,2)=opus(xiscratch, FormMtrx);
		
		ox = eta(:,2) ./ Errors.^2;
		%ox = eta(:,2) ./ row2vect(Errors).^2+Background'; % BLee added march 21. 2011
		
		xi(:,3) = tropus(ox,FormMtrx);
        
		ModelScratch = Model .* xi(:,3);
		ModelScratch2 = ModelScratch .* xi(:,3);
		a = sum(ModelScratch2);
		xi(:,3) = ModelScratch;
		
		a= 1/sqrt(a);
		xi(:,3) = a .* xi(:,3);
		xiscratch=xi(:,3);
		eta(:,3) = opus(xiscratch,FormMtrx);
		for i=1:3
			s1(i) = sum(xi(:,i) .* sgrad);
			c1(i) = sum(xi(:,i) .* cgrad);
        end

        c1=c1/Chisquare;
		
		s2=zeros(3);
		c2=zeros(3);
		for k=0:2
			for l=0:k
%                s2(k+1,l+1) = -1*sum(xi(:,k+1).*xi(:,l+1)./row2vect(Model));
                for i=0:(nbins-1)
					s2(k+1,l+1) = s2(k+1,l+1) - xi(i+1,k+1) * xi(i+1,l+1) / Model(i+1);
                end
%                c2(k+1,l+1) = sum(eta(:,k+1).*eta(:,l+1)./row2vect(Errors).^2);
                for j=0:(npnts-1)
					c2(k+1,l+1) = c2(k+1,l+1) + eta(j+1,k+1) * eta(j+1,l+1) / ((Errors(j+1))^2);
                end
            end
        end	
        
		s2 = s2 / InitialModelBckg(1);
		c2 = 2 * c2 /Chisquare;
        c2(1,2)=c2(2,1);
        c2(1,3)=c2(3,1);
        c2(2,3)=c2(3,2);
        s2(1,2)=s2(2,1);
        s2(1,3)=s2(3,1);
        s2(2,3)=s2(3,2);
%		c2 = c2';
%        s2 = s2';
        betaMX = [-0.5 * c1(1) / c2(1), 0, 0];
%        c1, c2, s1, s2
        if(iter>0)
			IR1R_Move(3); 
        end
%        Modify the current distribution (f-vector)
%        fSum = 0;              % find the sum of the f-vector
%        fChange = 0;          % and how much did it change?
        df = xi*betaMX';
        t = df < -Model;
        df(t) = 0.001*InitialModelBckg(t) - Model(t); % a patch
       	Model = Model + df;              % adjust the f-vector
       	fSum = sum(Model);
       	fChange = sum(df);
			
		ModelScratch= Model/fSum;		% fraction of Model(i) in this bin
		ModelScratch=ModelScratch .* log(ModelScratch);
		sEntropy = sEntropy - sum(ModelScratch);		% from Skilling and Brian eq. 1
		
		zscratch = opus(Model,FormMtrx);
        fitted = zscratch;
		zscratch = (MeasuredData - zscratch) ./ Errors;	%residuals
		Chisquare = sum(zscratch.^2);%new Chisquared
		
		CurrentEntropy = sEntropy;
		CurrentChiSq = Chisquare;
		CurChiSqMinusAlphaEntropy = Chisquare - MaxEntStabilityParam*sEntropy;
%		IR1R_DisplayDiagnostics(CurrentEntropy,CurrentChiSq, CurChiSqMinusAlphaEntropy,iter)		
%       display data in diagnostic graphs, if needed 
%       see, if we have reached solution
		ox = opus(Model,FormMtrx);
		%ox = opus(Model,FormMtrx) + Background';      % BLee add march 21. 2011
        figure(1);
        if iter == 0
            clf;
            subplot(2,1,1);
            loglog(qvec, MeasuredData, 'ro');
            fitlinehandle = line(qvec, ox);
            set(fitlinehandle, 'color', 'b')
            subplot(2,1,2);
            prhandle = plot(x_radius, Model);
        else
            set(fitlinehandle, 'ydata', ox);
            set(prhandle, 'ydata', Model);
        end
        drawnow
%		UpdateGraph()
		if (abs(Chisquare - chizer) < tolerance)
			if (test<tstLim)	% same solution limit
			%solution found
                y = iter;
				return
            end
        end
    end
    y = NaN;
	return
end
    
function model = tropus(data, FormMtrx)
    [m,n]=size(data);
    if m<n
        data = data';
    end
    model = FormMtrx'*data; % I(i) = FormMtrx(j,i) * model(j)
    %model = model';
end

function data = opus(model, FormMtrx)
    [m,n]=size(model);
    if m<n
        model = model';
    end
    data = FormMtrx*model;% * model; % I(i) = FormMtrx(j,i) * model(j)
end

function IR1R_Move(m)
global Chisquare;
global chtarg;
global chizer;
global fSum;
global blank;
global betaMX;
global c1;
global c2;
global s1;
global s2;

MxLoop=500;				%for no solution	
Passes=0.001;			%convergence test

    a1 = 0;                       % lower bracket  "a"
	a2 = 1;      
%    betaMX, c1, c2, s1, s2, 
%    disp('==================A')
	cmin = ChiNow (a1, m);		% get current chi
%    betaMX, c1, c2, s1, s2, 
%    disp('==================B')
    
	ctarg = 0;
	if ((cmin*Chisquare)>chizer) 
		ctarg = 0.5*(1 + cmin);
    end
	if ((cmin*Chisquare) <= chizer) 
		ctarg = chizer/Chisquare;
    end
	f1 = cmin - ctarg;
%	f2 = ChiNow (a2, m, betaMX, c1, c2, s1, s2, flWv, blWv) - ctarg;
	f2 = ChiNow (a2,m);
    f2 = f2 - ctarg;

    for i=0:(MxLoop-1)
		anew = 0.5 * (a1+a2); %          //! choose a new "a"
		fx = ChiNow (anew,m);
        fx = fx - ctarg;
		if (f1*fx >0) 
			a1 = anew;
			f1 = fx;
        end
		if (f2*fx > 0)
			a2 = anew;
			f2 = fx;
        end
		if (abs(fx) < Passes) 
			break
        end
    end

    if (i>=MxLoop-1)
		error(['No convergence in alpha chop (MOVE). Loop counter = ', num2str(MxLoop)])
    end
	w = -1*betaMX*s2*betaMX';

    if (w > 0.1*fSum/blank)
		for k=0:(m-1)
			betaMX(k+1) = betaMX(k+1) * sqrt(0.1 * fSum/(blank * w));
        end
    end
    chtarg = ctarg * Chisquare;
end

        
function res = ChiNow(ax,m)
%    bx = 1-ax;
%	aWv = bx*c2-ax*s2;
%    bWv = -(bx*c1-ax*s1);
%    betaMX = ChoSol(aWv,bWv,m,betaMX, flWv,blWv);
%    z = betaMX*c2';
%    w = sum(betaMX.*c1+0.5*z);
%	res = 1 + w;
global Chisquare;
global chtarg;
global chizer;
global fSum;
global blank;
global betaMX;
global c1;
global c2;
global s1;
global s2;
global blWv;

    aWv = zeros(3);
    bWv = zeros(1,3);
    bx = 1-ax;
    for k=0:(m-1)
        for l=0:(m-1)
            aWv(k+1,l+1) = bx*c2(k+1,l+1)-ax*s2(k+1,l+1);
        end
        bWv(k+1)=-(bx*c1(k+1)-ax*s1(k+1));
    end
    ChoSol(aWv, bWv, m);
    w=0;
    for k=0:(m-1)
        z=0;
        for l=0:(m-1)
            z=z+c2(k+1,l+1)*betaMX(l+1);
        end
        w = w+betaMX(k+1)*(c1(k+1)+0.5*z);
    end
    res = 1+w;
end

function ChoSol(a, b, m)
global betaMX;
global flWv;
global blWv;

	if (a(1,1) <= 0) 
		error(['Fatal error in CHOSOL: a(1,1) = ', num2str(a(1,1))]);
    end
	flWv(1,1) = sqrt(a(1,1));
%	variable i, j, z, k,i1
	for i =1:(m-1)
		flWv(i+1,1) = a(i+1,1) / flWv(1,1);
		for j = 1:i
			z = 0;
			for k = 0:(j-1)
				z = z + flWv(i+1,k+1) * flWv(j+1,k+1);
            end
			z = a(i+1,j+1) - z;
			if (j==i)
				flWv(i+1,j+1) = sqrt(z);
			else
				flWv(i+1,j+1) = z / flWv(j+1,j+1);
            end
        end
    end
	blWv(1) = b(1) / flWv(1,1);
	for i=1:(m-1)
		z = 0;
		for k = 0:(i-1)
			z = z + flWv(i+1,k+1) * blWv(k+1);
        end
%         flWv(i+1,j+1), disp('flWv')
		blWv(i+1) = (b(i+1) - z) / flWv(i+1,i+1);
    end
	betaMX(m) = blWv(m) / flWv(m, m);
	for i1=0:(m-2)
		i = m - 2 - i1;
		z = 0;
		for k = (i+1):(m-1)
			z = z + flWv(k+1,i+1) * betaMX(k+1);
        end
		betaMX(i+1) = (blWv(i+1) - z) / flWv(i+1,i+1);
    end
end

function a = row2vect(a)
[m,n]=size(a);
if m<n
    a = a';
end
end