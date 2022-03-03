function [fit2dpeak, img, imgc] = fit2Dgaussian(img0, X0, Y0, dXY, backoption, varargin)
% 
%    if nargin < 6
%        fighandle = [];
%    end
fighandle = [];
center = [];
isfit_pos = 1;
center_range_in_pixel = [];
rho_range = 1;
hkls = [];
for i=1:numel(varargin)
    if ~ischar(varargin{i})
        continue
    end
    switch varargin{i}
        case 'fighandle'
            fighandle = varargin{i+1};
        case 'center'
            center = varargin{i+1};
        case 'fix_center'
            isfit_pos = 0;
        case 'center_pix_range'
            center_range_in_pixel = varargin{i+1};
        case 'rho_range'
            rho_range = varargin{i+1};
        case 'hkls'
            hkls = varargin{i+1};
    end
end

    if isempty(center)
        Npeaks = 1;
    else
        Npeaks = size(center, 1);
    end
    
    fit2dpeak.A = [];
    fit2dpeak.X = [];
    fit2dpeak.Y = [];
    x1 = fix(X0-dXY(1)/2):fix(X0+dXY(1)/2);
    y1 = fix(Y0-dXY(2)/2):fix(Y0+dXY(2)/2);
    t = y1 > size(img0, 1) | y1 < 1;
    y1(t) = [];
    t = x1 > size(img0, 2) | x1 < 1;
    x1(t) = [];
    switch backoption
        case 1
            backindex = 7;
        case 2
            backindex = 7:9;
        case 3
            backindex = 7:11;
    end        
    NPparam = 6;
    Nparam = NPparam*Npeaks+numel(backindex);

    % fit parameters
    p = zeros(1, Nparam);
    LowerB = p;
    UpperB = p;
    if ~isempty(x1) && ~isempty(y1)
        [X, Y] = meshgrid(x1, y1);
        img = img0(y1, x1);
        
        Xd = X-X0;
        Yd = Y-Y0;

        for ind = 1:Npeaks
            if isempty(center)
                [maxv, indn] = max(img(:));
                [y0ind, x0ind] = ind2sub(size(img), indn);
                p((ind-1)*NPparam+1) = maxv;
                p((ind-1)*NPparam+2) = Xd(x0ind);
                p((ind-1)*NPparam+3) = Yd(y0ind);
            else
                maxv = interp2(X, Y, img, center(ind, 1), center(ind, 2));
                p((ind-1)*NPparam+1) = maxv;
                p((ind-1)*NPparam+2) = center(ind, 1)-X0;
                p((ind-1)*NPparam+3) = center(ind, 2)-Y0;
            end
            LowerB((ind-1)*NPparam+1) = 1;
            UpperB((ind-1)*NPparam+1) = p((ind-1)*NPparam+1)*1000;
            if isfit_pos
                if isempty(center_range_in_pixel)
                    LowerB((ind-1)*NPparam+2) = -size(img, 2)/4;
                    UpperB((ind-1)*NPparam+2) = size(img, 2)/4;
                    LowerB((ind-1)*NPparam+3) = -size(img, 1)/4;
                    UpperB((ind-1)*NPparam+3) = size(img, 1)/4;
                else
                    LowerB((ind-1)*NPparam+2) = -center_range_in_pixel(1);
                    UpperB((ind-1)*NPparam+2) = center_range_in_pixel(1);
                    LowerB((ind-1)*NPparam+3) = -center_range_in_pixel(2);
                    UpperB((ind-1)*NPparam+3) = center_range_in_pixel(2);
                end
                
            else
                LowerB((ind-1)*NPparam+2) = p((ind-1)*NPparam+2) ;
                UpperB((ind-1)*NPparam+2) = p((ind-1)*NPparam+2) ;
                LowerB((ind-1)*NPparam+3) = p((ind-1)*NPparam+3) ;
                UpperB((ind-1)*NPparam+3) = p((ind-1)*NPparam+3) ;
            end
%             %p((ind-1)*Nparam+2) = 9;
%             LowerB((ind-1)*Nparam+2) = -1*numel(x1);
%             UpperB((ind-1)*Nparam+2) = numel(x1);
%             %p((ind-1)*Nparam+3) = 9;
%             LowerB((ind-1)*Nparam+3) = -1*numel(y1);
%             UpperB((ind-1)*Nparam+3) = numel(y1);
            
            p((ind-1)*NPparam+4) = size(img, 2)/2;
            LowerB((ind-1)*NPparam+4) = 1;
            UpperB((ind-1)*NPparam+4) = size(img, 2);
            p((ind-1)*NPparam+5) = size(img, 1)/2;
            LowerB((ind-1)*NPparam+5) = 1;
            UpperB((ind-1)*NPparam+5) = size(img, 1);
    %         p(6) = 0.1;
    %         LowerB(6) = -1;
    %         UpperB(6) = 1;
            p((ind-1)*NPparam+6) = 0.01;
            LowerB((ind-1)*NPparam+6) = -rho_range;
            UpperB((ind-1)*NPparam+6) = rho_range;
        end
        switch backoption
            case 1
                p(Nparam) = 0.01;
                LowerB(Nparam) = -mean(img(:))/10;
                UpperB(Nparam) = mean(img(:))/3;
                backindex = Nparam;
            case 2
                
                back = mean([img(1, 1), img(1, end), img(end, 1), img(end, end)]);
                dYint = img(end, 1) - img(1,1);
                dXint = img(1, end) - img(1, 1);

                p((Nparam-2):Nparam) = [back, dXint/size(img, 2), dYint/size(img, 2)];
                if back == 0
                    back = 10;
                end
                LowerB((Nparam-2):Nparam) = [-abs(back)-maxv, -abs(p(end-1))*2, -abs(p(end))*2];
                UpperB((Nparam-2):Nparam) = [abs(back)+maxv, abs(p(end-1))*2, abs(p(end))*2];
                
                backindex = (Nparam-2):Nparam;
            case 3
                back = mean([img(1, 1), img(1, end), img(end, 1), img(end, end)]);
                dYint = img(end, 1) - img(1,1);
                dXint = img(1, end) - img(1, 1);
                p((Nparam-4):Nparam) = [back, dXint/size(img, 2), dYint/size(img, 2), 1, 1];
                if back == 0
                    back = 10;
                end
                LowerB((Nparam-4):Nparam) = [-abs(back)-maxv, ...
                    -abs(p(end-3))*10, -abs(p(end-2))*10, -100, -100];
                UpperB((Nparam-4):Nparam) = [abs(back)+maxv, ...
                    abs(p(end-3))*10, abs(p(end-2))*10, 100, 100];

%                 p(Nparam-4) = 0.01;
%                 LowerB(Nparam-4) = -mean(img(:))/10;
%                 UpperB(Nparam-4) = mean(img(:))/3;
%                 p(Nparam-3) = 0.01;
%                 LowerB(Nparam-3) = -10;
%                 UpperB(Nparam-3) = 10;
%                 p(Nparam-2) = 0.01;
%                 LowerB(Nparam-2) = -10;
%                 UpperB(Nparam-2) = 10;
%                 p(Nparam-1) = 0.01;
%                 LowerB(Nparam-1) = -10;
%                 UpperB(Nparam-1) = 10;
%                 p(Nparam) = 0.01;
%                 LowerB(Nparam) = -10;
%                 UpperB(Nparam) = 10;
                backindex = (Nparam-4):Nparam;
        end
        
        options = optimset('fminsearch');
        options = optimset(options, 'TolX',0.0000001);
        options = optimset(options, 'MaxIter',50000);


        fit2dpeak.xarr = Xd;
        fit2dpeak.yarr = Yd;
        disp('Data fitting ... wait a while.')
        INLP = fminsearchcon(@(x) fitwith2dgaussian(x,double(img), ...
            {Xd, Yd, numel(backindex)}),p,LowerB, UpperB, [], [], [], ...
            options);
        isfitdone = 1;
    else
        INLP = ones(1, Nparam)/0;
        isfitdone = 0;
    end

    fit2dpeak.A = [];
    fit2dpeak.X = [];
    fit2dpeak.Y = [];
    fit2dpeak.sigX = [];
    fit2dpeak.sigY = [];
    fit2dpeak.rho = [];
    fit2dpeak.background = [];
    ind = 1;
    while (1)
        fit2dpeak.A = [fit2dpeak.A, INLP((ind-1)*NPparam+1)];
        fit2dpeak.X = [fit2dpeak.X, INLP((ind-1)*NPparam+2)+X0];
        fit2dpeak.Y = [fit2dpeak.Y, INLP((ind-1)*NPparam+3)+Y0];
        fit2dpeak.sigX = [fit2dpeak.sigX, INLP((ind-1)*NPparam+4)];
        fit2dpeak.sigY = [fit2dpeak.sigY, INLP((ind-1)*NPparam+5)];
        fit2dpeak.rho = [fit2dpeak.rho, INLP((ind-1)*NPparam+6)];
        fit2dpeak.background = [INLP(backindex)];
        if numel(INLP) > (ind*NPparam+numel(backindex))
            ind = ind+1;
        else
            break
        end
    end

    fit2dpeak.param = INLP;
    if isfitdone
        [~,imgc]=fitwith2dgaussian(fit2dpeak.param, X(:),{X(:)-X0, Y(:)-Y0, numel(backindex)});
        bkg = background2d(X(:)-X0, Y(:)-Y0, fit2dpeak.background);
        imgc = reshape(imgc, size(X));
        if sum(bkg) < 0
            fprintf('Background for (%s) is negative.\n', hkls); 
        end
        bkg = reshape(bkg, size(X));
    end


    if ~isempty(fighandle) && isfitdone
        figure(fighandle)
        subplot(3,1,1)
        imagesc(img);axis xy;axis image;
        title(sprintf('X = %0.1f, Y = %0.1f. (%s)', X0, Y0, hkls))
        cl = get(gca, 'clim');
        subplot(3,1,2)
        imagesc(imgc);axis xy;axis image;
        title('Fit');
        set(gca, 'clim', cl);
        subplot(3,1,3)
        imagesc(abs(bkg));axis xy;axis image;
        title('Background');
        set(gca, 'clim', cl);
        drawnow;
        pause(0.5);
    end

end