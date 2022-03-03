function stop = BLFit_outfun(x,optimValues,state)
stop = false;
Iq = evalin('base', 'Iq');
% try
%     Iqfit = evalin('base', 'Iqfit');
% catch
%     Iqfit = [];
% end
fit = evalin('base', 'fit');
 switch state
     case 'init'
%         hold on
        disp('Fit starts')
     case 'iter'
     % Concatenate current point and objective function
     % value with history. x must be a row vector.
%       history.fval = [history.fval; optimValues.fval];
%       history.x = [history.x; x];
     % Concatenate current search direction with 
     % searchdir.
%       searchdir = [searchdir;... 
%                    optimValues.searchdirection'];
        if numel(fit.fitlineh) == 1
            if (size(Iq, 2) > 1) && (size(Iq, 1)>size(Iq, 2))
                set(fit.fitlineh, 'ydata', Iq(:,1));
            else
                set(fit.fitlineh, 'ydata', Iq);
            end
            drawnow
        else
%            eidx = 0;
            numpnts = fit.NdataSet;
            for m=1:numpnts
%                if m==1
%                    sidx = 1;
%                else
%                    sidx = eidx + 1;
%                end
%                eidx = eidx+numpnts(m);
%                set(fit.fitlineh(m), 'ydata', Iq(sidx:eidx));
                set(fit.fitlineh(m), 'ydata', Iq{m});
                drawnow
            end
        end

     % Label points with iteration number and add title.
%       text(x(1)+.15,x(2),... 
%            num2str(optimValues.iteration));
%       title('Sequence of Points Computed by fmincon');
     case 'done'
%         hold off
        disp('Fit done')
     otherwise
 end
