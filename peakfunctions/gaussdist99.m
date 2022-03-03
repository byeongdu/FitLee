function [gDist, xx] = gaussdist99(center, width, dmax)
% gaussian distribution generator, 99% 
% gDist, xx will be column vector
if nargin < 3
    xx = eps:1:40;xx = xx/40*(center+abs(width)*3);
    %xx = 0.01:center/100:(center+abs(width)*3);
    if width == 0
        gDist = 1;
        xx = center;
    else
        if abs(width) < center
            gDist = normpdf(xx, center, abs(width))';
            cutoff = find(gDist < 0.03*max(gDist));
            gDist(cutoff) = [];  xx(cutoff) = [];
            xx = xx';
            %gDist = gDist/sum(gDist);
        else
            gDist = -1;
            xx = center;
        end
    end
else
    %xx = 0.01:dmax/50:dmax;
    xx = eps:1:40;xx=xx/40*dmax;
    if width == 0
        gDist = 1;
        xx = center;
    else
            gDist = normpdf(xx, center, abs(width))';
            xx = xx';
            %gDist = gDist/sum(gDist);
    end
end    
    %interXX = 6*width/40;
    %if (center - 3*width) > 0
        %xx = (center-3*width):interXX:(center+3*width);
    %else
        %xx = 1:interXX:(center+3*width);
    %end
    
    %if width > 0
        %gDist = normpdf(xx, center, width)';
    %elseif width == 0
        %gDist = 1;
        %xx = center;
    %elseif width < 0 
        %gDist = normpdf(xx, center, -1*width)';
    %end
