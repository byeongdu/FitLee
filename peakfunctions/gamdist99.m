function [gDist, xx] = gamdist99(center, width)
% gaussian distribution generator, 99% 
% gDist, xx will be column vector
    xx = gaminv(0.001:0.01:0.999, center/abs(width), abs(width));
    if width == 0
        gDist = 1;
        xx = center;
    else
        if abs(width) < center
            gDist = gampdf(xx, center/abs(width), abs(width))';
            xx = xx';
            gDist = gDist/sum(gDist);
        else
            gDist = -1;
            xx = center;
        end
    end