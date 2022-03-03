function BLFit_drawdistribution
fit = evalin('base', 'fit');
mn = [get(gcbf, 'tag'), 'distr'];
mf = findobj('tag', mn);
if isempty(mf)
    mf =figure;
    set(mf, 'tag', mn);
end
plot(fit.distr(:,1), fit.distr(:,2:end));
title('Number distribution')
xlabel('radius (A)')