function [minidx,maxidx] = getScenarioIndex(scenario,N)
% return start and end indexes for observed variables

umin = 1;
umax = N/3;

vmin=N/3 + 1;
vmax = 2*N/3;

hmin=2*N/3 + 1;
hmax = N;

if scenario == 1      
    minidx = umin;
    maxidx = vmax;
elseif scenario == 2
    minidx = 1;
    maxidx = N;
else
    minidx = hmin;
    maxidx = hmax;
end

end