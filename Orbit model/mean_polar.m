%% Polar wobble variables
% iers 2010
function [xp_mean, yp_mean] = mean_polar(t)
t0 = 2000;
if t < 2010
    xp_mean = 55.974 + (t-t0)*1.8243 + (t-t0)^2*0.18413 + (t-t0)^3*0.007024;
    yp_mean = 346.346 + (t-t0)*1.7896 - (t-t0)^2*0.10729 - (t-t0)^3*0.000908;
else
    xp_mean = 23.513 + (t-t0)*7.6141 ;
    yp_mean = 358.891 - (t-t0)*0.6287 ;
end

xp_mean = xp_mean/1000; % arc-seconds
yp_mean = yp_mean/1000;