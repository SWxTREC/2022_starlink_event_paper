function yout = ode4(F,step,y0,parameters,Tsamp)
% ODE4  Classical Runge-Kutta ODE solver.
%   yout = ODE4(F,t0,h,tfinal,y0) uses the classical
%   Runge-Kutta method with fixed step size h on the interval
%      t0 <= t <= tfinal
%   to solve
%      dy/dt = F(t,y)
%   with y(t0) = y0.

%   Copyright 2014 - 2015 The MathWorks, Inc.

y = y0;
yout = y';
t0 = step(1);
tfinal = step(end);
h = step(2) - step(1);
kk = 1;
for t = t0 : h : tfinal-h
    kk = kk+1;
    s1 = F(t,y,parameters);
    s2 = F(t+h/2, y+h*s1/2,parameters);
    s3 = F(t+h/2, y+h*s2/2,parameters);
    s4 = F(t+h, y+h*s3,parameters);
    y = y + h*(s1 + 2*s2 + 2*s3 + s4)/6;
    yout(kk,:) = y';
end

yout = yout(1:Tsamp:kk,:);
