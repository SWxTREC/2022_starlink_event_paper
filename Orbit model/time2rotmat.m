% Julian Centuries past 1-Jan-2000 12:00 Terrestrial Time
function [rot_ECI2ECEF,dut1,xp,yp] = time2rotmat(eop, time_var,X_state, eqeterms, nut80)
ttt = (time_var + eop.dat + 32.184 - 86400*2451545.0)/86400/36525; % sec
% Interpolate Delta UT1 to TLE time stamps
dut1 = interp1(86400*eop.fmjd,eop.dut1,time_var-86400*2400000.5,'linear'); % still in sec
% Julian day in UT1
jdut1 = (time_var+dut1)/86400;
% Interpolate Earth Orientation parameters (EOP) to TLE time stamps (in UTC)
lod = 1e-3*interp1(86400*eop.fmjd,eop.rlod,time_var-86400*2400000.5,'linear'); % now in sec
xp = (pi/180/3600)*interp1(86400*eop.fmjd,eop.xp,time_var-86400*2400000.5,'linear'); % now in radians (from arc-sec)
yp = (pi/180/3600)*interp1(86400*eop.fmjd,eop.yp,time_var-86400*2400000.5,'linear'); % now in radians (from arc-sec)
dpsi = (pi/180/3600/1000)*interp1(86400*eop.fmjd,eop.dpsi,time_var-86400*2400000.5,'linear'); % now in radians (from milli-arc-sec)
deps = (pi/180/3600/1000)*interp1(86400*eop.fmjd,eop.deps,time_var-86400*2400000.5,'linear'); % now in radians (from milli-arc-sec)
[~,~, ~, rot_ECI2ECEF] = eci2ecef( 1e-3*X_state(1:3), 1e-3*X_state(4:6), [0;0;0], ttt, jdut1, lod, ...
        xp, yp, eqeterms, dpsi, deps, nut80 );