%% GMST from IERS 2010
% find gmst from date (radians)
% inputs jdTT_ref = TT centuries from J2000 TT epoch
% Tu = julian UT1  - 2451545 
function gmst = gmst_iers10(jdTT_ref, Tu)

era = mod(2*pi*(0.7790572732640 + 1.00273781191135448*Tu), 2*pi);   % Earth rotation angle in radians
t = jdTT_ref;
gmst = era+ (0.014506 + 4612.156534*t + 1.3915817*t^2 - 0.00000044*t^3 -0.000029956*t^4 - 0.000000036800*t^5)/3600*pi/180;
