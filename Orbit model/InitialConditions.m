%% Initial conditions
%% Starlink simulations
Hp = 250*1e3;           % perigee altitude
Ha = 250*1e3;
% a_sma= Re+H;          % only for circular
r_p = Re+Hp;            % perigee radius
r_a = Re+Ha;
e = (r_a-r_p)/(r_a+r_p);   % eccentricity
a_sma = r_p/(1-e);        % semi major axis from perigee altitude
inc = 0;                  %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 97.8
raan = 0;
w_arg = 40;         % random
true_ano = 0;      % random
u_arg = 0;          % argument of latitude
str_special = 'CE'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None
