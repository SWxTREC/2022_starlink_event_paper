%% Load data
load GravCoeff_EGM2008_zeroTide
load fes2010_oceantide
Cbar(3,1) = -0.48416948e-3 + 11.6e-12*18.934;
Cbar(4,1) = 0.9571612e-6 + 4.9e-12*18.934;
Cbar(5,1) = 0.5399659e-6 + 4.7e-12*18.934;

%% Furnishing spice kernels to use the spice functions
cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
cspice_furnsh('naif0012.tls')                    %% time leap seconds
cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
cspice_furnsh('pck00010.tpc')
cspice_furnsh('gm_de431.tpc')                    %% GM values
cspice_furnsh('earth_000101_210629_210407.bpc')
%% Force Models to be considered
parameters.sun       = 1;
parameters.moon      = 1;
parameters.srp       = 1;
parameters.drag      = 1;
parameters.tides     = 0;
parameters.relativity= 0;
flag_tides.SolidEarthTides = 0;
flag_tides.OceanTides = 0;
parameters.earth_rad = 0;
parameters.empirical = 0;
%% flags
% vec_est = [1:6, 8:10, 12:14, 18, 20, 31:33, 36:37];
% vec_est = [1:6, 7:10,28:29];
% vec_est = [1:9, 11,12,26,27];
vec_est = [1:6];
vec_meas = [1:6];                  % which measurements to consider
flag_srp = 'Cball';                % Cball - cannonball, Three = three constants, 'Panel'
case_run = 'Truth';                     % 'Truth', 'Estimation'
case_est = 'EKF';
area_factor = 1;
if strcmp(case_run, 'Truth')
    load('spire83_11_7_bff_avg_msis','Af_total_mat','Bf_total_mat','Cf_total_mat','Df_total_mat')
    Af_std = abs(Af_total_mat(:,:,1));
    Bf_std = abs(Bf_total_mat(:,:,1));
    Cf_std = abs(Cf_total_mat(:,:,1));
    Df_std = abs(Df_total_mat(:,:,1));
    K_ind = 3;                           % true drag-coeff selection
    estimated_coeff.Cr = 0;
    estimated_coeff.Cd = 0;
    estimated_coeff.rho_DMC = 0;
    estimated_coeff.CdTrue = 0;
    estimated_coeff.CdDiscrete = 0;
    estimated_coeff.CdDiscrete_est = 0;
    estimated_coeff.Cr_erp = 0;
    flag_rho = 'MSIS00';                  % 'MSIS00', 'JB08', 'Exp'
    flag_earthrad = 'Knocke';
    Kl = 0;    %1.44e6                      % Mehta and walker .  % . Sesam 5e6/133.322;
    del_r = [0;0;0];     % initial standard deviations
    del_v = [0;0;0];
    parameters.roll = 0;
    parameters.pitch = 0;
    parameters.yaw = 0;
else
    load('spire83_11_7_bff_avg_msis','Af_total_mat','Bf_total_mat','Cf_total_mat','Df_total_mat')
    estimated_coeff.Cr = 1;
    estimated_coeff.Cd = 0;
    estimated_coeff.rho_DMC = 1;
    estimated_coeff.CdTrue = 1;
    estimated_coeff.CdDiscrete = 0;
    estimated_coeff.CdDiscrete_est = 0;
    estimated_coeff.Cr_erp = 0;
    K_ind = 3;                           % true drag-coeff selection
    flag_rho = 'MSIS00';                  % 'MSIS00', 'JB08', 'Exp'
    flag_earthrad = 'Knocke';
    Kl = 5e7;    %1.44e6                      % Mehta and walker .  % . Sesam 5e6/133.
    Af_std = 0.1*ones(31,1);% abs(Af_total_mat(:,:,K_ind));
    Bf_std = 0.1*ones(31,1);% abs(Bf_total_mat(:,:,K_ind));
    Cf_std = 0.1*ones(31,1);% abs(Cf_total_mat(:,:,K_ind));
    Df_std = 0.1*ones(31,1);% abs(Df_total_mat(:,:,K_ind));
    
%     parameters.roll = 0;
%     parameters.pitch = 20;
%     parameters.yaw = 0;
end
flag_drag = 'Bod';                 % body model: Bod, orbit model: Orb; Body orbit model: Bod_orb
bo_mod = 'BOS';                   % BODF : don't forget to change both otheta and ophi
Order_b = 0;
Order_o = 0;
%% Planetary constants
c_light = 299792458.000000;                             %% Speed of light  [m/s]; DE430
mu_e    = cspice_bodvrd('EARTH', 'GM',3)*1e9;           %% mu of Earth in m3/s2 % 0.3986004415E+15; %  doesn't make a difference
mu_s    = cspice_bodvrd('SUN', 'GM',3)*1e9;             %% mu of Sun in m3/s2
mu_m    = cspice_bodvrd('MOON', 'GM',3)*1e9;            %% mu of Moon in m3/s2
Rad_e   = cspice_bodvrd('EARTH', 'RADII',3);            %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Rad_s   = cspice_bodvrd('SUN', 'RADII',3);              %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Re      = Rad_e(1)*1e3; %  6378136.3; %  % doesn't make a difference
Rs      = Rad_s(1)*1e3;
Psun    = 4.56e-6;                                      %% Radiation pressure at 1 au (Doornbos PhD) N/m2
AU      = 149597870660;                                 %% AU value m
omega_e = 7.2921150e-5;                                 %% angular velocity of Earth in rad/s
GPS_sec_offset = 19;                                    %% TAI to GPST offset (TAI-GPST)
deg_grav = 10;                                         %% degree of spherical harmonics
ord_grav = 10;                                          %% order of gravity harmonics

Tides_coefficients
parameters.coeff0 = coeff0;
parameters.coeff1 = coeff1;
parameters.coeff2 = coeff2;


%% Time definition
yyyy = 2003; %2021;%  
mon  =  02;%10; %10; %  %  1 ; % 
day_mon  = 02; % 29; %29; 04; %% 31; % 
n_days = 1;
del_T = 10;
T_total = n_days*86400;
t_sim = 0:del_T:T_total;

load nut80.dat
eop = read_finals_data;
eop.dat = read_tai_utc_dat; % reads only the most current value in the local data file
eqeterms = 1;


jdutc_sec = 86400*GREGORIANtoJD_vector(yyyy,mon,day_mon) + t_sim; %%%%%%%%%%%%
time_prop = jdutc_sec;

time_prop_utc = time_prop - time_prop(1);  %% use time_prop  to calculate theta
time_prop_utc_ekf = time_prop_utc;

jd_epoch = 86400*GREGORIANtoJD_vector(2000,1,1) + 3600*12 + 60*0 + 0;
tai_time_gps  = jdutc_sec -jd_epoch + eop.dat;
et_time_gps   = cspice_unitim(tai_time_gps, 'TAI', 'TDB');
%% Polar variables
[xp_mean, yp_mean] = mean_polar(yyyy(1)); % average in arc-secs
Cbar(3,2) = (sqrt(3)*xp_mean*Cbar(3,1) - xp_mean*Cbar(3,3) + yp_mean*Sbar(3,3))/3600*pi/180;
Sbar(3,2) = (-sqrt(3)*yp_mean*Cbar(3,1) - yp_mean*Cbar(3,3) - xp_mean*Sbar(3,3))/3600*pi/180;

%% Sun and Moon positions
[sun_pos, ~]  = cspice_spkpos('Sun', et_time_gps, 'J2000', 'NONE', 'EARTH');          %% positions of sun and moon in km
sun_pos = sun_pos*1e3;
[moon_pos, ~] = cspice_spkpos('Moon', et_time_gps, 'J2000', 'NONE', 'EARTH');
moon_pos = moon_pos*1e3;
[earth_state, ~] = cspice_spkezr('Earth', et_time_gps, 'J2000', 'NONE', 'Sun');
earth_state = earth_state*1e3;
earth_vel = earth_state(4:6,:);

%% GPM Parameters data for day 180 (June 29 2017)
mass = 300;                                       %% in kg before the maneuver of Aug 9
area = 1;                                         %% change it for order 0
srp_area = 1.4;                                   %% reference area GPM, might need to change later
sun_area_mass = srp_area/mass;                    %% reference area for solar radiation  pressure

epoch = 0;                                    %% no. of seconds since the beginning of the day in UTC for the first observation (nrlmsise-00)
td = datetime(yyyy,mon,day_mon);
doy =  day(td, 'dayofyear');  % 302; %                     %% day of year (nrlmsise-00)
year = yyyy(1);   %  2003; %                          %% year (nrlmsise00)

eps  = 1000;                                   %% altitude diff. in m
days_prev = 0;
%% DMC model parameters
X_s = [0;0;0];
Xs_std = [0.3, 0.03, 1]; %[1e-14, 1e-14,1e-14];
zeta =  0.5; %0.0093; %0.2019; % 0.01; %
omega = 0.0011; %0.0018; %2*pi/T_orb; % 2*2*pi/T_orb;
tau_inv = 1.77e-4;%1.2e-4; %1/100;
sig_dmc = sqrt(5.5e-22);
q_dmc = 1;
c1 = 0; %1e-17; % 1e-17; %sqrt(4*omega^3*zeta*sig_dmc^2/q_dmc);
c2 = 1e-8; %1e-18; %1.1e-17; %1e-16; %1e-16;
c3 = 1e-5; %1e-17; %1.1e-14; %1e-17;

% cd dmc parameters
X_cd = 0;
Xcd_std = 0.5;
tau_inv_cd = 0; %T_orb; %1/100;
c_cd = 0; % 0.5; %0.5e-3; %1e-3;
%% SRP/drag Parameters
if estimated_coeff.CdDiscrete_est
    Cd_nom = 1;
else
    Cd_nom = 6; %2.2*(0.681+0.426); %2.2*(2.77*8.05) %0.3648;
end
if strcmp(flag_srp, 'Cball') || strcmp(flag_srp, 'Panel') 
    Cr = 1; % 0.41;
    Cr_std = 0.5;
elseif strcmp(flag_srp, 'Three') 
    Cr = 1;
    A0 = -4.0;
    A1 = -1.6;
    A2 = -0.5;
    A0_std = 0.5;
    A1_std = 0.5;
    A2_std = 0.5;    
else
    Cr = 1;
end
Cr_erp = 1;
Cr_erp_std = 1;
%% Body properties for drag and SRP
R = 8314.47215 ;                         % Universal gas constant (SI)
Tw = 300;                             % Temperature of the satellite wall
Alpha = 0.93;
frac = 1;
area_vec(:,1) = [1;0;0];
area_vec(:,2) = [-1;0;0];
area_vec(:,3) = [0;1;0];
area_vec(:,4) = [0;-1;0];
area_vec(:,5) = [0;0;1];
area_vec(:,6) = [0;0;-1];
area_vec(:,7) = [ 0;1; 0];               % front solar panel half
area_vec(:,8) = [ 0;-1; 0];               % front solar panel half

Area_plates = [0.21, 0.21, 0.681+0.426, 0.681+0.426, ...
    1.27*2.82, 1.27*2.82, 2.77*8.05, 2.77*8.05]/area;
Ar = 1;

spec_ref = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05, 0.2];
diff_ref = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.4];
N_plates = numel(spec_ref);
flag_axis = 'x';

M_s = [26.981538, 26.981538, 26.981538, 26.981538, 26.981538, 26.981538, 60.0843, 26.981538];             % asymmetry test case
% M_s = [100];             % asymmetry test case
Ms = mean(M_s);
atm_mass = [4.002602,15.9994,28.0134,31.9988,39.948,1.0079,14.0067, 15.9994]; % amu's or molar mass:g/mol; He, O, N2, O2, Ar, H, N, O
amu = 1.66053904e-27; % kg
k_b = 1.38064852e-23; % m2 kg s-2 K-1
shape_model = 'plate_dria_multi';

angle_panels = zeros(numel(time_prop_utc),1,N_plates);
%% Angles for drag
% theta = repmat(theta,N_plates,1);
% phi = zeros(size(theta));
%% Measurement noise values, actual GPM values
sig_meas(1) = 10; %3;  ; found by using the evar function on the residuals of the real data for order 150
sig_meas(2) = 1e-1; % 0.01; %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% noise mismatch
sigma_pos = sig_meas(1);
sigma_vel = sig_meas(2);
R_aug = diag([sigma_pos,sigma_pos,sigma_pos,sigma_vel,sigma_vel,sigma_vel].^2);
%% Initial conditions
InitialConditions

[r_eci_true, v_eci_true] = coe2rv(a_sma,e,inc,raan,w_arg,true_ano,mu_e,str_special, u_arg);


if e > 1e-6
    str_special = 'NO';
else
    str_special = 'CI';
end
Fwind = (1-r_p*omega_e/norm(v_eci_true(:,1))*cosd(inc))^2;
r_circ = a_sma;
v_circ = sqrt(mu_e/r_circ);
n_mean = v_circ/r_circ;
T_orb = 2*pi/n_mean;
%% Initialize state

X_init = [r_eci_true(:,1);v_eci_true(:,1)] + [del_r;del_v];

X_ref_st(:,1) = X_init;
% y_meas = [reci;veci];     % measurements
% y_meas = y_meas(vec_meas,:);
%% Density calculations
g_acc = mu_e/r_p^2;          % gravity acceleration
he  = atm_mass(1);
oxa = atm_mass(2);
nitm= atm_mass(3);
oxm = atm_mass(4);
ar  = atm_mass(5);
h   = atm_mass(6);
nita= atm_mass(7);

%% Density model
switch flag_rho
    case 'Exp'
        g_acc = mu_e/r_p^2;          % gravity acceleration
        [T_atm,n_den] = atmosnrlmsise00(Hp, 0, 0, 2018, 1, 0, 69, 69, 4*ones(1,7),'Oxygen');
        rho0     = n_den(:,6);
        Talt = T_atm(2);
        he  = atm_mass(1);                        %%%%% need to change to amu
        oxa = atm_mass(2);
        nitm= atm_mass(3);
        oxm = atm_mass(4);
        ar  = atm_mass(5);
        h   = atm_mass(6);
        nita= atm_mass(7);
        rho0_he     = n_den(:,1)*he*amu;
        rho0_oxa     = n_den(:,2)*oxa*amu;
        rho0_nitm     = n_den(:,3)*nitm*amu;
        rho0_oxm     = n_den(:,4)*oxm*amu;
        rho0_ar     = n_den(:,5)*ar*amu;
        rho0_h     = n_den(:,7)*h*amu;
        rho0_nita     = n_den(:,8)*nita*amu;
        rho0_oxa_ano = n_den(:,9)*oxa*amu;
        rho0_all = [rho0_he,rho0_oxa,rho0_nitm,rho0_oxm,rho0_ar,rho0_h,rho0_nita,rho0_oxa_ano];
        
        H_he = k_b*Talt/(he*amu*g_acc);
        H_oxa = k_b*Talt/(oxa*amu*g_acc);
        H_nitm = k_b*Talt/(nitm*amu*g_acc);
        H_oxm = k_b*Talt/(oxm*amu*g_acc);
        H_ar = k_b*Talt/(ar*amu*g_acc);
        H_h = k_b*Talt/(h*amu*g_acc);
        H_nita = k_b*Talt/(nita*amu*g_acc);
        H_scale_all = [H_he,H_oxa,H_nitm,H_oxm,H_ar,H_h,H_nita,H_oxa];
        mass_sum = n_den(:,1)*he + n_den(:,2)*oxa + n_den(:,3)*nitm + n_den(:,4)*oxm + n_den(:,5)*ar + n_den(:,7)*h + n_den(:,8)*nita + n_den(:,9)*oxa;
        num = sum(n_den,2);
        M_mean = mass_sum./num*amu;
        
        H_scale = k_b*Talt/(M_mean*g_acc);
        parameters.rho0 = rho0;
        parameters.H_scale = H_scale;
        parameters.r0 = r_p;
        parameters.rho0_all = rho0_all;
        parameters.H_scale_all = H_scale_all;
        parameters.Talt = Talt;
    case 'MSIS00'
        load SOLFSMY
        load SOLRESAP
        [Ap_total, Ap_daily, F10_total, F81c] = msis_inputs(year, SOLdata, geomag_mat);
        parameters.Ap_total = Ap_total;
        parameters.Ap_daily = Ap_daily;
        parameters.F10_total = F10_total;
        parameters.F81c = F81c;
        parameters.days_year_prev = 0;
    case 'JB08'
        load SOLFSMY
        load DTCFILE
        earth_model = wgs84Ellipsoid;                  %% shape model for geodetic to geocentric latitude
        flattening = earth_model.Flattening;           %% flattening factor
        ind_sol = find(SOLdata(1,:) == year,1);        %% index to point towards data for doy
        ind_mag = find(DTCdata(1,:) == year,1);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        parameters.flattening = flattening;
        parameters.ind_sol = ind_sol;
        parameters.ind_mag = ind_mag;
        parameters.SOLdata = SOLdata;
        parameters.DTCdata = DTCdata;
    case 'HASDM'
        load 2003_HASDM_500-575KM
        ind_day = [];
        for nn = 1:n_days
            doy_hasdm = doy+nn-1;
            ind_day  = [ind_day;find(den_mat(:,2) == doy_hasdm)];
        end
        hasdm_mat = den_mat(ind_day,:);

%         ind_alt1 = find(den_mat(:,7) == 250);
%         ind_alt2 = find(den_mat(:,7) == 225);
%         ind_day_alt1 = intersect(ind_day, ind_alt1);
%         ind_day_alt2 = intersect(ind_day, ind_alt2);
%         hasdm_mat = [den_mat(ind_day_alt1,:);den_mat(ind_day_alt2,:)];
%         
% %         load 2007_HASDM_175-275KM
%         ind_day = [];
%         for nn = 1:n_days
%             doy_hasdm = doy+nn-1;
%             ind_day  = [ind_day;find(den_mat(:,2) == doy_hasdm)];
%         end
%         ind_alt1 = find(den_mat(:,7) == 275);
%         ind_day_alt1 = intersect(ind_day, ind_alt1);
%         hasdm_mat = [hasdm_mat;den_mat(ind_day_alt1,:)];
        
        clear den_mat
        vec_alt  = unique(hasdm_mat(:,7));
        vec_njd  = unique(hasdm_mat(:,6));
        vec_long = unique(hasdm_mat(:,10));
        vec_lat  = unique(hasdm_mat(:,9));
        hasdm_mat(:,11) = log(hasdm_mat(:,11));
        ind_long = [1:numel(vec_long)];
        den_grid = zeros(numel(vec_alt), numel(vec_njd),numel(vec_long), numel(vec_lat));
        
        N_njd = numel(vec_njd);
        N_long = numel(vec_long);
        N_lat = numel(vec_lat);
        N_alt = numel(vec_alt);
        for ii_alt = 1:N_alt
            ind_alt = find(hasdm_mat(:,7) == vec_alt(ii_alt));
            for ii_njd = 1:N_njd
                ind_njd = find(hasdm_mat(:,6) == vec_njd(ii_njd));
                for ii_lat = 1:N_lat
                    ind_alt_njd = intersect(ind_alt,ind_njd);
                    hasdm_small = hasdm_mat(ind_alt_njd,:);
                    
                    ind_lat = find(hasdm_small(:,9) == vec_lat(ii_lat));
                    hasdm_small = hasdm_small(ind_lat,:);
                    
                    den_ind = interp1(hasdm_small(:,10), hasdm_small(:,11), vec_long,'linear', 'extrap');
                    den_grid(ii_alt,ii_njd,ind_long,ii_lat) = den_ind;
                end
            end
        end
        gridVecs = {vec_alt,vec_njd,vec_long,vec_lat};
        parameters.F_hasdm = griddedInterpolant(gridVecs,den_grid);
        
        
        parameters.jd_ref = 86400*GREGORIANtoJD_vector(1950,1,0) + 3600*0 + 60*0 + 0;
end

%% Initialize constants
parameters.e = e;
parameters.a_sma = a_sma;
parameters.inc = inc;
parameters.raan = raan;
parameters.w_arg = w_arg;
parameters.true_ano = true_ano;
parameters.u_arg = u_arg;
parameters.str_special = str_special;
parameters.Fwind = Fwind;
parameters.Re = Re;
parameters.mu_e = mu_e;
parameters.R = R;
parameters.Kl = Kl;
parameters.atm_mass = atm_mass;
parameters.M_s = M_s;
parameters.Ar = Ar;
parameters.Area_plates = Area_plates;
parameters.Tw = Tw;
parameters.Alpha = Alpha;
parameters.Order_b = Order_b;
parameters.Order_o = Order_o;
parameters.flag_axis = flag_axis;
parameters.area_vec = area_vec;
parameters.amu = amu;
parameters.sun_pos = sun_pos;
parameters.epoch = epoch;
parameters.doy = doy;
parameters.year = year;
parameters.eps = eps;
parameters.n_mean = n_mean;
parameters.flag_rho = flag_rho;
parameters.k_b = k_b;
parameters.shape_model = shape_model;
parameters.omega_e = omega_e;
% parameters.rot_SBF2ECI = rot_SBF2ECI;
% parameters.rot_ECI2ECEF = rot_ECI2ECEF;
parameters.theta = zeros(N_plates,numel(time_prop_utc_ekf));  %repmat(theta, N_plates,1);
% parameters.angle_panels = angle_panels;
parameters.phi = zeros(N_plates,numel(time_prop_utc_ekf)); %repmat(phi, N_plates,1);
parameters.eqeterms = eqeterms;
parameters.eop = eop;
parameters.nut80 = nut80;
parameters.srp_area = srp_area;
parameters.mass  = mass;
parameters.deg_grav = deg_grav;
parameters.ord_grav = ord_grav;
parameters.P_sun = Psun;
parameters.R_sun = Rs;
parameters.rho_diff = diff_ref;
parameters.rho_spec = spec_ref;
parameters.drag_area = area;
parameters.flag_srp = flag_srp;
parameters.mu_sun  = mu_s;
parameters.mu_moon  = mu_m;
parameters.AU = AU;
parameters.c = c_light;
parameters.Cd_discrete = 0;
parameters.time_et = et_time_gps;
parameters.flag_tides = flag_tides;
parameters.doo_coeff = doo_coeff;
parameters.Cnmp_ot = Cnmp_ot;
parameters.Snmp_ot = Snmp_ot;
parameters.Cnmm_ot = Cnmm_ot;
parameters.Snmm_ot = Snmm_ot;
parameters.deg_do = deg_do;
parameters.ord_do = ord_do;
parameters.pole_mean =  [xp_mean, yp_mean];
parameters.case_est = case_est;
parameters.c_cr = 0;
parameters.altitude = 0;
parameters.latitude = 0;
parameters.time_jd = 0;
parameters.longitude = 0;
parameters.flag_earthrad = flag_earthrad;
% parameters.lat = lat';
% parameters.lon = lon';
% parameters.Psw = sw_flux/c_light;
% parameters.Plw = lw_flux/c_light;
parameters.Cr_erp = Cr_erp;
parameters.frac = frac;