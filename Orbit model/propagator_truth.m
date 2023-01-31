%% propagator
function [Xdot,Cd_est,rho, theta_v, phi_v,a_grav,a_sun,a_moon,a_drag,a_srp, a_earthrad] = propagator_truth(t,X,parameters)
% Cd_est,rho,theta_in,phi_in,bff_coeff,S,m_r,frac,r_ads
%% Parameters
N_st = parameters.N_st;
N_nom = parameters.N_nom;
Order_o = parameters.Order_o;
Order_b = parameters.Order_b;
Af_ind = parameters.Af_ind;
Bf_ind = parameters.Bf_ind;
Cf_ind = parameters.Cf_ind;
Df_ind = parameters.Df_ind;
time_prop = parameters.time_prop;  % et time is cspice is used, jd if vallado is used
t_et  = parameters.time_et;
tym = parameters.tym;
omega_e = parameters.omega_e;
epoch = parameters.epoch;
doy = parameters.doy;
eps = parameters.eps;
flag_srp = parameters.flag_srp;
earth_velocity = parameters.earth_vel;
estimated_coeff = parameters.estimated_coeff;
omega = parameters.omega;
zeta = parameters.zeta;
tau_inv = parameters.tau_inv;
tau_inv_cd = parameters.tau_inv_cd;
c1 = parameters.c1;
c2 = parameters.c2;
c3 = parameters.c3;
c_cd = parameters.c_cd;
Cd_est = parameters.Cd_est;
Cnm = parameters.Cnm;
Snm = parameters.Snm;
mu_e = parameters.mu_e;
% theta_max = parameters.theta_max;
theta = parameters.theta;
Cr_est = parameters.Cr_est;
eop = parameters.eop;
roll = parameters.roll;
pitch = parameters.pitch;
yaw = parameters.yaw;

X_state = X(1:6);


%% Interpolation


[~,t_ind] = min(abs(tym - t));
% if strcmp(flag_drag,'Bod')
%     %           phi_in = Input(2,t_ind);
%     %           theta_in = Input(3,t_ind);
%     mu_e = Force(1).Model.mu_e;
%     coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
%     theta_in = coe(7); % for circulr orbits, take the argument of latitude
%     phi_in = coe(7);
% end
% if strcmp(flag_drag,'Bod_orb')
%
%     mu_e = Force(1).Model.mu_e;
%     coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
%     theta_in = coe(7); % for circulr orbits, take the argument of latitude
%     phi_in = coe(7);
% end

if t_ind ~= numel(tym)
    if tym(t_ind) ~= tym(t_ind+1)
        t2 = tym(t_ind+1);
        
    else
        t2 = tym(t_ind+2);
    end
    t1 = tym(t_ind);
    time_var = (t-t2)*(time_prop(t_ind+1) - time_prop(t_ind))/(t2-t1) + time_prop(t_ind+1);
    time_et = (t-t2)*(t_et(t_ind+1) - t_et(t_ind))/(t2-t1) + t_et(t_ind+1);
    ind_theta = (theta(:,t_ind+1)./abs(theta(:,t_ind+1))).*(theta(:,t_ind)./abs(theta(:,t_ind))) > 0;
    theta_in(ind_theta) = (t-t2)*(theta(ind_theta,t_ind+1) - theta(ind_theta,t_ind))/(t2-t1) + theta(ind_theta,t_ind+1);
    theta_in(~ind_theta) = theta(~ind_theta,t_ind);
    coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
    phi_in = coe(7);
else
    t2 = tym(t_ind);
    if tym(t_ind) ~= tym(t_ind-1)
        t1 = tym(t_ind-1);
        
    else
        t1 = tym(t_ind-2);
    end
    time_var = (t-t2)*(time_prop(t_ind) - time_prop(t_ind-1))/(t2-t1) + time_prop(t_ind);
    time_et = (t-t2)*(t_et(t_ind) - t_et(t_ind-1))/(t2-t1) + t_et(t_ind);
    ind_theta = (theta(:,t_ind)./abs(theta(:,t_ind))).*(theta(:,t_ind-1)./abs(theta(:,t_ind-1))) > 0;
    theta_in(ind_theta) = (t-t2)*(theta(ind_theta,t_ind) - theta(ind_theta,t_ind-1))/(t2-t1) + theta(ind_theta,t_ind);
    theta_in(~ind_theta) = theta(~ind_theta,t_ind);
    coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
    phi_in = coe(7);
end



X_f = X(N_nom+1:N_st,1);
% if strcmp(flag_drag,'Orb') == 1
%     %%%%%%%%%%%%%%%%% Orbit Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     theta_in = 0;
%     mu_e = Force(1).Model.mu_e;
%     coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
%     phi_in = coe(7);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
earth_vel = earth_velocity(:,t_ind);

%% Drag coefficient calculation

% %%%%%%
% if theta_in < 0
%     theta_in = -theta_in;
% end
% theta_in = theta_in/theta_max;
% %%%%%%



%% ECEF to ECI rotation matrix
[rot_ECI2ECEF,dut1,xp,yp] = time2rotmat(eop, time_var,X_state, parameters.eqeterms, parameters.nut80); % parameters.rot_ECI2ECEF(:,:,t_ind);
% rot_ECI2ECEF = cspice_pxform( 'J2000', 'ITRF93', time_var );
parameters.TEI = rot_ECI2ECEF;

%% ECI TO NTW rotation matrix
rhat = X_state(1:3)/norm(X_state(1:3));
vhat = X_state(4:6)/norm(X_state(4:6));
That = vhat;
What = cross(vhat,rhat);
What = What/norm(What);
Nhat = cross(That,What);
Nhat = Nhat/norm(Nhat);
ECI2NTW = [ Nhat'; That'; What'];
theta_v = 0;
phi_v = 0;
Xhat = rhat;
Zhat = cross(rhat,vhat);
Yhat = cross(Zhat,Xhat);
Yhat = Yhat/norm(Yhat);
ECI2RIC = [Xhat';Yhat';Zhat'];
%% ECI TO BODY
% x_hat_body = parameters.sun_pos(:,t_ind)/norm(parameters.sun_pos(:,t_ind));
% h_ang = cross(X_state(1:3)/norm(X_state(1:3)), X_state(4:6)/norm(X_state(4:6)));
% h_ang = h_ang/norm(h_ang);
% z_hat_body = cross(x_hat_body,h_ang);
% z_hat_body = z_hat_body/norm(z_hat_body);
% y_hat_body = cross(z_hat_body,x_hat_body);
% y_hat_body = y_hat_body/norm(y_hat_body);
% SBF2ECI = [x_hat_body,y_hat_body,z_hat_body];
% ECI2SBF = SBF2ECI';
%% ECI TO LVLH
z_hat = -X_state(1:3)/norm(X_state(1:3));   % nadir
y_hat = cross(X_state(4:6), X_state(1:3));  % negative orbit normal
y_hat = y_hat/norm(y_hat);
x_hat = cross(y_hat, z_hat);
x_hat = x_hat/norm(x_hat);

ECI2LVLH = [x_hat,y_hat,z_hat]';

%% LVLH2SBF
roll_mat = [1,0,0; 0,cosd(roll),sind(roll); 0,-sind(roll),cosd(roll)];
pitch_mat = [cosd(pitch), 0, -sind(pitch); 0,1,0; sind(pitch), 0, cosd(pitch)];
yaw_mat = [cosd(yaw), sind(yaw), 0; -sind(yaw),cosd(yaw),0;0,0,1];
LVLH2SBF = roll_mat*pitch_mat*yaw_mat;

ECI2SBF = LVLH2SBF*ECI2LVLH;
SBF2ECI = ECI2SBF';

v_eci = X_state(4:6) - cross([0;0;omega_e],X_state(1:3));
v_sbf = ECI2SBF*v_eci;
v_sbf = v_sbf/norm(v_sbf);
theta_v = atan2(v_sbf(2),v_sbf(1));
phi_v = atan2(v_sbf(3),sqrt(v_sbf(1)^2+v_sbf(2)^2));


% [phi_v, theta_v, SBF2ECI] = theta_iner(X_state,parameters.sun_pos(:,t_ind),parameters.R_sun,parameters.Re,omega_e);
SBF2ECEF = rot_ECI2ECEF*SBF2ECI;
n_hat = SBF2ECI*parameters.area_vec;           % body to eci of plate vectors
%% Force models
N_est = 6; %numel(vec_est);
Ncr = 0;
Ba = zeros(6,1);
Fpar_dot = [];
Xpar_dot = [];


if parameters.tides
    
    tai_sec = time_var + eop.dat;
    tt_sec = tai_sec + 32.184;
    jdTT_ref = (tt_sec/86400 - 2451545)/36525;
    jd_ut1 = (time_var + dut1)/86400 - 2451545;
    % [Cnm, Snm] = tides(parameters.Cnm, parameters.Snm, time_et, jdTT_ref, jd_ut1, rot_ECI2ECEF, parameters.moon_pos(:,t_ind), parameters.sun_pos(:,t_ind), ...
    %     parameters.Re, mu_e,parameters.mu_moon, parameters.mu_sun, xp, yp, parameters.flag_tides, parameters.coeff0, parameters.coeff1, parameters.coeff2);
    [Cnm, Snm] = tides_fes(parameters.Cnm, parameters.Snm, time_et, jdTT_ref, jd_ut1, rot_ECI2ECEF, parameters.moon_pos(:,t_ind), parameters.sun_pos(:,t_ind),...
        parameters.Re, mu_e,parameters.mu_moon, parameters.mu_sun, xp, yp, parameters.flag_tides, parameters.coeff0, parameters.coeff1, parameters.coeff2, ...
        parameters.doo_coeff, parameters.Cnmp_ot, parameters.Snmp_ot, parameters.Cnmm_ot, parameters.Snmm_ot, parameters.deg_do, parameters.ord_do, parameters.pole_mean);
    
end


Fmod          = HighGrav1(parameters.deg_grav, parameters.ord_grav, mu_e, parameters.Re, Cnm, Snm, rot_ECI2ECEF);
[accel,Fpos,Fvel,Fa_par]  = Fmod.compAccelPar(X_state,t);
a_grav = accel;
if parameters.sun
    Fmod      = ThreeBody(parameters.mu_sun, parameters.sun_pos(:,t_ind));
    [acc_new,Fpos_new,Fvel_new,~]  = Fmod.compAccelPar(X_state,t);
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
    a_sun = acc_new;
else
    a_sun = [0;0;0];
end

if parameters.moon
    Fmod      = ThreeBody(parameters.mu_moon, parameters.moon_pos(:,t_ind));
    [acc_new,Fpos_new,Fvel_new,~]  = Fmod.compAccelPar(X_state,t);
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
    a_moon = acc_new;
else
    a_moon = [0;0;0];
end

if parameters.srp
    SRP_am    = parameters.srp_area/parameters.mass;
    if strcmp(flag_srp, 'Cball')
        if estimated_coeff.Cr
            Ncr = 1;
            Cr_est = X(7,1);
        end
        Fmod      = SRP(SRP_am, parameters.AU, Cr_est, parameters.P_sun, parameters.sun_pos(:,t_ind), parameters.Re, parameters.R_sun);
        fs = Fmod.shadow_func(X_state,t);
    elseif strcmp(flag_srp, 'Three')
        if estimated_coeff.Cr
            Ncr = 3;
            A0 = X(7,1);
            A1 = X(8,1);
            A2 = X(9,1);
        end
        Fmod      = SRP_3constants(SRP_am, parameters.AU, Cr_est, parameters.P_sun, parameters.sun_pos(:,t_ind), parameters.Re, parameters.R_sun, earth_vel,A0,A1,A2);
    elseif strcmp(flag_srp, 'Panel')
        Ncr = 1;
        % %         ang_pan = parameters.angle_panels(t_ind,:,:);
        % %         zero_mat = zeros(1,1,numel(ang_pan));
        % %         one_mat = ones(1,1,numel(ang_pan));
        % %         rot_pan2bod = [cosd(ang_pan), zero_mat, -sind(ang_pan); zero_mat,one_mat,zero_mat; sind(ang_pan), zero_mat, cosd(ang_pan)];
        % %         area_vec = reshape(parameters.area_vec,3,1,9);
        % %         area_vec_b = squeeze(pagemtimes(rot_pan2bod,area_vec));
        % %         SBF2ECI = parameters.rot_SBF2ECI(:,:,t_ind);
        % %         n_hat = SBF2ECI*area_vec_b;
        
        %         SBF2ECI = parameters.rot_SBF2ECI(:,:,t_ind);
        
        if estimated_coeff.Cr
            Cr_est = X(7,1);
        end
        Fmod      = SRP_GPM(parameters.mass, parameters.AU, parameters.rho_spec, parameters.rho_diff,parameters.P_sun, parameters.sun_pos(:,t_ind), parameters.Re,...
            parameters.R_sun, n_hat,parameters.Area_plates,numel(parameters.Area_plates),Cr_est);
        
    end
    [acc_new,Fpos_new,Fvel_new,Fpar_new]  = Fmod.compAccelPar(X_state,t);
    a_srp = acc_new;
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
    
    if estimated_coeff.Cr
        Fa_par = [Fa_par, Fpar_new];
        Fpar_dot = [Fpar_dot; zeros(Ncr,N_st)];
        Xpar_dot = [Xpar_dot; zeros(Ncr,1)];
        N_est = N_est+Ncr;
        Ba = [Ba;zeros(Ncr,1)];
    end
else
    a_srp = [0;0;0];
end



if parameters.drag
    %% density calculation
    omega_vec = [0;0;1]*omega_e;
    if estimated_coeff.Cd
        Cd_est = Cd_est+ X(N_nom,1);
        
        n_grid = repmat([0:Order_o],Order_b+1,1);
        m_grid = repmat([0:Order_b]',1,Order_o+1);
        Af_arg = cosd(n_grid*phi_in).*cosd(m_grid*theta_in(1));
        Bf_arg = cosd(n_grid*phi_in).*sind(m_grid*theta_in(1));
        Cf_arg = sind(n_grid*phi_in).*cosd(m_grid*theta_in(1));
        Df_arg = sind(n_grid*phi_in).*sind(m_grid*theta_in(1));
        % nominal entries corresponing to the indices
        Af_filter = Af_arg(Af_ind); Bf_filter = Bf_arg(Bf_ind); Cf_filter = Cf_arg(Cf_ind); Df_filter = Df_arg(Df_ind);
        % vectors from those entries
        argf_vec = [Af_filter(:); Bf_filter(:); Cf_filter(:); Df_filter(:)];
        argf_vec = argf_vec(2:end);
        
        Cd_est = Cd_est + sum(X_f.*argf_vec);
        [rho, delrho] = density_output(t, epoch, doy, time_var, X_state, parameters.sun_pos(:,t_ind), eps, parameters);
    elseif estimated_coeff.CdTrue
        parameters.t = t;
        parameters.time_jd = time_var;
        parameters.theta =  theta_v; % pi/180*parameters.theta(:,t_ind)'; %
        parameters.X_state = X_state;
        parameters.phi =  phi_v; % pi/180*parameters.phi(:,t_ind)';   %
        [rho,delrho,Cd_est, bff_coeff,S,m_r,frac,r_ads] = cd_truth_bff(t,parameters);
    elseif estimated_coeff.CdDiscrete
        [rho, delrho] = density_output(t, epoch, doy, time_var, X_state, parameters.sun_pos(:,t_ind), eps, parameters);
        Cd_est = parameters.Cd_discrete(t_ind);
    else
        [rho, delrho] = density_output(t, epoch, doy, time_var, X_state, parameters.sun_pos(:,t_ind), eps, parameters);
        bff_coeff = [zeros(Order_b+1,1),zeros(Order_b+1,1)];
    end
    if estimated_coeff.rho_DMC
        rho = rho + X(6+Ncr+1)+X(6+Ncr+3);
    end
    %     rho = rho*parameters.rho_step(t_ind);
    Fmod      = Drag_atm(parameters.mass, parameters.drag_area, Cd_est, omega_vec, rho, delrho, [0;0;0]);
    [acc_new,Fpos_new,Fvel_new,Fpar_new]  = Fmod.compAccelPar(X_state,t);
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
    
    if estimated_coeff.rho_DMC
        Fdmc_rho = [acc_new/rho, zeros(3,1), acc_new/rho];
        Fa_par = [Fa_par, Fdmc_rho];
        Flast_dmc_rho = [zeros(1,N_est), -tau_inv, 1,0, zeros(1,N_st-N_est-3); zeros(1,N_est), -omega^2, -2*zeta*omega,0, zeros(1,N_st-N_est-3);...
            zeros(1,N_est),0,0,0,zeros(1,N_st-N_est-3)];
        Fpar_dot = [Fpar_dot; Flast_dmc_rho];
        Ba = [Ba;c1;c2;c3];
        Xpar_dot = [Xpar_dot; -tau_inv*X(7+N_p)+X(8+N_p); -omega^2*X(7+N_p)-2*zeta*omega*X(8+N_p);0];
        N_est = N_est + 3;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if estimated_coeff.Cd
        if numel(argf_vec) > 0
            %         Fcd2 = Fcd.*argf_vec(vec_est((N_nom+1):end)-N_nom)';
            Fpar_new = Fpar_new.*argf_vec';
            Fa_par = [Fa_par, Fpar_new];
        end
        Fpar_dot = [Fpar_dot;zeros(1,N_nom-1),-tau_inv_cd, zeros(1,N_st-N_nom); zeros(N_st-N_nom,N_st)];
        Ba = [Ba;c_cd;zeros(N_st-N_nom,1)];
        Xpar_dot = [Xpar_dot;-tau_inv_cd*X(N_nom);zeros(N_st-N_nom,1)];
        N_est = N_est + (N_st-N_nom+1);
    end
    
    a_drag = acc_new;
   
%     Xpar_dot = [norm(a_drag)];
   
else
    a_drag = [0;0;0];
end

if parameters.relativity
    Fmod      = Relativity(parameters.c,parameters.mu_e);
    [acc_new,Fpos_new,Fvel_new,~]  = Fmod.compAccelPar(X_state,t);
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
end


if parameters.earth_rad
    if strcmp(parameters.flag_earthrad, 'Knocke')
        sod = mod(epoch+t, 86400);
        Fmod = bwalbedo(parameters.year - 2000,doy,sod,parameters.sun_pos(:,t_ind),parameters.area_vec,parameters.Area_plates,parameters.mass,...
            parameters.rho_spec,parameters.rho_diff,parameters.rho_spec_ir,parameters.rho_diff_ir,parameters.TEI',SBF2ECEF, parameters.Cr_erp);
        [acc_new,Fpos_new,Fvel_new,~]  = Fmod.compAccelPar(X_state,t);
        accel = accel + acc_new;
        Fpos = Fpos + Fpos_new;
        Fvel = Fvel + Fvel_new;
        a_earthrad = acc_new;
    else
        Fmod = ceres_erp(parameters.sun_pos(:,t_ind),n_hat,parameters.Area_plates,parameters.mass,parameters.rho_spec,parameters.rho_diff,...
            parameters.rho_spec_ir,parameters.rho_diff_ir,parameters.TEI', parameters.Psw, parameters.Plw, parameters.lon, parameters.lat, ...
            parameters.Cr_erp, parameters.Re);
        [acc_new,Fpos_new,Fvel_new,~]  = Fmod.compAccelPar(X_state,t);
        accel = accel + acc_new;
        Fpos = Fpos + Fpos_new;
        Fvel = Fvel + Fvel_new;
        a_earthrad = acc_new;
    end
else
    a_earthrad = [0;0;0];
end
X_state_dot = [X_state(4:6);accel;Xpar_dot];

%% Forming the dynamics vector
Xdot = [X_state_dot];                            % add dynamics of other states if there                          % add dynamics of other states if there

end


