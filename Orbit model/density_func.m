%% Density function call
function [rho,M,rho_oxa,T] = density_func(rho_inputs)
%#codegen
flag_rho = rho_inputs.flag_rho;

switch flag_rho
    case 'Exp'
        rho0 = rho_inputs.rho0;
        H_scale = rho_inputs.H_scale;
        r0 = rho_inputs.r0;
        r = norm(rho_inputs.X_state(1:3));
        rho = rho0*exp((r0-r)/H_scale);
        
        atm_mass = rho_inputs.atm_mass;
        H_scale_all = rho_inputs.H_scale_all;
        rho0_all = rho_inputs.rho0_all;
        n_den = rho0_all.*exp((r0 - r)./H_scale_all);
        M = 1/(n_den/sum(n_den)*(1./atm_mass'));
        rho_oxa = n_den(:,2) + n_den(:,8);
        T  = rho_inputs.Talt;
    case 'MSIS00'
        day      = rho_inputs.doy;
        altitude = rho_inputs.altitude;
        days_year_prev = rho_inputs.days_year_prev;
        UTsec    = rho_inputs.UTsec;
        year     = rho_inputs.year;
        latitude = rho_inputs.latitude;
        F10_total= rho_inputs.F10_total;
        Ap_daily = rho_inputs.Ap_daily;
        Ap_total = rho_inputs.Ap_total;
        longitude = rho_inputs.longitude;
        doy_ind = day + days_year_prev;
        F10d = F10_total(doy_ind-1);                                   %% see nrlmsise_inputs_processing.m daily index
%         F10  = mean(F10_total((doy_ind-40):(doy_ind+40)));               %% average of 81 days centered on current day
        F10  = rho_inputs.F81c(doy_ind);               %% average of 81 days centered on current day
        UThour  = UTsec/3600;
        ap_daily = Ap_daily(doy_ind);                                  %% see nrlmsise_inputs_processing.m  ap daily index
        ap_ind = (doy_ind-1)*8+ceil(UThour/3);
        ap      = [ap_daily, Ap_total(ap_ind), Ap_total(ap_ind-1), Ap_total(ap_ind-2), Ap_total(ap_ind-3),...
            mean(Ap_total((ap_ind-11):(ap_ind-4))), mean(Ap_total((ap_ind-19):(ap_ind-12)))];
        ap(~any(ap,1)) = ap_daily;
        
        flags = ones(23,1);
        flags(9) = -1;
        [Temp,n_den] = atmosnrlmsise00(altitude, latitude, longitude, year, day, UTsec, F10, F10d, ap,flags,'Oxygen');
        rho     = n_den(:,6);
        he  = rho_inputs.atm_mass(1);
        oxa = rho_inputs.atm_mass(2);
        nitm= rho_inputs.atm_mass(3);
        oxm = rho_inputs.atm_mass(4);
        ar  = rho_inputs.atm_mass(5);
        h   = rho_inputs.atm_mass(6);
        nita= rho_inputs.atm_mass(7);
        amu = rho_inputs.amu;
        mass_sum = n_den(:,1)*he + n_den(:,2)*oxa + n_den(:,3)*nitm + n_den(:,4)*oxm + n_den(:,5)*ar + n_den(:,7)*h + n_den(:,8)*nita + n_den(:,9)*oxa;
        num = sum(n_den,2);
        M = mass_sum./num;
        rho_oxa     = n_den(:,2)*oxa*amu + n_den(:,9)*oxa*amu;
        T = Temp(2);
    case 'HP'
        day      = rho_inputs.doy;
        altitude = rho_inputs.altitude;
        days_year_prev = rho_inputs.days_year_prev;
        X_state  = rho_inputs.X_state;
        sun_pos  = rho_inputs.sun_pos;
        F10_total= rho_inputs.F10_total;
        F10_vec = [65,75,100,125,150,175,200,225,250,275];
        doy_ind = day + days_year_prev;
        F10d = F10_total(doy_ind);
        [~,f_ind] = min(abs(F10d-F10_vec));
        F10_input = F10_vec(f_ind);
        rho = Density_HP_new(altitude,sun_pos,X_state(1:3),F10_input);
        M = [];
        rho_oxa = [];
    case 'JB08'
        SUN = [0,0];
        SAT = [0,0,0];
        day      = rho_inputs.doy;
        altitude = rho_inputs.altitude;
        X_state  = rho_inputs.X_state;
        sun_pos  = rho_inputs.sun_pos;
        UTsec    = rho_inputs.UTsec;
        year     = rho_inputs.year;
        latitude = rho_inputs.latitude;
        flattening = rho_inputs.flattening;
        SOLdata    = rho_inputs.SOLdata;
        DTCdata    = rho_inputs.DTCdata;
        ind_sol    = rho_inputs.ind_sol;
        ind_mag    = rho_inputs.ind_mag;
        UThour  = UTsec/3600;
        hour = floor(UThour);
        minute = floor(UTsec/60) - hour*60;
        sec = UTsec - hour*3600 - minute*60;
        [month, day_mon, ~,~, ~] = days2mdh (year, day);
        MJD = Mjday(year,month,day_mon,hour,minute,sec);
        
        ra_Sun  = atan2(sun_pos(2), sun_pos(1));
        dec_Sun = atan2(sun_pos(3), sqrt(sun_pos(1)^2+sun_pos(2)^2));
        SUN(1)  = ra_Sun;
        SUN(2)  = dec_Sun;
        
        
        SAT(1) = atan2(X_state(2), X_state(1));                   % right-ascension
        SAT(2) = geocentricLatitude(latitude*pi/180,flattening,'radians');         % F = flattening    latitude
        SAT(3) = altitude/1000;                                    % altitude in km
        
        
        ind_sol = ind_sol+day-2;                  % one day lag
        
        SOL = SOLdata(:,ind_sol);
        F10 = SOL(4);
        F10B = SOL(5);
        S10 = SOL(6);
        S10B = SOL(7);
        
        % USE 2 DAY LAG FOR M10 FOR JB2008
        SOL = SOLdata(:,ind_sol-1);
        XM10 = SOL(8);
        XM10B = SOL(9);
        
        % USE 5 DAY LAG FOR Y10 FOR JB2008
        SOL = SOLdata(:,ind_sol-4);
        Y10 = SOL(10);
        Y10B = SOL(11);
        
        ind_mag = ind_mag+day-1;
        DTC = DTCdata(:,ind_mag);
        ii = floor(hour)+3;
        DSTDTC = DTC(ii);
        
        [Temp,rho,M,rho_oxa] = JB2008(MJD,SUN,SAT,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC);
        T = Temp(2);
    case 'HASDM'
        altitude = 1e-3*rho_inputs.altitude;
        longitude = rho_inputs.longitude;
        if longitude < 0
            longitude = longitude + 360;
        end
        latitude = rho_inputs.latitude;
        njd = (rho_inputs.time_jd - rho_inputs.jd_ref)/86400;
        rho = exp(rho_inputs.F_hasdm(altitude, njd,longitude, latitude));
        M = [];
        rho_oxa = [];
        T = [];
        %     case [{'CHAMP'}, {'HASDM_dis'}]
    case 'CHAMP'
        njd = (rho_inputs.time_jd - rho_inputs.jd_ref);
        rho_log = ppval(rho_inputs.F_champ,njd);
        rho = exp(rho_log);
        M = [];
        rho_oxa = [];
        T = [];
    case 'HASDM_dis'
        njd = (rho_inputs.time_jd - rho_inputs.jd_ref);
        rho_log = ppval(rho_inputs.F_champ,njd);
        rho = exp(rho_log);
        M = [];
        rho_oxa = [];
        T = []; 
    case 'WAMIPE'
        altitude = 1e-3*rho_inputs.altitude;
        longitude = rho_inputs.longitude;
        if longitude < 0
            longitude = longitude + 360;
        end
        latitude = rho_inputs.latitude;
        njd = (rho_inputs.time_jd - rho_inputs.jd_ref)/86400;
        rho = exp(rho_inputs.F_wamipe(longitude, latitude,altitude, njd));
        M = [];
        rho_oxa = [];
        T = []; 
end