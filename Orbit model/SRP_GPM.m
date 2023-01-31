classdef SRP_GPM
    properties
        Mass                              % Area to mass ratio of satellite
        AU                                % Astronomical Unit
        rho_spec                          % specular reflectivity
        rho_diff                          % diffuse reflectivity
        P_sun                             % Radiation pressure
        X_sun                             % Position of sun w.r.t Earth in J2000 frame
        Re                                % Radius of Earth
        Rsun                              % Radius of Sun
        n_hat                             % area normals of each face
        area                              % area of each face (actual area)
        N_plates
        Cr
    end
    
    methods
          function obj =  SRP_GPM(mass, AU, rho_spec, rho_diff, P_sun, X_sun, R, Rsun, n_hat,area,N_plates,Cr)
            obj.Mass      = mass;
            obj.AU        = AU;
            obj.rho_spec  = rho_spec;
            obj.rho_diff  = rho_diff;
            
            obj.P_sun     = P_sun;
            obj.X_sun     = X_sun;
            obj.Re        = R;
            obj.Rsun      = Rsun;
            obj.n_hat     = n_hat;
            obj.area      = area;
            obj.N_plates  = N_plates;
            obj.Cr        = Cr;
        end       
        
        function [a, F_ar, F_av, Fa_Cr] = compAccelPar(obj,X_state,~)       % compute partial from positions
%             global sp_angle
            X = X_state(1:3,1);               % X_state = [vel;pos] w.r.t to primary
            X_c = obj.X_sun;
            X_rel = X_c - X;                 % third body w.r.t spacecraft
            r_rel = norm(X_rel);
            Ps  = obj.P_sun;
            mass = obj.Mass;
            au = obj.AU;
            fs = shadow_func(obj, X_state);
            sun_vec = X_rel/r_rel;
            sun_vec_mat = repmat(sun_vec,1,obj.N_plates);
            cos_theta = dot(sun_vec_mat, obj.n_hat);
%             sp_angle = acosd(cos_theta);
            cos_theta(cos_theta<0) = 0;
            e_coeff = sum(obj.area.*cos_theta.*(1-obj.rho_spec));
            n_coeff_mat = repmat(obj.area.*obj.rho_spec.*cos_theta.^2 + obj.area.*obj.rho_diff.*cos_theta/3,3,1);
            n_comp = sum(n_coeff_mat.*obj.n_hat,2);
            a = -Ps*fs*au^2*obj.Cr/(mass*r_rel^2)*(e_coeff*sun_vec + 2*n_comp);
            
            F_ar = -Ps*fs*obj.Cr*au^2/r_rel^3*( 3*(sun_vec*sun_vec') - eye(3))*e_coeff/mass;
            F_av = zeros(3,3);
            Fa_Cr = -Ps*fs*au^2/(mass*r_rel^2)*(e_coeff*sun_vec + 2*n_comp);
        end
        
        function a = compAccel(obj,X_state,~)       % compute partial from positions
%             global sp_angle
            X = X_state(1:3,1);               % X_state = [vel;pos] w.r.t to primary
            X_c = obj.X_sun;
            X_rel = X_c - X;                 % third body w.r.t spacecraft
            r_rel = norm(X_rel);
            Ps  = obj.P_sun;
            mass = obj.Mass;
            au = obj.AU;
            fs = shadow_func(obj, X_state);
            sun_vec = X_rel/r_rel;
            sun_vec_mat = repmat(sun_vec,1,obj.N_plates);
            cos_theta = dot(sun_vec_mat, obj.n_hat);
%             sp_angle = acosd(cos_theta);
            cos_theta(cos_theta<0) = 0;
            e_coeff = sum(obj.area.*cos_theta.*(1-obj.rho_spec));
            n_coeff_mat = repmat(obj.area.*obj.rho_spec.*cos_theta.^2 + obj.area.*obj.rho_diff.*cos_theta/3,3,1);
            n_comp = sum(n_coeff_mat.*obj.n_hat,2);
            a = -Ps*fs*au^2*obj.Cr/(mass*r_rel^2)*(e_coeff*sun_vec + 2*n_comp);
            
        end
        
        function fs = shadow_func(obj, X_state, ~)
            X = X_state(1:3,1);               % X_state = [vel;pos] w.r.t to primary
            X_c = obj.X_sun;
            X_rel = X - X_c;                  %  spacecraft w.r.t third body
            a = asin(obj.Rsun/norm(X_rel));
            b = asin(obj.Re/norm(X));
            c = acos(X'*(X_rel)/(norm(X)*norm(X_rel)));
            if c <= abs(a-b)
                fs = 0;
            elseif c >= a+b
                fs = 1;
            elseif c < a+b && c > abs(a-b)
                x = (c^2+a^2-b^2)/(2*c);
                y = sqrt(a^2-x^2);
                A = a^2*acos(x/a) + b^2*acos((c-x)/b) - c*y;
                fs = 1-A/(pi*a^2);
            else 
                fs = 0;
            end
        end
        
    end
    
    
end