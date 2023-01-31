classdef SRP_3constants
    properties
        SRP_am                            % Area to mass ratio of satellite
        AU                                % Astronomical Unit
        Cr                                % Radiation coefficient of satellite
        P_sun                             % Radiation pressure
        X_sun                             % Position of sun w.r.t Earth in J2000 frame
        Re                                % Radius of Earth
        Rsun                              % Radius of Sun
        earth_vel                         % velocity of earth
        A0                                % 3 constants
        A1
        A2
    end
    
    methods
        function obj =  SRP_3constants(SRP_am, AU, Cr, P_sun, X_sun, R, Rsun, earth_vel,A0,A1,A2)
            obj.SRP_am    = SRP_am;
            obj.AU        = AU;
            obj.Cr        = Cr;
            obj.P_sun     = P_sun;
            obj.X_sun     = X_sun;
            obj.Re        = R;
            obj.Rsun      = Rsun;
            obj.earth_vel = earth_vel;
            obj.A0        = A0;
            obj.A1        = A1;
            obj.A2        = A2;
        end
        function [a, F_ar, F_av,F_Cr] = compAccelPar(obj, X_state,~)
            X = X_state(1:3,1);               % X_state = [vel;pos] w.r.t to primary
            X_c = obj.X_sun;
            U_hat = X_c/norm(X_c);
            W_hat = cross(obj.earth_vel/norm(obj.earth_vel),U_hat);
            W_hat = W_hat/norm(W_hat);
            V_hat = cross(W_hat,U_hat);
            V_hat = V_hat/norm(V_hat);
            X_rel = X_c - X;                 % third body w.r.t spacecraft
            r_rel = norm(X_rel);
            Ps  = obj.P_sun;
            cr = obj.Cr;
            am = obj.SRP_am;
            au = obj.AU;
            fs = shadow_func(obj, X_state);
            
            a0 = obj.A0;
            a1 = obj.A1;
            a2 = obj.A2;
            
            K = Ps*cr*am*au^2*fs;
            a = K/r_rel^2*(a0*U_hat + a1*V_hat + a2*W_hat);
            
            F_ar = 2/r_rel^2*a*X_rel';
            
            F_av = zeros(3,3);
            
            K = Ps*cr*am*au^2*fs;
            F_Cr = K/r_rel^2*[U_hat,V_hat,W_hat];
            
        end
        
        function a = compAccel(obj,X_state,~)       % compute partial from positions
            X = X_state(1:3,1);               % X_state = [vel;pos] w.r.t to primary
            X_c = obj.X_sun;
            U_hat = X_c/norm(X_c);
            W_hat = [0;0;1];
            V_hat = cross(W_hat,U_hat);
            V_hat = V_hat/norm(V_hat);
            X_rel = X_c - X;                 % third body w.r.t spacecraft
            r_rel = norm(X_rel);
            Ps  = obj.P_sun;
            cr = obj.Cr;
            am = obj.SRP_am;
            au = obj.AU;
            fs = shadow_func(obj, X_state);
            
            a0 = obj.A0;
            a1 = obj.A1;
            a2 = obj.A2;
            
            K = Ps*cr*am*au^2*fs;
            a = K/r_rel^2*(a0*U_hat + a1*V_hat + a2*W_hat);
            
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