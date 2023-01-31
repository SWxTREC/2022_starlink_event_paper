classdef SRP %#codegen
    properties
        SRP_am                            % Area to mass ratio of satellite
        AU                                % Astronomical Unit
        Cr                                % Radiation coefficient of satellite
        P_sun                             % Radiation pressure
        X_sun                             % Position of sun w.r.t primary body in J2000 frame
        R                                 % Radius of primary body
        Rsun                              % Radius of Sun
    end
    
    methods
        function obj =  SRP(SRP_am, AU, Cr, P_sun, X_sun, R, Rsun)
            obj.SRP_am    = SRP_am;
            obj.AU        = AU;
            obj.Cr        = Cr;
            obj.P_sun     = P_sun;
            obj.X_sun     = X_sun;
            obj.R         = R;
            obj.Rsun      = Rsun;
        end 
        function [a, F_ar, F_av, F_aCr] = compAccelPar(obj, X_state,~)
            X             = X_state(1:3,1);                                % X_state = [vel;pos] w.r.t to primary
            X_c           = obj.X_sun;               
            X_rel         = X_c - X;                                       % sun w.r.t spacecraft
            r_rel         = norm(X_rel);
            Ps            = obj.P_sun;
            cr            = obj.Cr;
            am            = obj.SRP_am;
            au            = obj.AU;
            fs            = shadow_func(obj, X_state);
            a             = -Ps*cr*fs*am*X_rel/r_rel^3*au^2; 
            
            % Computation of partials
            r_cap         = X_rel/r_rel;
            F_ar          = -Ps*fs*cr*am*au^2/r_rel^3*( 3*(r_cap*r_cap') - eye(3)); 
            F_av          = zeros(3,3);
            F_aCr         = -Ps*fs*am*X_rel/r_rel^3*au^2; 
        end
        
        function a = compAccel(obj,X_state,~)                       
            X             = X_state(1:3,1);               
            X_c           = obj.X_sun;               
            X_rel         = X_c - X;                
            r_rel         = norm(X_rel);
            Ps            = obj.P_sun;
            cr            = obj.Cr;
            am            = obj.SRP_am;
            au            = obj.AU;
            fs            = shadow_func(obj, X_state);
            a             = -Ps*cr*fs*am*X_rel/r_rel^3*au^2;  
        end
        
        
        function fs = shadow_func(obj, X_state, ~)
            X             = X_state(1:3,1);             
            X_c           = obj.X_sun;               
            X_rel         = X - X_c;               
            a             = asin(obj.Rsun/norm(X_rel));
            b             = asin(obj.R/norm(X));
            c             = acos(X'*(X_rel)/(norm(X)*norm(X_rel)));
            if c <= abs(a-b)
                fs        = 0;
            elseif c >= a+b
                fs        = 1;
            elseif c < a+b && c > abs(a-b)
                x         = (c^2+a^2-b^2)/(2*c);
                y         = sqrt(a^2-x^2);
                A         = a^2*acos(x/a) + b^2*acos((c-x)/b) - c*y;
                fs        = 1-A/(pi*a^2);
            else
                fs = 1;
            end
        end 
            
    end
    
    
end