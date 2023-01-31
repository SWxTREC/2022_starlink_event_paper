classdef ThreeBody
    properties
        Mu                                                                 % gravitational constant of third body
        X_third                                                            % position of third body in J2000 frame
    end
    
    methods
        
        function obj = ThreeBody(Mu, X_third)
            obj.Mu       = Mu;
            obj.X_third  = X_third;
        end
        
        function [a, F_ar, F_av, F_aPar] = compAccelPar(obj, X_state,~)
            X            = X_state(1:3,1);                                 % X_state = [vel;pos] w.r.t to primary
            X_c          = obj.X_third;
            mu           = obj.Mu;
            X_rel        = X_c - X;                                        % third body w.r.t spacecraft
            r_c          = norm(X_c);
            r_rel        = norm(X_rel);
            a            = mu*(X_rel/r_rel^3 - X_c/r_c^3);
            
            r_cap        = X_rel/r_rel;
            F_ar         = mu/(r_rel)^3*( 3*(r_cap*r_cap') - eye(3));  
            F_av         = zeros(3,3);
            F_aPar       = [];                                             % jacobian w.r.t to any parameter that needs to be estimated
        end
        
        function a = compAccel(obj,X_state,~)                              
            X            = X_state(1:3,1);                                 % X_state = [vel;pos] w.r.t to primary
            X_c          = obj.X_third;
            mu           = obj.Mu;
            X_rel        = X_c - X;                                        % third body w.r.t spacecraft
            r_c          = norm(X_c);
            r_rel        = norm(X_rel);
            a            = mu*(X_rel/r_rel^3 - X_c/r_c^3);
        end 
    end
    
    
end