classdef Relativity
    properties
        c_light              % speed of light
        mu                   % gravitational parameter
    end
    
    methods
        function obj = Relativity(c_light,mu)
            obj.c_light   = c_light;
            obj.mu        = mu;
        end
        function [a, F_ar, F_av, F_aPar] = compAccelPar(obj, X_state,~)
            v             = X_state(4:6);
            r             = X_state(1:3);
            % Relative position vector of satellite w.r.t. point mass
            r_Sat         = norm(r);
            v_Sat         = norm(v);
            
            % Acceleration
            a             = -obj.mu/(obj.c_light^2*r_Sat^3)*((4*obj.mu/r_Sat-v_Sat^2)*r+4*dot(r,v)*v);
            
            % Jacobians
            F_ar          = (3*obj.mu/(obj.c_light^2*r_Sat^4)*((4*obj.mu/r_Sat - v_Sat^2)*r + 4*dot(r,v)*v ) ...
                            + 4*obj.mu^2/(obj.c_light^2*r_Sat^5)*r)*r'/r_Sat- obj.mu/(obj.c_light^2*r_Sat^3)*((4*obj.mu/r_Sat - v_Sat^2)*eye(3,3) + 4*(v*v'));
            
            F_av          = 2*obj.mu/(obj.c_light^2*r_Sat^3)*r*v' - 4*obj.mu/(obj.c_light^2*r_Sat^3)*(v*r' + (r'*v)*eye(3,3));
            F_aPar = [];   
        end

        function a = compAccel(obj, X_state,~)
            v             = X_state(4:6);
            r             = X_state(1:3);
            % Relative position vector of satellite w.r.t. point mass
            r_Sat         = norm(r);
            v_Sat         = norm(v);
            
            % Acceleration
            a             = -obj.mu/(obj.c_light^2*r_Sat^3)*((4*obj.mu/r_Sat-v_Sat^2)*r+4*dot(r,v)*v);
            
        end        
    end
    
    
end