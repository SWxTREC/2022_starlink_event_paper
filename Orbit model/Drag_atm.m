classdef Drag_atm
    properties
        Mass                                                               % Mass of satellite
        Area                                                               % Area (cross-sectional flow projected) of satellite
        Cd                                                                 % Drag coefficient of satellite
        omega_vec                                                          % angular velocity of Earth
        rho                                                                % density
        delrho                                                             % density partial
        wind                                                               % winds
    end
    
    methods
        function obj = Drag_atm(Mass, Area, Cd, omega_vec, rho, delrho, wind)    % constructor
            obj.Mass       = Mass;
            obj.Area       = Area;
            obj.Cd         = Cd;
            obj.omega_vec  = omega_vec;
            obj.rho        = rho;
            obj.delrho     = delrho;
            obj.wind       = wind;
        end 
        
        function [a, F_ar, F_av, F_aCd] = compAccelPar(obj, X_state,~)
            V              = X_state(4:6);
            R              = X_state(1:3);
            r              = norm(R);
            r_cap          = R/r;            
            OmegaE         = obj.omega_vec;
            Vr             = V -cross(OmegaE,R) - obj.wind;                % relative velocity to atmosphere
            v              = norm(Vr);            
            a              = -0.5*obj.rho*obj.Cd*obj.Area/obj.Mass*norm(Vr)*Vr;

            % Computation of partials
            B              = obj.Mass/(obj.Area*obj.Cd);                   % ballistic coefficient
            w_e            = OmegaE;
            om_mat         = [0,-w_e(3),w_e(2);w_e(3),0,-w_e(1);-w_e(2),w_e(1),0];   
            F_ar           = -1/2/B*obj.delrho*v*Vr*r_cap' + obj.rho/2/B*(Vr*Vr'/v + v*eye(3))*om_mat;
            F_av           = -obj.rho/2/B*(Vr*Vr'/v + v*eye(3));   
            F_aCd          = -0.5*obj.rho*obj.Area/obj.Mass*norm(Vr)*Vr;
        end
        
        function a = compAccel(obj, X_state,~)
            V              = X_state(4:6);
            R              = X_state(1:3);          
            OmegaE         = obj.omega_vec;
            Vr             = V -cross(OmegaE,R) - obj.wind;                % relative velocity to atmosphere           
            a              = -0.5*obj.rho*obj.Cd*obj.Area/obj.Mass*norm(Vr)*Vr;

        end        
    end
    
    
end