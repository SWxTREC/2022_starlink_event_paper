classdef Lift_atm
    properties
        Mass                                                               % Mass of satellite
        Area                                                               % Area (cross-sectional flow projected) of satellite
        Cl                                                                 % Drag coefficient of satellite
        omega_vec                                                          % angular velocity of Earth
        rho                                                                % density
        delrho                                                             % density partial
        wind                                                               % winds
        An                                                                 % Area vector in eci frame
        N_plates
    end
    
    methods
        function obj = Lift_atm(Mass, Area, Cl, omega_vec, rho, delrho, wind,An,N_plates)    % constructor
            obj.Mass       = Mass;
            obj.Area       = Area;
            obj.Cl         = Cl;
            obj.omega_vec  = omega_vec;
            obj.rho        = rho;
            obj.delrho     = delrho;
            obj.wind       = wind;
            obj.An         = An;
            obj.N_plates   = N_plates;
        end
        
        function [a, F_ar, F_av, F_aCl] = compAccelPar(obj, X_state,~)
            V              = X_state(4:6);
            R              = X_state(1:3);
            r              = norm(R);
            OmegaE         = obj.omega_vec;
            Vr             = V -cross(OmegaE,R) - obj.wind;                % relative velocity to atmosphere
            v_mag          = norm(Vr);
            Vr_hat         = Vr/v_mag;
            vderterm       = (v_mag^2*eye(3) - Vr'*Vr)/v_mag^3;
            a              = zeros(3,1);
            dUNdV          = zeros(3,3);
            dUNdR          = zeros(3,3);
             % Computation of partials
            w_e            = OmegaE;
            om_mat         = [0,-w_e(3),w_e(2);w_e(3),0,-w_e(1);-w_e(2),w_e(1),0];
            for ii = 1:obj.N_plates
                Ulhat          = cross(cross(Vr_hat, obj.An(:,ii)),Vr_hat);
                Ulnorm         = vecnorm(Ulhat);
                Ulhat          = Ulhat/Ulnorm;
                a              = a - obj.Cl(ii)*Ulhat;
                nmat           = [0,-obj.An(3,ii),obj.An(2,ii);obj.An(3,ii),0,-obj.An(1,ii);-obj.An(2,ii),obj.An(1,ii),0];
                vmat           = [0,-Vr_hat(3),Vr_hat(2);Vr_hat(3),0,-Vr_hat(1);-Vr_hat(2),Vr_hat(1),0];
                vncr           = cross(Vr_hat,obj.An(:,ii));
                vnmat          = [0,-vncr(3),vncr(2);vncr(3),0,-vncr(1);-vncr(2),vncr(1),0];
                UN             = cross(cross(Vr_hat, obj.An(:,ii)), Vr_hat);
                UN_norm        = norm(UN);
                dUNdV_in       = vderterm*nmat*vmat + vnmat*vderterm;
                dUNdV          = dUNdV + v_mag^2*(UN_norm^2*dUNdV_in - UN*UN'*dUNdV_in)/UN_norm^3  + UN/UN_norm*Vr'; 
                
                dVhatdR        =  (-v_mag^2*om_mat + Vr*Vr'*om_mat)/v_mag^3;
                dUNdR_in       = dVhatdR*nmat*vmat + vnmat*dVhatdR;
                dUNdR          = dUNdR + v_mag^2*(UN_norm^2*dUNdR_in - UN*UN'*dUNdR_in)/UN_norm^3 - 2*UN/UN_norm*Vr'*om_mat; 
            end
            a = a*0.5*obj.rho*obj.Area/obj.Mass*norm(Vr)^2;
            F_av = -0.5*obj.rho*obj.Area/obj.Mass*dUNdV;
            F_ar = -0.5*obj.rho*obj.Area/obj.Mass*dUNdR;
            F_aCl          = [];
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
