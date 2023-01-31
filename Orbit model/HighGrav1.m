classdef HighGrav1
    properties
        deg                                             % degree of model
        ord                                             % order of model
        mu_e                                            % grav coeff in m3/s2  
        Re                                              % radius of earth in m 
        Cbar                                            % Normalized grav. coeff. in 2d matrix, see GravCoeffArrange.m in data
        Sbar
        TEI
    end
    methods
        function obj = HighGrav1(deg,ord, mu, R, Cbar, Sbar, Rot_I2B)
            obj.deg           = deg;
            obj.ord           = ord;
            obj.mu_e            = mu;
            obj.Re             = R;
            obj.Cbar          = Cbar;
            obj.Sbar          = Sbar;
            obj.TEI       = Rot_I2B;      
        end 
        function [a, F_ar, F_av, F_aPar] = compAccelPar(obj,X_state,~)
        coder.extrinsic('GravAccelPartialsExterior_mex')
            P_ecef = obj.TEI*X_state(1:3);
            [g_ecef,Jac] = GravAccelPartialsExterior_mex(obj.deg, obj.ord, obj.Re/1000, obj.mu_e*1e-9, P_ecef(1:3)'/1000, obj.Cbar, obj.Sbar);
            a = obj.TEI'*g_ecef*1e3;
            
            F_ar = obj.TEI'*Jac*obj.TEI;
            
            F_av = zeros(3,3); 
            F_aPar = [];
        end            
        
        
        %%%%%%%%%%%%%%% Acceleration function hasnt been modified to take
        %%%%%%%%%%%%%%% order as an input yet
        function a = compAccel(obj,X_state,~)            % X_state = [vel;pos]
            P_ecef = obj.TEI*X_state(1:3);
            g_ecef = GravAccelExterior_mex(obj.deg, obj.ord, obj.Re/1000, obj.mu_e*1e-9, P_ecef(1:3)'/1000, obj.Cbar, obj.Sbar);
            a = obj.TEI'*g_ecef*1e3;
        end
    end
    
    
end