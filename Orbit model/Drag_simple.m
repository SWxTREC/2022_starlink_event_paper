classdef Drag_simple
    properties
        Mass                              % Mass of satellite
        Area                              % Area (cross-sectional flow projected) of satellite
        Cd                                % Drag coefficient of satellite
        rho                               % density
        F                                 % wind factor 
    end
    
    methods
        
        function a = compAccel(obj,X_state, ~)       % compute partial from positions
            V = X_state(4:6);
            a = -0.5*obj.rho*obj.Cd*obj.Area*obj.F/obj.Mass*norm(V)*V;
            
        end
        
        
    end
    
    
end