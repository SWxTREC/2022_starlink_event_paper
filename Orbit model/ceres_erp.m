classdef ceres_erp
    % Written by Vishal Ray 03/2022
    % Ref: Doornbos et al. air density models derived from multiple
    % satellite observations
    properties
        sunvector_j2000        % Sun coordinates in J2000 frame (m)
        nrm_eci                % panel normal vectors in eci frame
        plt_area               % panel areas (m^2)
        scmass                 % satellite mass (kg)
        refl_vis_spec          % visible specular reflectivity coefficient
        refl_vis_diff          % visible diffusive reflectivity coefficient
        refl_ir_spec           % infrared specular reflectivity coefficient
        refl_ir_diff           % infrared diffusive reflectivity coefficient
        ECEFtoJ2000            % rotation matrix from ECEF to J2000 frame
        Psw                    % outgoing shortwave radiation pressure (flux/c)
        Plw                    % outgoing longwave radiation pressure (flux/c)
        lam                    % longitude of grid (degrees)
        phi                    % latitude of grid (degrees)
        Cr_erp                 % scaling factor
        Re                     % radius of Earth
    end
    
    methods
        function obj =  ceres_erp(sunvector_j2000,nrm_eci,plt_area,scmass,refl_vis_spec,refl_vis_diff,refl_ir_spec,refl_ir_diff,ECEFtoJ2000, ...
                Psw, Plw, lam, phi, Cr_erp, Re)
            obj.sunvector_j2000     = sunvector_j2000;
            obj.nrm_eci    = nrm_eci;
            obj.plt_area   = plt_area;
            obj.scmass     = scmass;
            obj.refl_vis_spec  = refl_vis_spec;
            obj.refl_vis_diff  = refl_vis_diff;
            obj.refl_ir_spec   = refl_ir_spec;
            obj.refl_ir_diff   = refl_ir_diff;
            obj.ECEFtoJ2000    = ECEFtoJ2000;
            obj.Psw = Psw;
            obj.Plw = Plw;
            obj.lam = lam;
            obj.phi = phi;
            obj.Cr_erp = Cr_erp;
            obj.Re = Re;
        end
         function [a, F_ar, F_av, F_aPar] = compAccelPar(obj, X_state,~)
            
            
            a  = compAccel(obj, X_state);
            delp = 100; 
            F_ar(:,1) = (compAccel(obj, X_state+ [delp;0;0;0;0;0]) - compAccel(obj, X_state- [delp;0;0;0;0;0]))/(2*delp);
            F_ar(:,2) = (compAccel(obj, X_state+ [0;delp;0;0;0;0]) - compAccel(obj, X_state- [0;delp;0;0;0;0]))/(2*delp);
            F_ar(:,3) = (compAccel(obj, X_state+ [0;0;delp;0;0;0]) - compAccel(obj, X_state- [0;0;delp;0;0;0]))/(2*delp);
            
            F_av = zeros(3,3);
            F_aPar = a/obj.Cr_erp;       
        
        function a = compAccel(obj, X_state,~)
            
            a_alb = zeros(3,1);
            a_ir = zeros(3,1);
            
%             for m = 1:M
%                 phik = -90 +(m-0.5)*obj.dlam;                               % longitude of grid cell
            M = numel(obj.phi);
            L = numel(obj.lam);
            CfA_alb = zeros(3,L,M);
            CfA_ir  = zeros(3,L,M);
            m = [1:M];
            l = [1:L];
                Ak = 4*pi*obj.Re^2/L*sin(pi/2/M)*sin((m-0.5)*pi/M);         % area of grid cell [1xM]
                Ak = repmat(Ak, L, 1);                                      % [LxM]
%                 for l = 1:L
%                     lamk = 180 - (l-0.5)*obj.dphi;                          % latitude of grid cell  
                    rk_ecef(1,:,:) = cosd(obj.lam)'*cosd(obj.phi);             % grid cell outward normal
                    rk_ecef(2,:,:) = sind(obj.lam)'*cosd(obj.phi);
                    rk_ecef(3,:,:) = repmat(sind(obj.phi), L,1);
                    rk = pagemtimes(obj.ECEFtoJ2000,rk_ecef);                              % [3xLxM]
                    rksat = X_state(1:3) - rk;                                   % vector from grid cell to satellite   
                    rksun = obj.sunvector_j2000 - rk;                       % from grid cell to sun
                    
                    rk_norm_mat = repmat(vecnorm(rk,2,1),3,1,1);                % [3xLxM]
                    uk = rk./rk_norm_mat;  
                    rksat_norm = vecnorm(rksat,2,1);
                    rksat_norm_mat = repmat(rksat_norm,3,1,1);
                    uksat = rksat./rksat_norm_mat;
                    rksun_norm_mat = repmat(vecnorm(rksun,2,1),3,1,1);
                    uksun = rksun./rksun_norm_mat;
                    sinEksun = squeeze(dot(uk,uksun));                                % elevation angle from grid cell [LxM]
                    sinEksat = squeeze(dot(uk,uksat));
                    
                    sinEksat(sinEksat<0) = 0;                             
                    sinEksun(sinEksun>0) = 1;
                    sinEksun(sinEksun<=0) = 0;
                    
                        Pkalb = obj.Psw.*Ak.*sinEksat.*sinEksun./squeeze(rksat_norm).^2/pi; % radiation pressure [LxM]
                        Pkir = obj.Plw.*Ak.*sinEksat./squeeze(rksat_norm).^2/pi;
                    
                        % radiation coefficient vector calculation
                        for kk = 1:numel(obj.plt_area)
                            plate_vec = repmat(obj.nrm_eci(:,kk), 1, L, M);         % [3xLxM]
                    cos_theta = dot(-uksat, plate_vec);                     % [1xLxM]
                    cos_theta(cos_theta<0) = 0;
                    e_coeff = obj.plt_area(kk)*cos_theta*(1-obj.refl_vis_spec(kk)); % [1xLxM]
                    e_coeff_mat = repmat(e_coeff,3,1,1);
                    n_coeff_mat = repmat(obj.plt_area(kk)*obj.refl_vis_spec(kk)*cos_theta.^2 + obj.plt_area(kk)*obj.refl_vis_diff(kk)*cos_theta/3, 3,1,1);
                    % [3xLxM]
                    n_comp = n_coeff_mat.*repmat(obj.nrm_eci(:,kk),1,L,M);
                    CfA_alb = CfA_alb + (-e_coeff_mat.*uksat + 2*n_comp);
                    
                    e_coeff = obj.plt_area(kk)*cos_theta*(1-obj.refl_ir_spec(kk)); % [1xLxM]
                    e_coeff_mat = repmat(e_coeff,3,1,1);
                    n_coeff_mat = repmat(obj.plt_area(kk)*obj.refl_ir_spec(kk)*cos_theta.^2 + obj.plt_area(kk)*obj.refl_ir_diff(kk)*cos_theta/3, 3,1,1);
                    % [3xLxM]
                    n_comp = n_coeff_mat.*repmat(obj.nrm_eci(:,kk),1,L,M);
                    CfA_ir = CfA_ir + (-e_coeff_mat.*uksat + 2*n_comp);   %[3xLxM]
                        end
                    
                    a_alb(1,1) =  - sum(sum(squeeze(CfA_alb(1,:,:)).*Pkalb))/obj.scmass;
                    a_alb(2,1) =  - sum(sum(squeeze(CfA_alb(2,:,:)).*Pkalb))/obj.scmass;
                    a_alb(3,1) =  - sum(sum(squeeze(CfA_alb(3,:,:)).*Pkalb))/obj.scmass;
                    
                    a_ir(1,1) =  - sum(sum(squeeze(CfA_ir(1,:,:)).*Pkir))/obj.scmass;
                    a_ir(2,1) =  - sum(sum(squeeze(CfA_ir(2,:,:)).*Pkir))/obj.scmass;
                    a_ir(3,1) =  - sum(sum(squeeze(CfA_ir(3,:,:)).*Pkir))/obj.scmass;
                    
            a = obj.Cr_erp*(a_alb + a_ir);
        end
        end
    end
end