classdef bwalbedo
    % Written by Eric Sutton
    % Modified by Vishal Ray
    properties
        yy                     % 2-digit year
        doy                    % Day of Year
        sec                    % Second of the day
        sunvector_j2000         % Sun coordinates in J2000 frame (km)
        nrm_sbf                % panel normal vectors in satellite-body fixed (sbf) frame
        plt_area               % panel areas (m^2)
        scmass                 % satellite mass (kg)
        refl_vis_spec          % visible specular reflectivity coefficient
        refl_vis_diff          % visible diffusive reflectivity coefficient
        refl_ir_spec           % infrared specular reflectivity coefficient
        refl_ir_diff           % infrared diffusive reflectivity coefficient
        ECEFtoJ2000            % rotation matrix from ECEF to J2000 frame
        SBFtoECEF              % rotation matrix from satellite-body fixed (SBF) to Earth-centered Earth-fixed (ECEF) frame
        Cr_erp                 % global scaling factor
    end
    properties (Constant)
        num_rings   = 6;  %number of rings, not counting center spot
        Re          = 6378137;      %meters
        AUtoM      = 149597870.0e3;  %m/AU
        SFat1AU     = 1353;         %W/m^2
        c           = 2.99792458E8; %m/s
        IRflux      = 237; %W/m^2
    end
    
    
    % Outputs:
    % albedo = albedo in EFECF frame
    % infrared = infrared in ECEF frame
    methods
        function obj = bwalbedo(yy,doy,sec,sunvector_j2000,nrm_sbf,plt_area,scmass,refl_vis_spec,refl_vis_diff,refl_ir_spec,refl_ir_diff,ECEFtoJ2000,SBFtoECEF,Cr_erp)
            obj.yy     = yy;
            obj.doy    = doy;
            obj.sec    = sec;
            obj.sunvector_j2000     = sunvector_j2000;
            obj.nrm_sbf     = nrm_sbf;
            obj.plt_area     = plt_area;
            obj.scmass     = scmass;
            obj.refl_vis_spec  = refl_vis_spec;
            obj.refl_vis_diff  = refl_vis_diff;
            obj.refl_ir_spec  = refl_ir_spec;
            obj.refl_ir_diff  = refl_ir_diff;
            obj.ECEFtoJ2000  = ECEFtoJ2000;
            obj.SBFtoECEF  = SBFtoECEF;
            obj.Cr_erp     = Cr_erp;
        end
        
        
        function [a, F_ar, F_av, F_aPar] = compAccelPar(obj, X_state,~)
            
            
            a  = compAccel(obj, X_state);
            delp = 100; 
            F_ar(:,1) = (compAccel(obj, X_state+ [delp;0;0;0;0;0]) - compAccel(obj, X_state- [delp;0;0;0;0;0]))/(2*delp);
            F_ar(:,2) = (compAccel(obj, X_state+ [0;delp;0;0;0;0]) - compAccel(obj, X_state- [0;delp;0;0;0;0]))/(2*delp);
            F_ar(:,3) = (compAccel(obj, X_state+ [0;0;delp;0;0;0]) - compAccel(obj, X_state- [0;0;delp;0;0;0]))/(2*delp);
            
            F_av = zeros(3,3);
            F_aPar = a/obj.Cr_erp;
            
            
            
            function accel = compAccel(obj, X_state)
            satpos_ecef = obj.ECEFtoJ2000'*X_state(1:3);
            % t-to is days since dec. 21 0.0h 1981 (doy=355)
%             t_to = yydoyTOj2000(obj.yy,obj.doy) - yydoyTOj2000(100-81,355) + obj.sec/86400;
            t_to = yydoyTOj2000(obj.yy,obj.doy) + yydoyTOj2000(100-82,10.25)  + obj.sec/86400;
            cosTIME = cos(2*pi*mod(t_to,365.25)/365.25);
            
            sundist = norm(obj.sunvector_j2000/obj.AUtoM);   %has to be in AU
            Sflux = obj.SFat1AU/sundist^2;
            
            
            
            % rotate sununit_j2000 into ecef
            sununit_j2000 = obj.sunvector_j2000/(sundist*obj.AUtoM);
            sununit_ecef = obj.ECEFtoJ2000'*sununit_j2000;
            % rotate nrmls into ecef
            ecef_nrm = obj.SBFtoECEF*obj.nrm_sbf;
            
            % in ECEF frame
            albedo = zeros(3,1);
            infrared = zeros(3,1);
            
            R = norm(satpos_ecef);
            Ttotal = acos(obj.Re/R);
            dphi = Ttotal/(obj.num_rings+0.5);
            
            normxy = norm(satpos_ecef(1:2));
            c2 = satpos_ecef(1)/normxy;
            s2 = satpos_ecef(2)/normxy;
            c1 = satpos_ecef(3)/R;
            s1 = normxy/R;
            rot23 = [c2*c1,-s2,c2*s1;s2*c1,c2,s2*s1;-s1,0,c1];
            
            for n = 0:obj.num_rings
                
                phin = pi/2-n*dphi;
                cpn  = cos(phin);
                spn  = sin(phin);
                if n > 0
                    nseg = round(2*sin(dphi/2)*cpn/(1-cos(dphi/2)));
                    dlam = 2*pi/nseg;
                    areaseg = 2*dlam*obj.Re^2*sin(dphi/2)*cpn;
                    Dn   = sqrt(R^2 - 2*R*obj.Re*spn + obj.Re^2);
                else
                    nseg = 1;
                    dlam = 2*pi;
                    areaseg = 2*dlam*obj.Re^2*sin(dphi/4)*cos(phin-dphi/4);
                    Dn      = R - obj.Re;
                end
                
                for m = 1:nseg
                    lamm    = (m-1)*dlam;     %in reference to satellite lat/lon
                    xyz     = obj.Re*[cpn*cos(lamm);cpn*sin(lamm);sin(phin)];   %in reference to satellite lat/lon
                    spot_ecef = rot23*xyz;
                    spotlat = asin(spot_ecef(3)/obj.Re);                          %in reference to ecef
                    sinspotlat = sin(spotlat);
                    spotlon = atan2(spot_ecef(2),spot_ecef(1));                 %in reference to ecef
                    u       = (spot_ecef - satpos_ecef)/Dn;
                    rotZEN = [cos(spotlon)*cos(spotlat);sin(spotlon)*cos(spotlat);sin(spotlat)];
                    % compute satellite zenith angle from ring element n,m
                    cosZsat = rotZEN'*(-u);
                    % compute solar zenith angle from ring element n,m
                    cosZsun = rotZEN'*sununit_ecef;
                    
                    for i = 1:length(obj.nrm_sbf)
                        
                        cosT = dot(u,ecef_nrm(:,i)); % dotting two unit vectors
                        
                        if cosT > 0 % satellite panel i is facing toward surface element j
                            
                            % compute IR acceleration from earth (regardless of earth orientation to sun)
                            emm = 0.68 - 0.07*cosTIME*sinspotlat - 0.18*0.5*(3*sinspotlat^2 - 1);
                            term1 = (1-obj.refl_ir_spec(i))*u + 2*((1/3)*obj.refl_ir_diff(i) + obj.refl_ir_spec(i)*cosT)*ecef_nrm(:,i);
                            infrared = infrared - (Sflux*(emm/4)*cosZsat*cosT*obj.plt_area(i)*areaseg/(pi*obj.scmass*obj.c*Dn^2))*term1;
                            
                            if cosZsun > 0 %sun is above plane of earth surface element j
                                
                                % apply default albedo model
                                A = 0.34 + 0.1*cosTIME*sinspotlat + 0.29*0.5*(3*sinspotlat^2 - 1);
                                term2 = (1-obj.refl_vis_spec(i))*u + 2*((1/3)*obj.refl_vis_diff(i) + obj.refl_vis_spec(i)*cosT)*ecef_nrm(:,i);
                                albedo = albedo - (A*Sflux*cosZsun*cosZsat*cosT*obj.plt_area(i)*areaseg/(pi*obj.scmass*obj.c*Dn^2))*term2;
                                
                            end % cosZsun>0
                            
                        end % cosT>0
                        
                    end % i
                    
                    
                end
                
            end
            accel = obj.ECEFtoJ2000*(albedo + infrared)*obj.Cr_erp;
            end 
            
            
            
            function [j2000] = yydoyTOj2000(yy,doy)
                
                leapnum = floor((3+yy)/4);
                nonleapnum = yy - leapnum;
                
                j2000 = doy + 366*(leapnum) + 365*(nonleapnum) - 1.5;
            end
        end
    end
end