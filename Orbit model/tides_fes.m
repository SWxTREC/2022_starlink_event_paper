%% Solid Earth and Ocean Tides
% Changing the spherical harmonic coefficients
%#codegen
function [C, S] = tides_fes(C, S, time_et, jdTT_ref, mjd_ut1, TEI, pos_Moon1, pos_Sun1, r_ref, gm, mu_m,mu_s, x_pole, y_pole, flag_tides, ...
    coeff0, coeff1, coeff2, doo_coeff, Cnmp_ot, Snmp_ot,Cnmm_ot, Snmm_ot, deg_do, ord_do, pole_mean)

% convert moon and sun eci positions to ecef
pos_Moon = [0;0;0];
pos_Moon = TEI*pos_Moon1;
[lM, phiM, rM] = CalcPolarAngles(pos_Moon);
pos_Sun = [0;0;0];
pos_Sun = TEI*pos_Sun1;
[lS, phiS, rS] = CalcPolarAngles(pos_Sun);

T  = time_et/(86400*36525); % modified julian centuries of et from J2000

    % Mean arguments of luni-solar motion
    %
    %   l   mean anomaly of the Moon
    %   l'  mean anomaly of the Sun
    %   F   mean argument of latitude
    %   D   mean longitude elongation of the Moon from the Sun
    %   Om  mean longitude of the ascending node of the Moon
    l  = (134.96340251 + (1717915923.2178*T + 31.879200*T^2 +0.051635*T^3 - 0.0002447*T^4)/3600)*pi/180;
    lp = (357.52910918 + (129596581.0481*T - 0.553200*T^2+0.000136*T^3 - 0.00001149*T^4)/3600)*pi/180;
    F  = (93.27209062 + (1739527262.847800*T- 12.7512*T^2 - 0.001037*T^3 + 0.00000417*T^4)/3600)*pi/180;
    D  = (297.85019547 + (1602961601.209*T - 6.3706*T^2 + 0.00659300*T^3 -0.00003169*T^4)/3600)*pi/180;
    Om = (125.04455501 -  (6962890.543100*T + 7.472200*T^2 +0.007702*T^3 - 0.00005939*T^4)/3600)*pi/180;
    theta_g = gmst_iers10(jdTT_ref, mjd_ut1);
if (flag_tides.SolidEarthTides)
    % Effect of Solid Earth Tides (anelastic Earth)
    
 
    
    % STEP1 CORRECTIONS
    [lgM, ~] = Legendre(4,4,phiM);
    [lgS, ~] = Legendre(4,4,phiS);
    dCnm20 = (0.30190/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,1)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = (0.29830/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) )...
        +(-0.00144/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*(sin(lM))...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*(sin(lS)) );
    dSnm21 = (0.00144/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*(cos(lM))...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*(cos(lS)) )...
        + (0.29830/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*sin(lM)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*sin(lS) );
    dCnm22 = (0.30102/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) )...
        +(-0.00130/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*(sin(2*lM))...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*(sin(2*lS)) );
    dSnm22 = (0.00130/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*(cos(2*lM))...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*(cos(2*lS)) )...
        + (0.30102/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*sin(2*lM)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*sin(2*lS) );
    dCnm30 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,1)...
        + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,1) );
    dCnm31 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
        + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
    dSnm31 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
        + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
    dCnm32 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
        + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
    dSnm32 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
        + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
    dCnm33 = (0.094/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
        + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
    dSnm33 = (0.094/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
        + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
    dCnm40 = (-0.00089/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,1)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm41 = (-0.00080/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm41 = (-0.00080/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*sin(lM)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*sin(lS) );
    dCnm42 = (-0.00057/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
    dSnm42 = (-0.00057/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*sin(2*lM)...
        + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*sin(2*lS) );

%%%%%%%%% Elastic Earth%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dCnm20 = (0.30190/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,1)...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,1) );
%     dCnm21 = (0.29830/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
%     dSnm21 = (0.29830/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*(sin(lM))...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*(sin(lS)) );
%     dCnm22 = (0.30102/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
%     dSnm22 = (0.30102/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*(sin(2*lM))...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*(sin(2*lS)) );
%     dCnm30 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,1)...
%            + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,1) );
%     dCnm31 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
%            + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
%     dSnm31 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
%            + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
%     dCnm32 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
%            + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
%     dSnm32 = (0.093/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
%            + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
%     dCnm33 = (0.094/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
%            + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
%     dSnm33 = (0.094/7)*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
%            + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
%     dCnm40 = (-0.00087/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,1)...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,1) );
%     dCnm41 = (-0.00079/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
%     dSnm41 = (-0.00079/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*sin(lM)...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*sin(lS) );
%     dCnm42 = (-0.00057/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
%     dSnm42 = (-0.00057/5)*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*sin(2*lM)...
%            + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*sin(2*lS) );  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % STEP2 CORRECTIONS
    theta_f = -coeff1(:,1:5)*[l lp F D Om]';
    dC20 = 1e-12*(coeff1(:,6)'*cos(theta_f)-coeff1(:,7)'*sin(theta_f));
    dCnm20 = dCnm20 + dC20;
    

    
    theta_f = (theta_g+pi)-coeff0(:,1:5)*[l lp F D Om]';
    dC21    = 1e-12*(coeff0(:,6)'*sin(theta_f)+coeff0(:,7)'*cos(theta_f));
    dS21    = 1e-12*(coeff0(:,6)'*cos(theta_f)-coeff0(:,7)'*sin(theta_f));

    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    

    
    theta_f = 2*(theta_g+pi)-coeff2(:,1:5)*[l lp F D Om]';
    dC22    = 1e-12*coeff2(:,6)'*cos(theta_f);
    dS22    = - 1e-12*coeff2(:,6)'*sin(theta_f);

    dCnm22 = dCnm22 + dC22;
    dSnm22 = dSnm22 + dS22;
%     
    % Treatment of the Permanent Tide (anelastic Earth)
    dC20 = 4.4228e-8*(-0.31460)*0.30190;
%     dC20 = 4.4228e-8*(-0.31460)*0.29525;
    dCnm20 = dCnm20 - dC20;

    % Effect of Solid Earth Pole Tide (anelastic Earth)
    x_pole = x_pole*180/pi*3600;
    y_pole = y_pole*180/pi*3600;
    m1 = x_pole - pole_mean(1);
    m2 = y_pole - pole_mean(2);
    dC21 = -1.333e-9*(m1-0.0115*m2);
    dS21 = -1.333e-9*(m2+0.0115*m1);
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    

    C(3,1) = C(3,1) + dCnm20;
    C(3,2) = C(3,2) + dCnm21;
    C(3,3) = C(3,3) + dCnm22;
    S(3,2) = S(3,2) + dSnm21;
    S(3,3) = S(3,3) + dSnm22;
    
    C(4,1) = C(4,1) + dCnm30;
    C(4,2) = C(4,2) + dCnm31;
    C(4,3) = C(4,3) + dCnm32;
    C(4,4) = C(4,4) + dCnm33;
    S(4,2) = S(4,2) + dSnm31;
    S(4,3) = S(4,3) + dSnm32;
    S(4,4) = S(4,4) + dSnm33;
    
    C(5,1) = C(5,1) + dCnm40;
    C(5,2) = C(5,2) + dCnm41;
    C(5,3) = C(5,3) + dCnm42;
    S(5,2) = S(5,2) + dSnm41;
    S(5,3) = S(5,3) + dSnm42;
end

if (flag_tides.OceanTides)
    % Ocean Tides

    s = F+Om;
    beta = [theta_g+pi-s, s, s-D, s-l, -Om, s-D-lp]';
    theta_f = doo_coeff*beta;
    dCnm = 1e-12*((Cnmp_ot+Cnmm_ot)'*cos(theta_f) + (Snmp_ot+Snmm_ot)'*sin(theta_f));
    dSnm = 1e-12*((Cnmm_ot-Cnmp_ot)'*sin(theta_f) + (Snmp_ot-Snmm_ot)'*cos(theta_f));
    
    for ii= 1:numel(deg_do)
    C(deg_do(ii)+1, ord_do(ii)+1) = C(deg_do(ii)+1, ord_do(ii)+1) + dCnm(ii);
    S(deg_do(ii)+1, ord_do(ii)+1) = S(deg_do(ii)+1, ord_do(ii)+1) + dSnm(ii);
    end
     
end






