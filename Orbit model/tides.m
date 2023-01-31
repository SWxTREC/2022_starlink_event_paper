%% Solid Earth and Ocean Tides
% Changing the spherical harmonic coefficients

function [C, S] = tides(C, S, time_et, Mjd_UTC, UT1_UTC, TEI, r_Moon, r_Sun, r_ref, gm, mu_m,mu_s, x_pole, y_pole, flag_tides, coeff0, coeff1, coeff2)

% convert moon and sun eci positions to ecef
r_Moon = TEI*r_Moon;
[lM, phiM, rM] = CalcPolarAngles(r_Moon);
r_Sun = TEI*r_Sun;
[lS, phiS, rS] = CalcPolarAngles(r_Sun);

Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;
T  = time_et/36525; % modified julian centuries of et

if (flag_tides.SolidEarthTides)
    % Effect of Solid Earth Tides (anelastic Earth)
    
    % Mean arguments of luni-solar motion
    %
    %   l   mean anomaly of the Moon
    %   l'  mean anomaly of the Sun
    %   F   mean argument of latitude
    %   D   mean longitude elongation of the Moon from the Sun
    %   Om  mean longitude of the ascending node of the Moon
    l  = (134.96340251 + 1717915923.2178*T + 31.879200*T^2 +0.051635*T^3 - 0.0002447*T^4)*pi/180;
    lp = (357.52910918 + 129596581.0481*T - 0.553200*T^2+0.000136*T^3 - 0.00001149*T^4)*pi/180;
    F  = (93.27209062 + 1739527262.847800*T- 12.7512*T^2 - 0.001037*T^3 + 0.00000417*T^4)*pi/180;
    D  = (297.85019547 + 1602961601.209*T - 6.3706*T^2 + 0.00659300*T^3 -0.00003169*T^4)*pi/180;
    Om = (125.04455501 -  6962890.543100*T + 7.472200*T^2 +0.007702*T^3 - 0.00005939*T^4)*pi/180;
    
    
    
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
    
    % STEP2 CORRECTIONS
    theta_f = -coeff1(:,1:5)*[l lp F D Om]';
    dC20 = 1e-12*(coeff1(:,6)'*cos(theta_f)-coeff1(:,7)'*sin(theta_f));

    dCnm20 = dCnm20 + dC20;
    
    theta_g = gmst(Mjd_UT1);

    
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
    
    % Treatment of the Permanent Tide (anelastic Earth)
    dC20 = 4.4228e-8*(-0.31460)*0.30190;
    dCnm20 = dCnm20 - dC20;
    
    % Effect of Solid Earth Pole Tide (anelastic Earth)
    dC21 = -1.348e-9*(x_pole+0.0112*y_pole);
    dS21 = 1.348e-9*(y_pole-0.0112*x_pole);
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
    [lgM, ~] = Legendre(6,6,phiM);
    [lgS, ~] = Legendre(6,6,phiS);
    
    dCnm20 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,1)...
           + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
           + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm21 = -0.3075/5*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,2)*sin(lM)...
           + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,2)*sin(lS) );
    dCnm22 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
    dSnm22 = -0.3075/5*( (mu_m/gm)*((r_ref/rM)^3)*lgM(3,3)*sin(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^3)*lgS(3,3)*sin(2*lS) );
    dCnm30 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,1)...
           + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,1) );
    dCnm31 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
           + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
    dSnm31 = -0.195/7*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
           + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
    dCnm32 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
    dSnm32 = -0.195/7*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
    dCnm33 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
           + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
    dSnm33 = -0.195/7*( (mu_m/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
           + (mu_s/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
    dCnm40 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,1)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,1) );
    dCnm41 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,2)*cos(lM)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,2)*cos(lS) );
    dSnm41 = -0.132/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,2)*sin(lM)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,2)*sin(lS) );
    dCnm42 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,3)*cos(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,3)*cos(2*lS) );
    dSnm42 = -0.132/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,3)*sin(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,3)*sin(2*lS) );
    dCnm43 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,4)*cos(3*lM)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,4)*cos(3*lS) );
    dSnm43 = -0.132/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,4)*sin(3*lM)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,4)*sin(3*lS) );
    dCnm44 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,5)*cos(4*lM)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,5)*cos(4*lS) );
    dSnm44 = -0.132/9*( (mu_m/gm)*((r_ref/rM)^5)*lgM(5,5)*sin(4*lM)...
           + (mu_s/gm)*((r_ref/rS)^5)*lgS(5,5)*sin(4*lS) );
    dCnm50 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,1)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,1) );
    dCnm51 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,2)*cos(lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,2)*cos(lS) );
    dSnm51 = -0.1032/9*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,2)*sin(lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,2)*sin(lS) );
    dCnm52 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,3)*cos(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,3)*cos(2*lS) );
    dSnm52 = -0.1032/9*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,3)*sin(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,3)*sin(2*lS) );
    dCnm53 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,4)*cos(3*lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,4)*cos(3*lS) );
    dSnm53 = -0.1032/9*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,4)*sin(3*lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,4)*sin(3*lS) );
    dCnm54 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,5)*cos(4*lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,5)*cos(4*lS) );
    dSnm54 = -0.1032/9*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,5)*sin(4*lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,5)*sin(4*lS) );
    dCnm55 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,6)*cos(5*lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,6)*cos(5*lS) );
    dSnm55 = -0.1032/9*( (mu_m/gm)*((r_ref/rM)^6)*lgM(6,6)*sin(5*lM)...
           + (mu_s/gm)*((r_ref/rS)^6)*lgS(6,6)*sin(5*lS) );
    dCnm60 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,1)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,1) );
    dCnm61 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,2)*cos(lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,2)*cos(lS) );
    dSnm61 = -0.0892/9*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,2)*sin(lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,2)*sin(lS) );
    dCnm62 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,3)*cos(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,3)*cos(2*lS) );
    dSnm62 = -0.0892/9*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,3)*sin(2*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,3)*sin(2*lS) );
    dCnm63 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,4)*cos(3*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,4)*cos(3*lS) );
    dSnm63 = -0.0892/9*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,4)*sin(3*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,4)*sin(3*lS) );
    dCnm64 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,5)*cos(4*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,5)*cos(4*lS) );
    dSnm64 = -0.0892/9*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,5)*sin(4*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,5)*sin(4*lS) );
    dCnm65 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,6)*cos(5*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,6)*cos(5*lS) );
    dSnm65 = -0.0892/9*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,6)*sin(5*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,6)*sin(5*lS) );
    dCnm66 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,7)*cos(6*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,7)*cos(6*lS) );
    dSnm66 = -0.0892/9*( (mu_m/gm)*((r_ref/rM)^7)*lgM(7,7)*sin(6*lM)...
           + (mu_s/gm)*((r_ref/rS)^7)*lgS(7,7)*sin(6*lS) );
    
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
    C(5,4) = C(5,4) + dCnm43;
    C(5,5) = C(5,5) + dCnm44;
    S(5,2) = S(5,2) + dSnm41;
    S(5,3) = S(5,3) + dSnm42;
    S(5,4) = S(5,4) + dSnm43;
    S(5,5) = S(5,5) + dSnm44;
    
    C(6,1) = C(6,1) + dCnm50;
    C(6,2) = C(6,2) + dCnm51;
    C(6,3) = C(6,3) + dCnm52;
    C(6,4) = C(6,4) + dCnm53;
    C(6,5) = C(6,5) + dCnm54;
    C(6,6) = C(6,6) + dCnm55;
    S(6,2) = S(6,2) + dSnm51;
    S(6,3) = S(6,3) + dSnm52;
    S(6,4) = S(6,4) + dSnm53;
    S(6,5) = S(6,5) + dSnm54;
    S(6,6) = S(6,6) + dSnm55;
    
    C(7,1) = C(7,1) + dCnm60;
    C(7,2) = C(7,2) + dCnm61;
    C(7,3) = C(7,3) + dCnm62;
    C(7,4) = C(7,4) + dCnm63;
    C(7,5) = C(7,5) + dCnm64;
    C(7,6) = C(7,6) + dCnm65;
    C(7,7) = C(7,7) + dCnm66;
    S(7,2) = S(7,2) + dSnm61;
    S(7,3) = S(7,3) + dSnm62;
    S(7,4) = S(7,4) + dSnm63;
    S(7,5) = S(7,5) + dSnm64;
    S(7,6) = S(7,6) + dSnm65;
    S(7,7) = S(7,7) + dSnm66;    
end






