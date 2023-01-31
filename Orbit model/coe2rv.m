%% coe2rv
function [r_eci, v_eci] = coe2rv(a,e,i,raan,w,theta,mu,string, special)
if strcmp(string, 'CE')                % circular equatorial
    w = 0;
    raan = 0;
    theta = special;            % true longitude
end
if strcmp(string, 'CI')              % circular inclined
    w = 0;
    theta = special;            % argument of latitude
end
if strcmp(string, 'EE')                % elliptical equatorial
    raan = 0;
    w =special;                 % longitude of periapsis
end

p = a*(1-e^2);
r_pqw = [p*cosd(theta)/(1+e*cosd(theta)); p*sind(theta)/(1+e*cosd(theta)); 0];
v_pqw = [-sqrt(mu/p)*sind(theta); sqrt(mu/p)*(e+cosd(theta)); 0 ];

R_pqw2eci = [cosd(raan)*cosd(w) - sind(raan)*sind(w)*cosd(i), -cosd(raan)*sind(w) - sind(raan)*cosd(w)*cosd(i), sind(raan)*sind(i);
             sind(raan)*cosd(w) + cosd(raan)*sind(w)*cosd(i), -sind(raan)*sind(w) + cosd(raan)*cosd(w)*cosd(i),-cosd(raan)*sind(i);
             sind(w)*sind(i),                                 cosd(w)*sind(i),                                  cosd(i)];
r_eci = R_pqw2eci*r_pqw;
v_eci = R_pqw2eci*v_pqw;

% h = cross(r_eci,v_eci);
% z_hat = [0;0;1];
% n_hat = cross(z_hat,h);
% n_hat = n_hat/norm(n_hat);
% e_vec = cross(v_eci,h)/mu - r_eci/norm(r_eci);
% e_mag = norm(e_vec);
% theta2 = acos(dot(-n_hat, e_vec/e_mag));
