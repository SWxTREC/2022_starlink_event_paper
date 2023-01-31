%% ASEN 5050 HW 3 Calculating COE from radius and velocity vectors
function coe = rv2coe_E(r,v,mu)

r_mag = norm(r);
v_mag = norm(v);
r_hat = r/r_mag;
%% Semi major axis
epsilon = v_mag^2/2 - mu/r_mag;
a = -mu/(2*epsilon);

%% Eccentricity
h = cross(r,v);
h_mag = norm(h);
p = h_mag^2/mu;
% e = sqrt(1+2*epsilon*h_mag^2/mu^2);
e = sqrt(abs(1-p/a));
%% Inclination
i = acosd(h(3)/h_mag);

%% RAAN
z_hat = [0;0;1];
n_hat = cross(z_hat, h/h_mag);
% angle = asind(norm(cross(z_hat,h/h_mag)));
n_hat = n_hat/norm(n_hat);
RAAN = acosd(n_hat(1));
if n_hat(2)<0
    RAAN = -RAAN;
end
%%
%% Argument of periapsis
e_vec = cross(v,h)/mu - r/r_mag;
e_vec = e_vec/norm(e_vec);
cos_omg = dot(n_hat,e_vec);
if cos_omg > 1
    cos_omg=1;
elseif cos_omg < -1
    cos_omg = -1;
end
omega = acosd(cos_omg);
if e_vec(3)<0
    omega = -omega;
end
%% True anomaly
cos_theta = 1/e*(a*(1-e^2)/r_mag - 1); %dot(r_hat/norm(r_hat),e_vec/norm(e_vec));
if cos_theta > 1
    cos_theta = 1;
elseif cos_theta < -1
    cos_theta = -1;
end
theta = acosd(cos_theta);
rv = dot(r,v);
if rv < 0
    theta = -theta;
end
%% Argument of latitude
cos_u = dot(r/r_mag, n_hat/norm(n_hat));
if cos_u > 1
    cos_u = 1;
elseif cos_u < -1
    cos_u = -1;
end
u = acosd(cos_u);
if r(3) < 0
    u = -u;
end

%% Eccentric anomaly
if e < 1  %  e> 1e-6 &&
    sin_E = sind(theta)*sqrt(1-e^2)/(1+e*cosd(theta));
    cos_E = (e+cosd(theta))/(1+e*cosd(theta));
    E = atan2d(sin_E,cos_E);
%     M = E*pi/180-e*sind(E);
    M = u*pi/180;
else
    E = 0;
end
coe = [a,e,i,RAAN,omega,theta,u,E,M];