%% Function to calculate drag coeffcients based on model 1
function [Cd,Cl] = sentman(Vi, A, Ar,Y,Ri,Alpha,Tw,S,l)
%% Equations (Doornbos)
G = 1./(2*S.*S);
P = exp(-Y.*Y.*S.*S)./S;
Q = 1+G;
Z = 1+erf(Y*S);
Vr = Vi.*sqrt(0.5*(1+Alpha*(4*Ri*Tw./Vi.^2 -1)));

Cd = (P/sqrt(pi) + Y.*Q.*Z + Y/2.*Vr./Vi.*(Y*sqrt(pi).*Z + P))*A/Ar;
Cl = (l.*G.*Z + l/2.*Vr./Vi.*(Y*sqrt(pi).*Z + P))*A/Ar;