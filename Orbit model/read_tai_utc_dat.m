function [dat] = read_tai_utc_dat();
% Get the most up-to-date leap second in the file 'tai-utc.dat': 
%   dat (= tai - utc)
fid = fopen('tai-utc.dat','r');
str = fgets(fid);
while ~feof(fid),str = fgets(fid);end
dat = str2num(str(37:42));

return

