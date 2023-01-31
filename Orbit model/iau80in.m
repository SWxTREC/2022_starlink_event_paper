% ----------------------------------------------------------------------------
%
%                           function iau80in
%
%  this function initializes the nutation matricies needed for reduction
%    calculations. the routine needs the filename of the files as input.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    iar80       - integers for fk5 1980
%    rar80       - reals for fk5 1980             rad
%
%  outputs       :
%    none
%
%  locals        :
%    convrt      - conversion factor to degrees
%    i,j         - index
%
%  coupling      :
%    none        -
%
%  references    :
%
% [iar80,rar80] = iau80in();
% ----------------------------------------------------------------------------- }

function [iar80,rar80] = iau80in(nut80);

% ------------------------  implementation   -------------------
% 0.0001" to rad
convrt= 0.0001 * pi / (180*3600.0);



iar80 = nut80(:,1:5);
rar80 = nut80(:,6:9);


rar80(1:106,1:4) = rar80(1:106,1:4)*convrt;
