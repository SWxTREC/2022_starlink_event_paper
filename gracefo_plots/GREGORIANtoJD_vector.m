function [JD] = GREGORIANtoJD(Year,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [JD] = GREGORIANtoJD(Year,Mon,Day)
% [JD] = GREGORIANtoJD(Year,Mon,Day,h,m,s)
% [JD] = GREGORIANtoJD(Year,DayofYr)
% [JD] = GREGORIANtoJD(Year,DayofYr,h,m,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any( nargin == [2,5] ),
	edges      = [1,32,60,91,121,152,182,213,244,274,305,335,366];
	edges_Leap = [1,32,61,92,122,153,183,214,245,275,306,336,367];

	% Assign input Variables
	DayofYr = floor( varargin{1} );

	% Find Month and Day
	[Mon,Day] = deal( zeros(size(DayofYr)) );

	% Find Non-Leap Years
	i = find(mod(Year,4)~=0);
	if ~isempty(i),
		[n,bin] = histc(DayofYr(i),edges);
		Mon(i) = bin;
		Day(i) = DayofYr(i) - (edges( bin )-1);
	end

	% Find Leap Years
	i = find(mod(Year,4)==0);
	if ~isempty(i),
		[n,bin] = histc(DayofYr(i),edges_Leap);
		Mon(i) = bin;
		Day(i) = DayofYr(i) - (edges_Leap( bin )-1);
	end
else,
	[Mon,Day] = deal( varargin{1},varargin{2} );
end

if any( nargin == [5,6] ),
	% Assign input Variables
	[h,m,s] = deal( varargin{end-2},varargin{end-1},varargin{end} );
else,
	[h,m,s] = deal( zeros(size(Year)) );
end

% Calculate the JD from Year, Mon, Day
JD = 367*Year - floor( 7*(Year+floor((Mon+9)/12))/4 ) + floor(275*Mon/9) + Day + 1721013.5 + ((s/60+m)/60+h)/24;

return
