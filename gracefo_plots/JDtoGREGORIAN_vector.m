function [Year,varargout] = JDtoGREGORIAN_vector(JD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Year,Mon,Day] = JDtoGREGORIAN(JD)
% [Year,Mon,Day,h,m,s] = JDtoGREGORIAN(JD)
% [Year,Doy] = JDtoGREGORIAN(JD)
% [Year,Doy,h,m,s] = JDtoGREGORIAN(JD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges      = [1,32,60,91,121,152,182,213,244,274,305,335,366];
edges_Leap = [1,32,61,92,122,153,183,214,245,275,306,336,367];

T1900 = ( JD - 2415019.5 )/365.25;
Year = 1900 + floor( T1900 );
LeapYrs = floor( (Year - 1900 - 1)*0.25 );
Days = ( JD - 2415019.5 ) - ( (Year - 1900)*365.0 + LeapYrs );

% Find Days from previous Year
i = find(Days < 1.0);
Year(i) = Year(i) - 1;
LeapYrs(i) = floor( (Year(i) - 1900 - 1)*0.25 );
Days(i) = ( JD(i) - 2415019.5 ) - ( (Year(i) - 1900)*365.0 + LeapYrs(i) );

% Compute the Day of Year
DayofYr = floor( Days );

% Find Month and Day
if any( nargout == [3,6] )
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

	% Assign Month and Day to varargout
	varargout{1} = Mon;
	varargout{2} = Day;
else
	varargout{1} = DayofYr;
end

if nargout > 3,
	% Find Hours, Minutes, Seconds
	Temp = ( Days - DayofYr )*24;
	% Hours
	varargout{nargout-3} = floor( Temp );
	% Minutes
	varargout{nargout-2} = floor(mod(Temp*60,60));
	% Seconds
	varargout{nargout-1} = mod(Temp*3600,60);
end

return
