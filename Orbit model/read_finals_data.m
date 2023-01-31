function [eop] = read_finals_data()
%% Formatting of finals.data file
%fprintf('      ');fprintf('         %d',[1:9,0:9]);fprintf('\n      ');fprintf('%d',repmat([1:9,0],1,19));dbtype finals.data 1
%         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9
%1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
%92 1 1 48622.00 I  0.182987 0.000672  0.168775 0.000345  I-0.1251659 0.0000207  1.8335 0.0201  I   -10.437     .507     -.917     .165   .182400   .167900  -.1253000    -9.900    -1.700
%3I2   ,F8.2    , , F9.6    ,        , F9.6    ,          , F10.7    ,         , F7.4  ,         , F9.3    ,        , F9.3    ,
%yymmdd,--fmjd--,X,----xp---,XXXXXXXX,----yp---,XXXXXXXXXX,---dut1---,XXXXXXXXX,--rlod-,XXXXXXXXX,---dpsi--,XXXXXXXX,---deps--,
%
% note: LOD drops off about 1-7 days before the date that finals.data was created
%       predictions of dpsi and deps drop off about 75 days after the date that finals.data was created
%       predictions of xp, yp, and dut1 drop off about 1year+1week after the date that finals.data was created

% Get number of lines in the file
[s,w] = system('wc finals1.data');
leng = min(sscanf(w,'%d'));

fid = fopen('finals1.data','r');

str = fgets(fid);

[yy,mm,dd,fmjd,xp,yp,dut1,rlod,dpsi,deps] = deal(zeros(1,leng));
i=0;
while  i < 11047
	i=i+1;
	yy(i)=str2num(str(1:2));
	mm(i)=str2num(str(3:4));
	dd(i)=str2num(str(5:6));
	fmjd(i)=str2num(str(8:15));
	xp(i)=str2num(str(19:27));
	yp(i)=str2num(str(38:46));
	dut1(i)=str2num(str(59:68));
	temp=str2num(str(80:86));
	if ~isempty(temp),rlod(i)=temp;end
	temp=str2num(str(98:106));
	if isempty(temp),break;end
	dpsi(i)=temp;
	deps(i)=str2num(str(117:125));

	str = fgets(fid);
end
i = 1:i;
[eop.yy,eop.mm,eop.dd,eop.fmjd,eop.xp,eop.yp,eop.dut1,eop.rlod,eop.dpsi,eop.deps] = ...
	deal( yy(i),mm(i),dd(i),fmjd(i),xp(i),yp(i),dut1(i),rlod(i),dpsi(i),deps(i) );

return

