function [data] = read_gracefo_densities(fn)

% test case
%fn = '2022/GF_OPER_DNS1ACC_2__20220101T000000_20220101T235959_0001.cdf';
a = cdfread(fn);
cinfo = cdfinfo(fn);

% assign data to variables
for i = 1:size(cinfo.Variables,1),
	if strcmp(cinfo.Variables{i,1},'time');
		% convert CDF epoch to understandable time system
		data.time = zeros(1,size(a,1));
		for j = 1:size(a,1),
			data.t(j) = todatenum(a{j,i});
		end % j
		data.sod = round(86400*(data.t-floor(data.t)));
		data.t = data.t + 1721058.5; % converts datenums into Julian Date; corresponds to GREGORIANtoJD_vector(0000,1,1)+14
		[data.year,data.doy] = JDtoGREGORIAN_vector(data.t);
	else,
		eval(sprintf('data.%s = cell2mat(a(:,i))'';',cinfo.Variables{i,1}));
	end
end
		
return
