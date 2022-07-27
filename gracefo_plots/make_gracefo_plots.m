set(0,'defaultaxesfontsize',14);

% constants for making data patches
dtime = 0.035;
dlat = 0.65/2;

% plot x-axis limits
xlimit = [30,40];

% Define time period of files to load
doy = 29:40;
yy = 22*ones(size(doy));
[yyyy,mm,dd] = JDtoGREGORIAN_vector(GREGORIANtoJD_vector(2000+yy,doy));

% Pre-allocate memory
leng = length(doy);
[longD,longLat,longLon,longTime,longHeight,longValid] = deal( zeros(1,8640*leng) );

% Load data from multiple days
tleng = 0;
for idoy = 1:length(doy)

        % See if '.cdf' File Exists
        fn = sprintf('20%02d/GF_OPER_DNS1ACC_2__20%02d%02d%02dT000000_20%02d%02d%02dT235959_0001.cdf',...
		yy(idoy),yy(idoy),mm(idoy),dd(idoy),yy(idoy),mm(idoy),dd(idoy));
	disp(fn);
	[s,w] = system( sprintf('ls %s',fn) );

        % load data if file exists
        if s ~= 0, continue; end
        data = read_gracefo_densities(fn);

	% Append current day to timeseries
        leng = length(data.sod);
        i = tleng+1:tleng+leng;
        longD(i)     = data.density;
        longLat(i)   = data.latitude;
        longLon(i)   = data.longitude;
	longHeight(i) = data.altitude;
        longTime(i)  = 365*(yy(idoy)-yy(1)) + data.doy + data.sod/86400;
	if isfield(data,'validity_flag'),
        	longValid(i)   = data.validity_flag;
	end

	% increment size of timeseries arrays
        tleng = tleng + leng;
end

% truncate timeseries to total data points read in
longD = longD(1:tleng);
longLat = longLat(1:tleng);
longLon = longLon(1:tleng);
longTime = longTime(1:tleng);
longHeight = longHeight(1:tleng);
longValid = longValid(1:tleng);
longSlt = mod(longTime*24+longLon/15,24);


%% figure 1: Ascending Density
figure;
% find ascending portions of the orbits
i = find(diff(longLat) > 0);
i = union(i,i+1);
% Check validity
i = intersect(i,find( longD < 1 & longD > 0 & longValid == 0 ));

% find approximate local time of equatorial crossing
j = find( longLat(i) == 0 );
eqslt = mean(longSlt(i));

px = [longTime(i)-dtime;longTime(i)-dtime;longTime(i)+dtime;longTime(i)+dtime];
py = [longLat(i)-dlat;longLat(i)+dlat;longLat(i)+dlat;longLat(i)-dlat];

% plot ascending data
h = patch(px,py,longD(i),'edgecolor','none');


%% figure 1: Descending Density
% find descending portions of the orbits
i = find(diff(longLat) < 0);
i = union(i,i+1);
% Check validity
i = intersect(i,find( longD < 1 & longD > 0 & longValid == 0 ));

px = [longTime(i)-dtime;longTime(i)-dtime;longTime(i)+dtime;longTime(i)+dtime];
py = -[longLat(i)-dlat;longLat(i)+dlat;longLat(i)+dlat;longLat(i)-dlat] + 180;

% plot descending data
h = patch(px,py,longD(i),'edgecolor','none');

% Adjust, Annotate, and output plot
axis tight 
xlim(xlimit);
set(gca,'clim',[.1e-13,0.75e-12]);
cb = colorbar;
cb = colorbar;
set(gca,'layer','top','box','on')                    
hold on;
plot(xlim,90*ones(1,2),'k','linewidth',1);
text(30.45,108,sprintf('Descending (%d:%02d LT)',mod(floor(eqslt+12),24),floor(60*(eqslt-floor(eqslt)))),'fontsize',14,'color','w');
text(30.45,73,sprintf('Ascending (%d:%02d LT)',floor(eqslt),floor(60*(eqslt-floor(eqslt)))),'fontsize',14,'color','w');
ylim([-90,270]);
set(gca,'ytick',-90:45:270,'yticklabel',{'-90','-45','0','45','90','45','0','-45','-90'});
xlabel(sprintf('Day of Year (%04d)',yyyy(1)));
ylabel('Latitude (deg)');
set(get(cb,'ylabel'),'string','Density (kg/m^3)','rotation',-90,'verticalalignment','bottom','fontsize',14)
set(gcf,'paperposition',[0.25,3.5,8,4]);
print(gcf,'-depsc2','-painters','gracefo_LatVsTime_ascdesc');
%print(gcf,'-djpeg90','-r150','gracefo_LatVsTime_ascdesc');



%% Figure 2: Timeseries at 22 deg lat during ascending portions of orbit
figure

% sample data set at this latitude
sample_lat = 22;

% filter data
i = find( diff(longLat) > 0 );i = union(i,i+1); % find ascending data
i = intersect( i, find( (longLat(1:end-1) - sample_lat).*(longLat(2:end) - sample_lat) < 0 ) ); % find desired latitude
i = intersect( i, find( longD < 1 & longD > 0 & longValid == 0 ) ); % filter bad data
i = intersect( i, find( longTime >= 30 & longTime <= 40 ) ); % time period

% plot timeseries
plot(longTime(i),longD(i),'k','linewidth',1.5);

% Adjust, Annotate, and output plot
title(sprintf('GRACE-FO, %d^\\circ N, %d:%02d LT, ~500 km',sample_lat,floor(eqslt),floor(60*(eqslt-floor(eqslt)))));
xlabel(sprintf('Day of Year (%04d)',yyyy(1)));
ylabel('Neutral Density (kg/m^3)');
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'PaperPosition',[0.3611    4    7.7778    3])     
xlim(xlimit);
ylim([1,8]*1e-13)
print(gcf,'-depsc2',sprintf('gracefo_timeseries_%02ddeglat',sample_lat));
