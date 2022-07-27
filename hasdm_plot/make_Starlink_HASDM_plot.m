load StarlinkFeb2022_hasdm.mat
set(0,'defaultaxesfontsize',14)
figure
plot(Doy+Sec/86400,hasdm_density,'k','linewidth',1.5)
title('HASDM, 22^\circ N, 12:45 LT, 210 km');
xlabel('Day of Year 2022');ylabel('Neutral Density (kg/m^3)');
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'PaperPosition',[0.3611    4    7.7778    3])     
ylim([1.7,2.5]*1e-10)
print(gcf,'-depsc2','StarlinkFeb2022_hasdm');
