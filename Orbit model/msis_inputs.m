%% format msis inputs
function [Ap_total, Ap_daily, F10_total, F81c] = msis_inputs(year, SOLdata, geomag_mat)
F10_total = SOLdata(4, SOLdata(1,:) == year);
F81c = SOLdata(5, SOLdata(1,:) == year);
Ap_mat = geomag_mat(geomag_mat(:,12) == year, 4:11);
Ap_mat = Ap_mat';
Ap_total = Ap_mat(:);
Ap_daily = mean(Ap_mat,1);
end
