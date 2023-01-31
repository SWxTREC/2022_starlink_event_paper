%% read indices
clc
clearvars
data_mat = readmatrix('SOLFSMY_2022.txt');  
SOLdata = data_mat(:,1:11)';
% year, doy, jd, F10, F81c, S10, S81c, M10, M81c, Y10, Y81c

data_mat = readmatrix('DTCFILE_2022.txt');
DTCdata = data_mat(:,2:end)';

%%
data_mat = readmatrix('SOLRESAP_2022.txt');
% Ap index = 4:11