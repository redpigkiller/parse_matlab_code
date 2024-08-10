sim_name = 'thesis_final';
exp_name = 'P8_Nr1_ML_TIKH_PINV';

Mc = 5000;
snr = -10 : 10 : 40;
Q = 1 : 5;
N = 1024;
M = 128;
Lch = 23;

dir_fig = fullfile('figure');

%% add path
addpath(genpath('src'));

src_data_path = fullfile('_gendata', ...
    sim_name);

% est_data_path = fullfile('_simResult_est_dfo_chl', ...
%     sim_name, ...
%     exp_name);
est_data_path = fullfile('_dump', ...
    sim_name, ...
    exp_name);

n_snr   = length(snr);
n_Q     = length(Q);

r_MSE_Chlmtx_Mc = zeros(n_Q, n_snr, Mc);

% parfor_progress(Mc);

MSE_DFO = zeros(Mc, 1);

parfor i_Mc = 1 : Mc
% for i_Mc = 1 : Mc
    i_Mc
    % load src data
    filename = sprintf('%s_%03d.mat', sim_name, i_Mc);
    data = load(fullfile(src_data_path, filename));

    dfoVal = data.dfoVal;
    aoaVal = data.aoaVal;
    chlVal = data.chlVal;

    % load est data
    filename = sprintf('%s_mse_%03d.est.mat', sim_name, i_Mc);
    data = load(fullfile(est_data_path, filename));

    est_MSEdfo = data.est_MSE_dfo;
%     est_chlVal = data.est_chlVal;
    MSE_DFO(i_Mc) = est_MSEdfo(3, 6);
    
end
histogram(MSE_DFO);
return
%% plot

MSE_Chlmtx = sum(r_MSE_Chlmtx_Mc, 3) ./ Mc;

cmap = ['b', 'k', 'r', 'g', 'c', 'm', 'y'];
makr = ['o', '+', '*', '.', 'x', '-'];
lins = ["-", "--", ":", "-."];

fig = figure(1);
for q = 1 : length(Q)
	semilogy(snr, MSE_Chlmtx(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(3), ...
		'LineStyle', lins(1));
	hold on;
end

hold off;
xlabel('snr');
ylabel('MSE Channel Matrix');
grid on;

title_str = sprintf('N=%d, M=%d, Mc=%d', N, M, Mc);
title(title_str);

legend('Location', 'southwest');

% %%%%%%%%%% save figure %%%%%%%%%%
% add margin
a = annotation('rectangle', [0 0 1 1], 'Color', 'w');
% save to files
fname = sprintf('%s_%s_mse_chlmtx', sim_name, exp_name);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

