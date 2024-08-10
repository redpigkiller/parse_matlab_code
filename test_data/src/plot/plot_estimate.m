
%% flags & variables
ylim_dfo = [1e-12 inf];
ylim_aoa = [-inf inf];
ylim_chl = [1e-5 inf];

cmap = ['b', 'k', 'r', "#3E9651", "#7E2F8E", "#D95319", "#4DBEEE"];
makr = ["*", "diamond", "o", "square", "x", "^"];
lins = ["-", "--", ":", "-."];

%% add path
addpath(genpath('src'));

%% user input
[answer, IsCancelled] = request_user_input_plot_estimate();

if IsCancelled;	return;	end

sim_name    = answer.sim_name;
sub_trial   = answer.exp_name;
dir_sim     = fullfile(sim_name, sub_trial);
dir_fig     = fullfile('figure', dir_sim);
if ~exist(dir_fig, 'dir');	mkdir(dir_fig);	end

%% load data
data_name   = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name);
data        = load(fullfile('results', dir_sim, data_name));

config_name = sprintf("%s_channel_config.est.mat", sim_name);
config      = load(fullfile('results', dir_sim, config_name)).config;

config_name = sprintf("%s_sim_config.est.mat", sim_name);
sim_config  = load(fullfile('results', dir_sim, config_name));

% unpack
Mc  = sim_config.Mc;
snr = sim_config.snr;
Q   = sim_config.Q;
P   = sim_config.P;
Nr  = sim_config.Nr;
N   = config.N;
M   = N / P;

mse_eps_q   = data.MSE_dfo;
mse_theta_q = data.MSE_aoa;
mse_h_q     = data.MSE_chl;

%% plot

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot DFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure(1);

for i_Q = length(Q) : -1 : 1
	semilogy(snr, mse_eps_q(i_Q, :), ...
		'DisplayName', ['$Q=', num2str(Q(i_Q)), '$'], ...
		'Color', cmap(i_Q), ...
		'Marker', makr(i_Q), ...
		'LineStyle', lins(1), ...
        'LineWidth', 1);
	hold on;
end
hold off;
xlabel('snr (dB)', 'interpreter', 'latex');
ylabel('MSE DFO', 'interpreter', 'latex');
ylim(ylim_dfo);
grid on;
	
if Nr > 1
    title_str = sprintf('$P=%d$, $N_r=%d$', P, Nr);
else
    title_str = sprintf('$P=%d$', P);
end
title(title_str, 'interpreter', 'latex');

legend('Interpreter', 'latex','Location', 'southwest');
		
% %%%%%%%%%% save figure %%%%%%%%%%
% add margin
a = annotation('rectangle', [0 0 1 1], 'Color', 'w');
% save to files
fname = sprintf('%s_%s_mse_dfo', sim_name, sub_trial);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot AOA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure(2);

for i_Q = length(Q) : -1 : 1
	semilogy(snr, mse_theta_q(i_Q, :), ...
		'DisplayName', ['$Q=', num2str(Q(i_Q)), '$'], ...
		'Color', cmap(i_Q), ...
		'Marker', makr(i_Q), ...
		'LineStyle', lins(1), ...
        'LineWidth', 1);
	hold on;
end
hold off;
xlabel('snr (dB)', 'interpreter', 'latex');
ylabel('MSE AoA', 'interpreter', 'latex');
ylim(ylim_aoa);
grid on;

if Nr > 1
    title_str = sprintf('$P=%d$, $N_r=%d$', P, Nr);
else
    title_str = sprintf('$P=%d$', P);
end
title(title_str, 'interpreter', 'latex');

legend('Interpreter', 'latex','Location', 'southwest');
		
% %%%%%%%%%% save figure %%%%%%%%%%
% add margin
a = annotation('rectangle', [0 0 1 1], 'Color', 'w');
% save to files
fname = sprintf('%s_%s_mse_aoa', sim_name, sub_trial);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot Channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure(3);

for i_Q = length(Q) : -1 : 1
	semilogy(snr, mse_h_q(i_Q, :), ...
		'DisplayName', ['$Q=', num2str(Q(i_Q)), '$'], ...
		'Color', cmap(i_Q), ...
		'Marker', makr(i_Q), ...
		'LineStyle', lins(1), ...
        'LineWidth', 1);
	hold on;
end
hold off;
xlabel('snr (dB)', 'interpreter', 'latex');
ylabel('MSE Channel', 'interpreter', 'latex');
ylim(ylim_chl);
grid on;
	
if Nr > 1
    title_str = sprintf('$P=%d$, $N_r=%d$', P, Nr);
else
    title_str = sprintf('$P=%d$', P);
end
title(title_str, 'interpreter', 'latex');

legend('Interpreter', 'latex', 'Location', 'southwest');

% %%%%%%%%%% print info %%%%%%%%%%
fprintf(1, "Plot info:\n");
fprintf(1, "\tsim_src:\t%s\n", sim_name);
fprintf(1, "\texp_name:\t%s\n", sub_trial);
fprintf(1, "\tMc\t=\t%d\n", Mc);
fprintf(1, "\tQ\t=\t[%s]\n", join(string(Q), ','));
fprintf(1, "\tsnr\t=\t[%s]\n", join(string(snr), ','));
fprintf(1, "\tP\t=\t%d\n", P);
fprintf(1, "\tNr\t=\t%d\n", Nr);

% %%%%%%%%%% save figure %%%%%%%%%%
% add margin
a = annotation('rectangle', [0 0 1 1], 'Color', 'w');
% save to files
fname = sprintf('%s_%s_mse_chl', sim_name, sub_trial);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

fprintf('figures are saved in <a href="matlab: winopen(''%s'')">%s</a>.\n', dir_fig, dir_fig);
