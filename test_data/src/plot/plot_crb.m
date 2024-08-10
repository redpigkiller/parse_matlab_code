
%% flags & variables

cmap = ['b', 'k', 'r', "#3E9651", "#7E2F8E", "#D95319", "#4DBEEE"];
makr = ["*", "diamond", "o", "square", "x", "^"];
lins = ["-", "--", ":", "-."];

%% add path
addpath(genpath('src'));

%% user input
[answer, IsCancelled] = request_user_input_plot_crb();

if IsCancelled;	return;	end

sim_name_crb    = answer.sim_name_crb;
sim_name_est    = answer.sim_name_est;
exp_name_crb    = answer.exp_name_crb{1};
exp_name_est    = answer.exp_name_est{1};
dataMetric      = answer.dataMetric;

dir_fig = fullfile('figure_crb', sim_name_crb);
if ~exist(dir_fig, 'dir');	mkdir(dir_fig);	end

%% load data
config_name = sprintf("%s_channel_config.crb.mat", sim_name_crb);
config = load(fullfile('results_crb', sim_name_crb, exp_name_crb, config_name)).config;
% CRB
data_name = sprintf("%s_avg_crb.crb.mat", sim_name_crb);
data_crb = load(fullfile('results_crb', sim_name_crb, exp_name_crb, data_name));
config_name = sprintf("%s_sim_config.crb.mat", sim_name_crb);
sim_config_crb = load(fullfile('results_crb', sim_name_crb, exp_name_crb, config_name));

% Estimated
flag_est_data = 1;

data_name = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name_est);
try
    data_est = load(fullfile('results', sim_name_est, exp_name_est, data_name));
catch
    flag_est_data = 0;
end
if flag_est_data == 1
    config_name = sprintf("%s_sim_config.est.mat", sim_name_est);
    sim_config_est = load(fullfile('results', sim_name_est, exp_name_est, config_name));
end

% unpack
Mc  = sim_config_crb.Mc;
snr = sim_config_crb.snr;
Q   = sim_config_crb.Q;
P   = sim_config_crb.P;
Nr  = sim_config_crb.Nr;
N   = config.N;

% check simulation config
if flag_est_data == 1
    if P ~= sim_config_est.P;       warning("unmatched P");     end
    if Nr ~= sim_config_est.Nr;     warning("unmatched Nr");    end
end

% retrieve data
if strcmpi(dataMetric, 'DFO')
    mse_crb = data_crb.avg_dfo_crb;
    if flag_est_data == 1
        mse_est = data_est.MSE_dfo;
    end
else
    mse_crb = data_crb.avg_aoa_crb;
    if flag_est_data == 1
        mse_est = data_est.MSE_aoa;
    end
end

% align the data if snr is different
[l_snr, unmatched_snr]   = intersects(snr, sim_config_est.snr);
crb_snr_idx    = arrcmp(l_snr, snr, 1);
est_snr_idx    = arrcmp(l_snr, sim_config_est.snr, 1);
snr            = l_snr;
mse_crb        = mse_crb(:, crb_snr_idx);
if flag_est_data == 1
    mse_est        = mse_est(:, est_snr_idx);
end

%% plot
fig = figure(1);

for q = length(Q) : -1 : 1
    % if q ~= 4
	semilogy(snr, mse_crb(q, :), ...
		'DisplayName', ['CRB ($Q=', num2str(Q(q)), '$)'], ...
		'Color', cmap(q), ...
		'Marker', makr(q), ...
		'LineStyle', lins(1), ...
        'LineWidth', 1);
    % else
    %     semilogy(snr,mse_crb(q, :),'Color', "w", 'DisplayName', '','visible','off')
    % end
	hold on;
    if flag_est_data == 1
        semilogy(snr, mse_est(q, :), ...
		    'DisplayName', ['ML ($Q=', num2str(Q(q)), '$)'], ...
		    'Color', cmap(q), ...
		    'Marker', makr(q), ...
		    'LineStyle', lins(2), ...
            'LineWidth', 1);
    end
end
hold off;
xlabel('snr (dB)', 'interpreter', 'latex');
ylabel('MSE DFO', 'interpreter', 'latex');
grid on;

if Nr > 1
    title_str = sprintf('%s (N=%d, P=%d, Nr=%d)', dataMetric, N, P, Nr);
else
    title_str = sprintf('%s (N=%d, P=%d)', dataMetric, N, P);
end
title_str = "$P=8$";
title(title_str, 'interpreter', 'latex');

lgd = legend('Interpreter', 'latex', 'Location', 'southwest', 'NumColumns', 2, 'Orientation', 'horizontal');

% %%%%%%%%%% print info %%%%%%%%%%
fprintf(1, "Plot info:\n");
fprintf(1, "\tsim_src:\t%s\n", sim_name);
fprintf(1, "\texp_name:\t%s\n", exp_name_est);
fprintf(1, "\tMc\t=\t%d\n", Mc);
fprintf(1, "\tQ\t=\t[%s]\n", join(string(Q), ','));
fprintf(1, "\tsnr\t=\t[%s]\n", join(string(snr), ','));
fprintf(1, "\tP\t=\t%d\n", P);
fprintf(1, "\tNr\t=\t%d\n", Nr);

% %%%%%%%%%% save figure %%%%%%%%%%
% add margin
a = annotation('rectangle', [0 0 1 1], 'Color', 'w');
% save to files
if ~exist(dir_fig, 'dir');	mkdir(dir_fig); end
fname = sprintf('%s_%s_crb_%s', sim_name_crb, exp_name_crb, dataMetric);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

fprintf('figures are saved in <a href="matlab: winopen(''%s'')">%s</a>.\n', dir_fig, dir_fig);

