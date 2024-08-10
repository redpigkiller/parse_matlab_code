
%% flags & variables
cmap = ['b', 'k', 'r', "#3E9651", "#7E2F8E", "#D95319", "#4DBEEE"];
makr = ["*", "diamond", "o", "square", "x", "^"];
lins = ["-", "--", ":", "-."];

%% add path
addpath(genpath('src'));

%% user input
[answer, IsCancelled] = request_user_input_plot_ber();

if IsCancelled;	return;	end

sim_name = answer.sim_name;
sub_trial = answer.exp_name;
dir_sim = fullfile(sim_name, sub_trial);
dir_fig = fullfile('figure', dir_sim);
if ~exist(dir_fig, 'dir');	mkdir(dir_fig);	end

%% load data
data_name = sprintf("%s_avg_ber.ber.mat", sim_name);
data = load(fullfile('results', dir_sim, data_name));
config_name = sprintf("%s_channel_config.ber.mat", sim_name);
config = load(fullfile('results', dir_sim, config_name)).config;
config_name = sprintf("%s_sim_config.ber.mat", sim_name);
sim_config = load(fullfile('results', dir_sim, config_name));

% unpack
Mc  = sim_config.Mc;
snr = sim_config.snr;
Q   = sim_config.Q;
Nr  = sim_config.Nr;
N   = config.N;
avg_ber = data.avg_ber;

%% plot
fig = figure(1);

for q = length(Q) : -1 : 1
	semilogy(snr, avg_ber(q, :), ...
		'DisplayName', ['$Q=', num2str(Q(q)), '$'], ...
		'Color', cmap(q), ...
		'Marker', makr(q), ...
		'LineStyle', lins(1), ...
        'LineWidth', 1);
	hold on;
end
hold off;
xlabel('snr (dB)', 'interpreter', 'latex');
ylabel('BER', 'interpreter', 'latex');
grid on;

if Nr > 1
    title_str = sprintf('$N=%d$, $Nr=%d$', N, Nr);
else
    title_str = sprintf('$N=%d$', N);
end
title_str = "$P=8$, MMSE equalizer";
title(title_str, 'interpreter', 'latex');

lgd = legend('Interpreter', 'latex', 'Location', 'southwest');

% %%%%%%%%%% print info %%%%%%%%%%
fprintf(1, "Plot info:\n");
fprintf(1, "\tsim_src:\t%s\n", sim_name);
fprintf(1, "\texp_name:\t%s\n", sub_trial);
fprintf(1, "\tMc\t=\t%d\n", Mc);
fprintf(1, "\tQ\t=\t[%s]\n", join(string(Q), ','));
fprintf(1, "\tsnr\t=\t[%s]\n", join(string(snr), ','));
fprintf(1, "\tNr\t=\t%d\n", Nr);

% %%%%%%%%%% save figure %%%%%%%%%%
% add margin
a = annotation('rectangle', [0 0 1 1], 'Color', 'w');
% save to files
fname = sprintf('%s_%s_ber', sim_name, sub_trial);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

fprintf('figures are saved in <a href="matlab: winopen(''%s'')">%s</a>.\n', dir_fig, dir_fig);
