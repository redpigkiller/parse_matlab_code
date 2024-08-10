N       = 1024;
plot_Q  = 2;

sim_name = 'thesis_final_6000';
exp_name = { "P8_Nr1_uniESPRIT_PINV_PINV_ignore_mmse_Nratio8"
    "P8_Nr1_uniESPRIT_PINV_PINV_noise1_estimate_MMSE_Nratio8"
    "P8_Nr1_uniESPRIT_PINV_PINV_ideal_mmse_Nratio8"
    };
lgd_str  = {"ignore", "estimation", "ideal"};

cmap = ['b', 'k', 'r', "#3E9651", "#7E2F8E", "#D95319", "#4DBEEE"];
makr = ["*", "diamond", "o", "square", "x", "^"];
lins = ["-", "--", ":", "-."];

%% add path
addpath(genpath('src'));

%% load data
n_data = size(exp_name, 1);

plot_data    = cell(n_data, 1);
plot_config  = cell(n_data, 1);

for i_data = 1 : n_data
    dir_sim     = fullfile(sim_name, exp_name{i_data});
    
    data_name = sprintf("%s_avg_ber.ber.mat", sim_name);
    data      = load(fullfile('results', dir_sim, data_name));
    config_name = sprintf("%s_sim_config.ber.mat", sim_name);
    sim_config  = load(fullfile('results', dir_sim, config_name));
    plot_data{i_data}    = data.avg_ber;

    plot_config{i_data}  = sim_config;
end

dir_fig = fullfile('figure', 'thesis_ex3.6.4');
if ~exist(dir_fig, 'dir');	mkdir(dir_fig);	end


%% plot
fig = figure(1);

for i_data = 1 : n_data
    % unpack
    avg_ber     = plot_data{i_data};
    sim_config  = plot_config{i_data};

    Mc  = sim_config.Mc;
    snr = sim_config.snr;
    Q   = sim_config.Q;
    Nr  = sim_config.Nr;

    for q = length(Q) : -1 : 1
        if ~ismember(plot_Q, q)
            continue;
        end

	    semilogy(snr, avg_ber(q, :), ...
		    'DisplayName', [lgd_str{i_data}], ...
		    'Color', cmap(i_data), ...
		    'Marker', makr(i_data), ...
		    'LineStyle', lins(1), ...
            'LineWidth', 1);
	    hold on;
    end
end
hold off;
xlabel('snr (dB)', 'interpreter', 'latex');
ylabel('BER', 'interpreter', 'latex');    
grid on;

title(sprintf('Q = %d, MMSE equalizer', plot_Q), 'interpreter', 'latex');

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
fname = sprintf('%s_Q=%d_%s', sim_name, plot_Q, plot_metric);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

fprintf('figures are saved in <a href="matlab: winopen(''%s'')">%s</a>.\n', dir_fig, dir_fig);


