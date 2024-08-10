N       = 1024;
plot_Q  = 4;
% plot_metric = 'DFO';
% plot_metric = 'AoA';
% plot_metric = 'CHL';
plot_metric = 'BER';

sim_name = 'thesis_final';
if strcmpi(plot_metric, 'DFO') || strcmpi(plot_metric, 'CHL')
    exp_name = { "P8_Nr1_uniESPRIT_PINV_PINV"
        "P4_Nr1_uniESPRIT_PINV_PINV"
        "P8_Nr2_uniESPRIT_PINV_PINV"
        "P4_Nr2_uniESPRIT_PINV_PINV"
        "P8_Nr4_uniESPRIT_PINV_PINV"
        "P4_Nr4_uniESPRIT_PINV_PINV"
        "P8_Nr8_uniESPRIT_PINV_PINV"
        "P4_Nr8_uniESPRIT_PINV_PINV"
        "P8_Nr16_uniESPRIT_PINV_PINV"
        "P4_Nr16_uniESPRIT_PINV_PINV"
        };
elseif strcmpi(plot_metric, 'AoA')
    exp_name = { "P8_Nr2_uniESPRIT_PINV_PINV"
        "P4_Nr2_uniESPRIT_PINV_PINV"
        "P8_Nr4_uniESPRIT_PINV_PINV"
        "P4_Nr4_uniESPRIT_PINV_PINV"
        "P8_Nr8_uniESPRIT_PINV_PINV"
        "P4_Nr8_uniESPRIT_PINV_PINV"
        "P8_Nr16_uniESPRIT_PINV_PINV"
        "P4_Nr16_uniESPRIT_PINV_PINV"
        };
else
    exp_name = { "P8_Nr1_uniESPRIT_PINV_PINV_noise1_estimate_MMSE"
        "P4_Nr1_uniESPRIT_PINV_PINV_noise1_estimate_MMSE"
        "P8_Nr2_uniESPRIT_PINV_PINV_ZF"
        "P4_Nr2_uniESPRIT_PINV_PINV_ZF"
        "P8_Nr4_uniESPRIT_PINV_PINV_ZF"
        "P4_Nr4_uniESPRIT_PINV_PINV_ZF"
        "P8_Nr8_uniESPRIT_PINV_PINV_noise1_estimate_ZF"
        "P4_Nr8_uniESPRIT_PINV_PINV_noise1_estimate_ZF"
        "P8_Nr16_uniESPRIT_PINV_PINV_noise1_estimate_ZF"
        "P4_Nr16_uniESPRIT_PINV_PINV_noise1_estimate_ZF"
        };
end


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
    
    if strcmpi(plot_metric, 'DFO')
        data_name = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name);
        data      = load(fullfile('results', dir_sim, data_name));
        config_name = sprintf("%s_sim_config.est.mat", sim_name);
        sim_config  = load(fullfile('results', dir_sim, config_name));
        plot_data{i_data}    = data.MSE_dfo;    
    elseif strcmpi(plot_metric, 'AoA')
        data_name = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name);
        data      = load(fullfile('results', dir_sim, data_name));
        config_name = sprintf("%s_sim_config.est.mat", sim_name);
        sim_config  = load(fullfile('results', dir_sim, config_name));
        plot_data{i_data}    = data.MSE_aoa;    
    elseif strcmpi(plot_metric, 'CHL')
        data_name = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name);
        data      = load(fullfile('results', dir_sim, data_name));
        config_name = sprintf("%s_sim_config.est.mat", sim_name);
        sim_config  = load(fullfile('results', dir_sim, config_name));
        plot_data{i_data}    = data.MSE_chl;  
    else
        data_name = sprintf("%s_avg_ber.ber.mat", sim_name);
        data      = load(fullfile('results', dir_sim, data_name));
        config_name = sprintf("%s_sim_config.ber.mat", sim_name);
        sim_config  = load(fullfile('results', dir_sim, config_name));
        plot_data{i_data}    = data.avg_ber;
    end

    plot_config{i_data}  = sim_config;
end

dir_fig = fullfile('figure', 'thesis_ex4.5.2');
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

        % plot color settings
        if mod(i_data, 2) == 0; i_lins = 2; else;   i_lins = 1; end
        i_color = floor((i_data-1)/2)+1;

        if mod(i_data, 2) == 0; P = 4; else;   P = 8; end
    
	    semilogy(snr, avg_ber(q, :), ...
		    'DisplayName', ['$P=', num2str(P), '$, $N_r=', num2str(Nr), '$'], ...
		    'Color', cmap(i_color), ...
		    'Marker', makr(i_color), ...
		    'LineStyle', lins(i_lins), ...
            'LineWidth', 1);
	    hold on;
    end
end
hold off;
xlabel('snr (dB)', 'interpreter', 'latex');
if strcmpi(plot_metric, 'DFO')
    ylabel('MSE DFO', 'interpreter', 'latex');
elseif strcmpi(plot_metric, 'CHL')
    ylabel('MSE channel', 'interpreter', 'latex');
elseif strcmpi(plot_metric, 'AoA')
    ylabel('MSE AoA', 'interpreter', 'latex');
else
    ylabel('BER', 'interpreter', 'latex');    
end
grid on;

title(sprintf('Q = %d', plot_Q), 'interpreter', 'latex');

lgd = legend('Interpreter', 'latex', 'Location', 'southwest', 'NumColumns', 2, 'Orientation', 'horizontal');

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


