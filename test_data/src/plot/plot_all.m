
%% flags & variables
cmap = ['b', 'k', 'r', "#3E9651", "#7E2F8E", "#D95319", "#4DBEEE"];
makr = ["*", "diamond", "o", "square", "x", "^"];
lins = ["-", "--", ":", "-."];

%% add path
addpath(genpath('src'));

%% user input
[answer, IsCancelled] = request_user_input1();

if IsCancelled == 1
    return
end

sim_src     = answer.sim_src{1};
sim_exp     = answer.sim_exp;
plot_metric = answer.dataMetric;

[name_answer, IsCancelled] = request_user_input2(sim_exp);
plot_Q      = name_answer.plot_Q;

% pre-load
est_data_name = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_src);
ber_data_name = sprintf("%s_avg_ber.ber.mat", sim_src);

est_config_name     = sprintf("%s_channel_config.est.mat", sim_src);
est_sim_config_name = sprintf("%s_sim_config.est.mat", sim_src);

ber_config_name     = sprintf("%s_channel_config.ber.mat", sim_src);
ber_sim_config_name = sprintf("%s_sim_config.ber.mat", sim_src);

dir_fig     = fullfile('figure', 'comparison');
if ~exist(dir_fig, 'dir');	mkdir(dir_fig);	end

%% load data
n_exp = length(sim_exp);

sim_exp_name    = cell(n_exp, 1);
data            = cell(n_exp, 1);
config          = cell(n_exp, 1);
sim_config      = cell(n_exp, 1);

for i_exp = 1 : n_exp
    sim_exp_name{i_exp} = getfield(name_answer, sprintf('sim_exp_name_%d', i_exp));

    switch true
        case (strcmpi(plot_metric, 'MSE DFO') == 1)
            dir_sim = fullfile(sim_src, sim_exp{i_exp});
            data{i_exp}         = load(fullfile('results', dir_sim, est_data_name)).MSE_dfo;
            config{i_exp}       = load(fullfile('results', dir_sim, est_config_name)).config;
            sim_config{i_exp}   = load(fullfile('results', dir_sim, est_sim_config_name));

        case (strcmpi(plot_metric, 'MSE AOA') == 1)
            dir_sim = fullfile(sim_src, sim_exp{i_exp});
            data{i_exp} = load(fullfile('results', dir_sim, est_data_name)).MSE_aoa;
            config{i_exp}       = load(fullfile('results', dir_sim, est_config_name)).config;
            sim_config{i_exp}   = load(fullfile('results', dir_sim, est_sim_config_name));

        case (strcmpi(plot_metric, 'MSE Channel') == 1)
            dir_sim = fullfile(sim_src, sim_exp{i_exp});
            data{i_exp} = load(fullfile('results', dir_sim, est_data_name)).MSE_chl;
            config{i_exp}       = load(fullfile('results', dir_sim, est_config_name)).config;
            sim_config{i_exp}   = load(fullfile('results', dir_sim, est_sim_config_name));

        case (strcmpi(plot_metric, 'BER') == 1)
            dir_sim = fullfile(sim_src, sim_exp{i_exp});
            data{i_exp} = load(fullfile('results', dir_sim, ber_data_name)).avg_ber;
            config{i_exp}       = load(fullfile('results', dir_sim, ber_config_name)).config;
            sim_config{i_exp}   = load(fullfile('results', dir_sim, ber_sim_config_name));
    end

end

%% plot
fig = figure(1);

% plot data
for i_exp = 1 : n_exp
    Q   = sim_config{i_exp}.Q;
    snr = sim_config{i_exp}.snr;
    for i_Q = length(Q) : -1 : 1
        if ~ismember(Q(i_Q), plot_Q)
            continue;
        end
	    semilogy(snr, data{i_exp}(i_Q, :), ...
		    'DisplayName', [sim_exp_name{i_exp}, ' (Q=', num2str(Q(i_Q)), ')'], ...
		    'Color', cmap(i_Q), ...
		    'Marker', makr(i_Q), ...
		    'LineStyle', lins(i_exp), ...
            'LineWidth', 1);
	    hold on;
    end
end

hold off;
xlabel('snr');
ylabel(plot_metric);
grid on;

title_str = "comparison";

title(title_str);

lgd = legend('Location', 'southwest');

% %%%%%%%%%% save figure %%%%%%%%%%
% add margin
a = annotation('rectangle', [0 0 1 1], 'Color', 'w');
% save to files
fname = sprintf('%s_comparison_on_%s', sim_src, plot_metric);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

%% functions
function [answer, IsCancelled] = request_user_input1()
    % %%%%%%%%%%%%%%%%%%%%
    % get dump
    % %%%%%%%%%%%%%%%%%%%%
    dir_dmp_now = get_subfolders_r('results', 2);

    if dir_dmp_now.count == 0
		answer = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation.\n");
		return;
    end

    % %%%%%%%%%%%%%%%%%%%%
    % verify the results
    % %%%%%%%%%%%%%%%%%%%%
    dir_now = struct;

    count = 0;
    name = {};
    array = [];
    for i_arr = 1 : dir_dmp_now.count
        exp_struct = dir_dmp_now.array(i_arr);
        if exp_struct.count ~= 0
            count = count + 1;
            name{end+1} = dir_dmp_now.name{i_arr}; %#ok<AGROW> 
            array = [array; exp_struct]; %#ok<AGROW> 
        end
    end
    
    dir_now.count = count;
    dir_now.name = name;
    dir_now.array = array;

    % %%%%%%%%%%%%%%%%%%%%
    % default answer
    % %%%%%%%%%%%%%%%%%%%%
	defAns = struct;
	defAns.sim_src = dir_now.name(1);
	defAns.sim_exp = dir_now.array(1).name(1);

	% ask for user input
	prompt = {};
	title = 'Remove simulation data';
	formats = {};
	options = struct;

	prompt(1,:) = {'Simulation source:', 'sim_src', []};
	formats(1,1).type = 'list';
	formats(1,1).style = 'listbox';
	formats(1,1).size = [200 75];
	formats(1,1).format = 'text';
	formats(1,1).items = dir_now.name;
	formats(1,1).callback = @(~,~,h,k)(set(h(k+1), ...
		'String', dir_now.array(get(h(k),'Value')).name, ...
		'Value', 1));
	
	prompt(2,:) = {'Simulation experiment:', 'sim_exp',[]};
	formats(2,1).type = 'list';
	formats(2,1).style = 'listbox';
	formats(2,1).size = [200 150];
	formats(2,1).format = 'text';
	formats(2,1).items = dir_now.array(1).name;
    formats(2,1).limits = [0 2]; % multi-select

    prompt(3,:) = {'Enter data type: ', 'dataMetric', []};
    formats(3,1).type = 'list';
    formats(3,1).style = 'radiobutton';
	formats(3,1).size = -1;
    formats(3,1).format = 'text';
    formats(3,1).items = {'MSE DFO' 'MSE AOA' 'MSE Channel' 'BER'};
%     formats(3,1).span = [1, 2];

    prompt(4,:) = {'Continue to set the plot Q and name of legends...', [], []};
	formats(4,1).type = 'text';
    formats(4,1).size = [-1 0];

	options.Resize = 'off';
	options.AlignControls = 'on';
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
end

function [answer, IsCancelled] = request_user_input2(sim_exp)
    n_exp = length(sim_exp);

	% ask for user input
	prompt  = cell(n_exp+1, 3);
	title   = 'sim_exp name';
	formats = struct;
	defAns  = struct;
	options = struct;

    prompt(1,:) = {"plot Q", 'plot_Q',[]};
	formats(1,1).type = 'edit';
	formats(1,1).size = [-1, 20];
	formats(1,1).format = 'vector';

    for i_exp = 1 : n_exp
        data_name = sprintf('sim_exp_name_%d', i_exp);
        prompt(i_exp+1,:) = {sprintf("Input name for %s:", sim_exp{i_exp}), data_name,[]};
        formats(i_exp+1,1).type = 'edit';
	    formats(i_exp+1,1).size = -1;
	    formats(i_exp+1,1).format = 'text';
    end

	options.Resize          = 'off';
	options.AlignControls   = 'on';
	options.Interpreter     = 'none';
    options.CancelButton    = 'off';
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
end