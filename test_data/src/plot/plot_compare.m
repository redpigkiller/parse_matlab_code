%% user input
% lgd_str1 = 'ideal';
% lgd_str2 = 'estimation';
lgd_str1 = 'ML';
lgd_str2 = 'U-ESPRIT';
% lgd_str1 = 'ZF';
% lgd_str2 = 'MMSE';
% lgd_str1 = '$N_r=8$';
% lgd_str2 = '$N_r=1$';
plot_Q   = [1 2 3 4 5];
% plot_Q = 2;


%% flags & variables
ylim_dfo = [-inf inf];
ylim_aoa = [-inf inf];
ylim_chl = [1e-5 inf];
ylim_ber = [-inf inf];

cmap = ['b', 'k', 'r', "#3E9651", "#7E2F8E", "#D95319", "#4DBEEE"];
makr = ["*", "diamond", "o", "square", "x", "^"];
lins = ["-", "--", ":", "-."];

%% add path
addpath(genpath('src'));

%% user input
[answer, IsCancelled] = request_user_input();

if IsCancelled;	return;	end

sim_name1   = answer.sim_name1;
sim_name2   = answer.sim_name2;
exp_name1   = answer.exp_name1;
exp_name2   = answer.exp_name2;
plot_metric = answer.dataMetric;

dir_sim1    = fullfile(sim_name1, exp_name1);
dir_sim2    = fullfile(sim_name2, exp_name2);
dir_fig     = fullfile('figure', 'comparison');
if ~exist(dir_fig, 'dir');	mkdir(dir_fig);	end

%% check inputs
if ~strcmp(sim_name1, sim_name2)
    warning("The source of simulations are different! Make sure the " + ...
        "basic settings are the same.");
    fpritnf(1, "Use %s config by default...", sim_name1);
end

%% load data
if strcmp(plot_metric, 'MSE DFO')
    data_name   = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name1);
    data1       = load(fullfile('results', dir_sim1, data_name)).MSE_dfo;
    data_name   = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name2);
    data2       = load(fullfile('results', dir_sim2, data_name)).MSE_dfo;

elseif strcmp(plot_metric, 'MSE AoA')
    data_name   = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name1);
    data1       = load(fullfile('results', dir_sim1, data_name)).MSE_aoa;
    data_name   = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name2);
    data2       = load(fullfile('results', dir_sim2, data_name)).MSE_aoa;

elseif strcmp(plot_metric, 'MSE Channel')
    data_name   = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name1);
    data1       = load(fullfile('results', dir_sim1, data_name)).MSE_chl;
    data_name   = sprintf("%s_mse_dfo_aoa_channel.est.mat", sim_name2);
    data2       = load(fullfile('results', dir_sim2, data_name)).MSE_chl;

elseif strcmp(plot_metric, 'BER')
    data_name   = sprintf("%s_avg_ber.ber.mat", sim_name1);
    data1       = load(fullfile('results', dir_sim1, data_name)).avg_ber;
    data_name   = sprintf("%s_avg_ber.ber.mat", sim_name2);
    data2       = load(fullfile('results', dir_sim2, data_name)).avg_ber;

end

if strcmp(plot_metric, 'BER')
    config_name = sprintf("%s_channel_config.ber.mat", sim_name1);
    config1     = load(fullfile('results', dir_sim1, config_name)).config;
    config_name = sprintf("%s_channel_config.ber.mat", sim_name2);
    config2     = load(fullfile('results', dir_sim2, config_name)).config;
    
    config_name = sprintf("%s_sim_config.ber.mat", sim_name1);
    sim_config1 = load(fullfile('results', dir_sim1, config_name));
    config_name = sprintf("%s_sim_config.ber.mat", sim_name2);
    sim_config2 = load(fullfile('results', dir_sim2, config_name));
    
else
    config_name = sprintf("%s_channel_config.est.mat", sim_name1);
    config1     = load(fullfile('results', dir_sim1, config_name)).config;
    config_name = sprintf("%s_channel_config.est.mat", sim_name2);
    config2     = load(fullfile('results', dir_sim2, config_name)).config;
    
    config_name = sprintf("%s_sim_config.est.mat", sim_name1);
    sim_config1 = load(fullfile('results', dir_sim1, config_name));
    config_name = sprintf("%s_sim_config.est.mat", sim_name2);
    sim_config2 = load(fullfile('results', dir_sim2, config_name));
end

% unpack
Mc1     = sim_config1.Mc;
Mc2     = sim_config2.Mc;
snr1    = sim_config1.snr;
snr2    = sim_config2.snr;
Q1      = sim_config1.Q;
Q2      = sim_config2.Q;
if isfield(sim_config1, 'P');   P1      = sim_config1.P;
else;                           P1      = -1;end
if isfield(sim_config2, 'P');   P2      = sim_config2.P;
else;                           P2      = -1;end
Nr1     = sim_config1.Nr;
Nr2     = sim_config2.Nr;
N1      = config1.N;
N2      = config2.N;

%% plot
fig = figure(1);

% plot data1
if Nr1 ~= 1 || ~strcmpi(plot_metric, 'MSE AoA')
    for q1 = length(Q1) : -1 : 1
        if ~ismember(q1, plot_Q)
            continue;
        end
	    semilogy(snr1, data1(q1, :), ...
		    'DisplayName', [lgd_str1, ' ($Q=', num2str(Q1(q1)), '$)'], ...
		    'Color', cmap(q1), ...
		    'Marker', makr(q1), ...
		    'LineStyle', lins(1), ...
            'LineWidth', 1);
	    hold on;
    end
end

% semilogy(snr2,data1(5, :),'Color', "w", 'DisplayName', '','visible','off')

% plot data2
if Nr2 ~= 1 || ~strcmpi(plot_metric, 'MSE AoA')
    for q2 = length(Q2) : -1 : 1
        if ~ismember(q2, plot_Q)
            continue;
        end
        semilogy(snr2, data2(q2, :), ...
	        'DisplayName', [lgd_str2, ' ($Q=', num2str(Q2(q2)), '$)'], ...
	        'Color', cmap(q2), ...
	        'Marker', makr(q2), ...
	        'LineStyle', lins(2), ...
            'LineWidth', 1);
	    hold on;
    end
end
hold off;
xlabel('snr (dB)', 'interpreter', 'latex');
ylabel(plot_metric, 'interpreter', 'latex');
if strcmp(plot_metric, 'MSE DFO')
    ylim(ylim_dfo);
elseif strcmp(plot_metric, 'MSE AoA')
    ylim(ylim_aoa);
elseif strcmp(plot_metric, 'MSE Channel')
    ylim(ylim_chl);
elseif strcmp(plot_metric, 'BER')
    ylim(ylim_ber);
end
grid on;

title_str = "";
if N1 == N2
    title_str = title_str + "N=" + num2str(N1);
else
    title_str = title_str + "N1=" + num2str(N1) + ...
        " N2=" + num2str(N2);
end
if P1 == P2
    if P1 > 0
        title_str = title_str + ", " + "P=" + num2str(P1);
    end
else
    if P1 > 0
        title_str = title_str + ", " + "P1=" + num2str(P1);
    end
    if P2 > 0
        title_str = title_str + ", " + "P2=" + num2str(P2);
    end
end
if Nr1 == Nr2
    if Nr1 ~= 1
        title_str = title_str + ", " + "Nr=" + num2str(Nr1);
    end
else
    title_str = title_str + ", " + "Nr1=" + num2str(Nr1) + ...
        " Nr2=" + num2str(Nr2);
end

title_str = "$P=8$";
title(title_str, 'interpreter', 'latex');

% lgd = legend('Interpreter', 'latex', 'Location', 'northeast', 'NumColumns', 2, 'Orientation', 'vertical');
lgd = legend('Interpreter', 'latex', 'Location', 'southwest', 'NumColumns', 2, 'Orientation', 'vertical');
% lgd = legend('Interpreter', 'latex', 'Location', 'southwest');

% %%%%%%%%%% save figure %%%%%%%%%%
% add margin
a = annotation('rectangle', [0 0 1 1], 'Color', 'w');
% save to files
fname = sprintf('%s_%s_vs_%s_%s_%s', sim_name1, exp_name1, sim_name2, exp_name2, plot_metric);
path_full = fullfile(dir_fig, fname);
saveas(fig, strcat(path_full, '.fig'));
exportgraphics(fig, strcat(path_full, '.png'), 'Resolution', 300);
% delete margin
delete(a);

fprintf('figures are saved in <a href="matlab: winopen(''%s'')">%s</a>.\n', dir_fig, dir_fig);


%% functions
function [answer, IsCancelled] = request_user_input()
    % %%%%%%%%%%%%%%%%%%%%
    % get results
    % %%%%%%%%%%%%%%%%%%%%
    dir_res_now = get_subfolders_r('results', 2);

    if dir_res_now.count == 0
		answer = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation result.\n");
		return;
    end

    % %%%%%%%%%%%%%%%%%%%%
    % verify the results
    % %%%%%%%%%%%%%%%%%%%%
    dir_now = struct;

    count = 0;
    name = {};
    array = [];
    for i_arr = 1 : dir_res_now.count
        exp_struct = dir_res_now.array(i_arr);
        if exp_struct.count ~= 0
            count = count + 1;
            name{end+1} = dir_res_now.name{i_arr}; %#ok<AGROW> 
            array = [array; exp_struct]; %#ok<AGROW> 
        end
    end
    
    dir_now.count = count;
    dir_now.name = name;
    dir_now.array = array;

    % %%%%%%%%%%%%%%%%%%%%
    % default answer 1
    % %%%%%%%%%%%%%%%%%%%%
    defInd = struct;
    defInd.i_sim1 = 1;
    defInd.i_sim2 = 1;
    defInd.i_exp1 = 1;
    defInd.i_exp2 = 1;

	defAns = struct;
    defAns.dataMetric = 'MSE DFO';

    % %%%%%%%%%%%%%%%%%%%%
    % get MD5 of dir_now
    % %%%%%%%%%%%%%%%%%%%%
    md5 = GetMD5(dir_now, 'Array');

    % try to load last user input if nothing has changed
	defVal_path = fullfile('_tmp', 'defaultInput_plot_compare.mat');
    if isfile(defVal_path)
		defVal = load(defVal_path).defVal;

        md5_prev = defVal.md5;

        if strcmp(md5, md5_prev)
            defInd = defVal.defInd;
            defAns = defVal.defAns;
        else
            defAns.dataMetric = defVal.defAns.dataMetric;
        end
    end

    % %%%%%%%%%%%%%%%%%%%%
    % default answer 2
    % %%%%%%%%%%%%%%%%%%%%
	defAns.sim_name1 = dir_now.name(defInd.i_sim1);
	defAns.exp_name1 = dir_now.array(defInd.i_sim1).name(defInd.i_exp1);
	defAns.sim_name2 = dir_now.name(defInd.i_sim2);
	defAns.exp_name2 = dir_now.array(defInd.i_sim2).name(defInd.i_exp2);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%% ask for user input %%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	prompt = {};
	title = 'Plot setting';
	formats = {};
	options = struct;

    prompt(1,:) = {'Enter data type: ', 'dataMetric', []};
    formats(1,1).type = 'list';
    formats(1,1).style = 'radiobutton';
	formats(1,1).size = -1;
    formats(1,1).format = 'text';
    formats(1,1).items = {'MSE DFO' 'MSE AoA' 'MSE Channel' 'BER'};
    formats(1,1).span = [1, 2];

	prompt(2,:) = {'Enter simulation data name 1:', 'sim_name1', []};
	formats(2,1).type = 'list';
	formats(2,1).style = 'popupmenu';
	formats(2,1).size = -1;
	formats(2,1).format = 'text';
	formats(2,1).items = dir_now.name;
	formats(2,1).callback = @(~,~,h,k)(set(h(k+2), ...
		'String', dir_now.array(get(h(k),'Value')).name, ...
		'Value', 1));

    prompt(3,:) = {'Enter simulation data name 2:', 'sim_name2', []};
	formats(2,2).type = 'list';
	formats(2,2).style = 'popupmenu';
	formats(2,2).size = -1;
	formats(2,2).format = 'text';
	formats(2,2).items = dir_now.name;
	formats(2,2).callback = @(~,~,h,k)(set(h(k+2), ...
		'String', dir_now.array(get(h(k),'Value')).name, ...
		'Value', 1));

	prompt(4,:) = {'Experiment 1', 'exp_name1',[]};
	formats(3,1).type = 'list';
	formats(3,1).style = 'popupmenu';
	formats(3,1).size = -1;
	formats(3,1).format = 'text';
	formats(3,1).items = dir_now.array(defInd.i_sim1).name;

	prompt(5,:) = {'Experiment 2', 'exp_name2',[]};
	formats(3,2).type = 'list';
	formats(3,2).style = 'popupmenu';
	formats(3,2).size = -1;
	formats(3,2).format = 'text';
	formats(3,2).items = dir_now.array(defInd.i_sim2).name;
	
	options.Resize = 'off';
	options.AlignControls = 'on';
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
	if IsCancelled == 0
		answer.sim_name1 = answer.sim_name1{1};
		answer.sim_name2 = answer.sim_name2{1};
		answer.exp_name1 = answer.exp_name1{1};
		answer.exp_name2 = answer.exp_name2{1};

        % record defAns
        defVal = struct;
        defVal.md5 = md5;
        % compute ind
        ind = struct;
        ind.i_sim1 = find(ismember(dir_now.name, answer.sim_name1));
        ind.i_sim2 = find(ismember(dir_now.name, answer.sim_name2));
        ind.i_exp1 = find(ismember(dir_now.array(ind.i_sim1).name, answer.exp_name1));
        ind.i_exp2 = find(ismember(dir_now.array(ind.i_sim2).name, answer.exp_name2));
        % save to defVal
        defVal.defInd = ind;
        defVal.defAns = answer;
        if ~exist('_tmp', 'dir');	mkdir('_tmp'); end
        save(defVal_path, 'defVal');
	end
end
