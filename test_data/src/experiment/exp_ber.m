%
% Generate simulation trials for BER
%
clc;
clear;
close;

%% flags & variables
SHOW_PROGRESS = 0;

SEED = 506;

glob_flag = struct;	% global flags
algo_flag = struct;	% algorithm flags

% %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
glob_flag.add_noise = 1;
glob_flag.assume_known_channel = 0;
glob_flag.use_correct_dfo = 0;
glob_flag.use_correct_chl = 0;

% %%%%%%%%%%%%%%%%%%%% flags for BER %%%%%%%%%%%%%%%%%%%%
% ==================== simulation flags ====================
% the position of null subcarriers
% 0: (1 : null_subcarrier) and (end - null_subcarrier : end) set to 0
% 1: (end/2 - null_subcarrier : end/2 + null_subcarrier) set to 0
algo_flag.position = 1;
algo_flag.null_subcarrier = 0;

% ==================== equalizer flags ====================
% equalizer used to estimate data:
%   ZF
%   MMSE
algo_flag.equalizer = 'MMSE';
algo_flag.assume_no_dfo = 0;

% -------------------- [equalizer] options --------------------
if strcmpi(algo_flag.equalizer, 'ZF')
    % no options

elseif strcmpi(algo_flag.equalizer, 'MMSE')
    % 0: estimate snr
    % 1: assume RX knows snr
    algo_flag.equ_opt.use_correct_snr = 0;
    % 0: estimate snr for each data block
    % 1: estimate snr once based on all data blocks
    algo_flag.equ_opt.compute_snr_each_datablock = 0;

else
    error("wrong equalizer");
end

if glob_flag.assume_known_channel == 1
	glob_flag.use_correct_dfo = 1;
	glob_flag.use_correct_chl = 1;
end

%% add path
addpath(genpath('src'));

%% simulation settings
% no. of Monte Carlo
default_sim_Mc = 10;
% range of snr
default_sim_snr = 40;
% no. of path
default_sim_Q = 1:4;
% default experiment name
default_exp_name = 'first try';

%% user input
[answer, IsCancelled] = request_user_input(default_sim_Mc, default_exp_name);

if IsCancelled;	return;	end

% 1. load flag
glob_flag.use_gpu = answer.flag_use_gpu;
glob_flag.use_parfor = answer.flag_use_parfor;
% skip simulated i_Mc runs
glob_flag.skip_simulated = answer.flag_skip_simulated;
% only collect, no simulation
glob_flag.collect_mode = answer.flag_collect_mode;

% 2. get infomation
sim_name = answer.sim_src_name;
exp_name = answer.exp_src_name;
default_sim_Mc = answer.Mc;
sub_trial = answer.exp_name;
dir_exp = fullfile(sim_name, exp_name);
dir_sim = fullfile(sim_name, sub_trial);

%% load data
% get config
config_name = sprintf('%s_config.mat', sim_name);
config = load(fullfile('_gendata', sim_name, config_name)).config;
gen_config_name = sprintf('%s_gen_config.mat', sim_name);
gen_config = load(fullfile('_gendata', sim_name, gen_config_name));
est_config_name = sprintf('%s_sim_config.chl.mat', sim_name);
est_config = load(fullfile('_simResult_est_dfo_chl', dir_exp, est_config_name));

%	get Mc
Mc = min([default_sim_Mc, gen_config.Mc, est_config.Mc]);

%	get snr
[snr, unmatched_snr] = intersects(default_sim_snr, gen_config.snr, est_config.snr);
if unmatched_snr;	warning("different simulation snr setting");	end
snr_gendata_idx = cmp_arr(snr, gen_config.snr, 1);
snr_estchal_idx = cmp_arr(snr, est_config.snr, 1);

%	get Q
[Q, unmatched_Q] = intersects(default_sim_Q, gen_config.Q, est_config.Q);
if unmatched_Q;	warning("different simulation Q setting");	end
Q_gendata_idx = cmp_arr(Q, gen_config.Q, 1);
Q_estchal_idx = cmp_arr(Q, est_config.Q, 1);

%% save directory
dir_dump = fullfile('_exp_dump', dir_sim);
if ~exist(dir_dump, 'dir');	mkdir(dir_dump); end

%% auxiliary variables
n_snr = length(snr);
n_Q = length(Q);

%% pack variables
sim_Mc_info = struct;
% simulation info
sim_Mc_info.snr = snr;
sim_Mc_info.Q = Q;
sim_Mc_info.config = config;
% simulation flags
sim_Mc_info.glob_flag = glob_flag;
sim_Mc_info.algo_flag = algo_flag;
% directory of data
sim_Mc_info.dir_gendata = fullfile('_gendata', sim_name, sim_name);
sim_Mc_info.dir_estdata = fullfile('_simResult_est_dfo_chl', dir_exp, sim_name);
% indexes mapping
sim_Mc_info.snr_gendata_idx = snr_gendata_idx;
sim_Mc_info.snr_estchal_idx = snr_estchal_idx;
sim_Mc_info.Q_gendata_idx = Q_gendata_idx;
sim_Mc_info.Q_estchal_idx = Q_estchal_idx;

%% declare record var.
r_avg_ber_Mc = zeros(n_Q, n_snr, Mc);

%% run simulation
if glob_flag.use_parfor ~= 0
	poolobj = startparpool();
	parforArg = poolobj.NumWorkers;
else
	parforArg = 0;
end

% set random seed: rng(SEED); for each worker
% otherwise, the random number generated for each worker may not be same
prngsc = parallel.pool.Constant(RandStream('Threefry', 'Seed', SEED));
psimset = parallel.pool.Constant(sim_Mc_info);

% determine the sim_mode: collect, skip simulated
sim_mode = 2 * glob_flag.collect_mode + glob_flag.skip_simulated;
valid_Mc = ones(Mc, 1);

% start simulation
fprintf(1, "Start simulation with exp_name: %s\n", sub_trial);
if SHOW_PROGRESS == 1;  parfor_progress(Mc);    end
tic;
% parfor (i_Mc = 1 : Mc, parforArg)
for (i_Mc = 1 : Mc)
    % set random seed
    stream = prngsc.Value;
    stream.Substream = i_Mc;
    RandStream.setGlobalStream(stream);

    % path name for i_Mc
    path_full_name = sprintf('%s_r_avg_ber_%03d.ber.mat', sim_name, i_Mc);
	path_full = fullfile(dir_dump, path_full_name);
	
	% check mode
	if sim_mode >= 2
		if isfile(path_full)
		else
			valid_Mc(i_Mc) = 0;
		end

    elseif sim_mode >= 1 && isfile(path_full)
	else
		% run simulation
		[r_eZF_Mc, r_eMMSE_Mc, r_eCORR_Mc] = exp_one_Mc_ber(i_Mc, psimset.Value);

		% record
% 		r_avg_ber_Mc(:, :, i_Mc) = r_avg_ber;

		% dump record
		parsave(path_full, r_eZF_Mc, r_eMMSE_Mc, r_eCORR_Mc, Q, snr);
	end
	
	if SHOW_PROGRESS == 1;   parfor_progress;    end
end

if SHOW_PROGRESS == 1;  parfor_progress(0); end
fprintf(1, "done!\n");

%% show info
sim_exe_time = toc;
fprintf(1, "\nexecution time: %s\n", show_duration(sim_exe_time));



%% Functions
function [answer, IsCancelled] = request_user_input(default_Mc, default_exp_name)
    % %%%%%%%%%%%%%%%%%%%%
    % get gen_cases
    % %%%%%%%%%%%%%%%%%%%%
    dir_gen_now = get_subfolders_r('_gendata', 2);

	if dir_gen_now.count == 0
		answer = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation source.\n" + ...
            "Please execute gen_cases.m first!");
		return;
	end

    % %%%%%%%%%%%%%%%%%%%%
    % get simulation result for estimated DFO and Channel
    % %%%%%%%%%%%%%%%%%%%%
    dir_sim_now = get_subfolders_r('_simResult_est_dfo_chl', 3);

    if dir_sim_now.count == 0
		answer = struct;
		IsCancelled = 1;
        fprintf(1, "Can not find any simulation results.\n" + ...
            "Please execute sim_estimate_channel.m first!");
		return;
    end

    % %%%%%%%%%%%%%%%%%%%%
    % combine two source and merge them (also delete invalid folder)
    % %%%%%%%%%%%%%%%%%%%%
    dir_now = struct;

    tmp_names = intersect(dir_gen_now.name, dir_sim_now.name);
    tmp_name_idx = zeros(length(tmp_names), 1);
    for i_name = 1 : length(tmp_names)
        name = tmp_names(i_name);
        tmp_name_idx(i_name) = find(strcmp(dir_sim_now.name, name));
    end
    tmp_name_idx(tmp_name_idx == 0) = [];
    tmp_name_idx = sort(tmp_name_idx);
    
    count = 0;
    name = {};
    array = [];
    for i_arr = 1 : length(tmp_name_idx)
        exp_struct = dir_sim_now.array(tmp_name_idx(i_arr));
        if exp_struct.count ~= 0
            count = count + 1;
            name{end+1} = dir_sim_now.name{tmp_name_idx(i_arr)}; %#ok<AGROW> 
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
    defInd.i_sim_src = 1;
    defInd.i_exp_src = 1;

	defAns = struct;
	defAns.Mc = default_Mc;
	defAns.exp_name = default_exp_name;
	defAns.flag_collect_mode = false;
	defAns.flag_skip_simulated = true;
	defAns.flag_use_parfor = true;
	defAns.flag_use_gpu = false;

    % %%%%%%%%%%%%%%%%%%%%
    % get MD5 of dir_now
    % %%%%%%%%%%%%%%%%%%%%
    md5 = GetMD5(dir_now, 'Array');

	defVal_path = fullfile('_tmp', 'defaultInput_sim_ber.mat');
    if isfile(defVal_path)
		defVal = load(defVal_path).defVal;

        md5_prev = defVal.md5;

        if strcmp(md5, md5_prev)
            defInd = defVal.defInd;
			defAns = defVal.defAns;
		else
			defAns.exp_name = defVal.defAns.exp_name;
	        defAns.flag_collect_mode = defVal.defAns.flag_collect_mode;
	        defAns.flag_skip_simulated = defVal.defAns.flag_skip_simulated;
			defAns.flag_use_parfor = defVal.defAns.flag_use_parfor;
			defAns.flag_use_gpu = defVal.defAns.flag_use_gpu;
        end
    end

    % %%%%%%%%%%%%%%%%%%%%
    % default answer 2
    % %%%%%%%%%%%%%%%%%%%%
	defAns.sim_src_name = dir_now.name(defInd.i_sim_src);
	defAns.exp_src_name = dir_now.array(defInd.i_sim_src).name(defInd.i_exp_src);
	
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%% ask for user input %%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	prompt = {};
	title = 'Simulation setting';
	formats = {};
	options = struct;
	
	prompt(1,:) = {'Enter simulation data source:', 'sim_src_name', []};
	formats(1,1).type = 'list';
	formats(1,1).style = 'popupmenu';
	formats(1,1).size = -1;
	formats(1,1).format = 'text';
	formats(1,1).items = dir_now.name;
	formats(1,1).callback = @(~,~,h,k)(set(h(k+1), ...
		'String', dir_now.array(get(h(k),'Value')).name, ...
		'Value', 1));
	formats(1,1).span = [1 4];
	
	prompt(2,:) = {'Enter experiment source:', 'exp_src_name',[]};
	formats(2,1).type = 'list';
	formats(2,1).style = 'popupmenu';
	formats(2,1).size = -1;
	formats(2,1).format = 'text';
	formats(2,1).items = dir_now.array(defInd.i_sim_src).name;
	formats(2,1).span = [1 4];

	prompt(3,:) = {'number of simulation trials', 'Mc',[]};
	formats(3,1).type = 'edit';
	formats(3,1).size = -1;
	formats(3,1).format = 'integer';
	formats(3,1).span = [1 4];

	prompt(4,:) = {'experiment', 'exp_name',[]};
	formats(4,1).type = 'edit';
	formats(4,1).size = -1;
	formats(4,1).format = 'text';
	formats(4,1).span = [1 4];

	prompt(5,:) = {'collect mode', 'flag_collect_mode',[]};
	formats(5,1).type = 'check';

	prompt(6,:) = {'skip simulated', 'flag_skip_simulated',[]};
	formats(5,2).type = 'check';

	prompt(7,:) = {'use parfor', 'flag_use_parfor',[]};
	formats(5,3).type = 'check';

	prompt(8,:) = {'use GPU', 'flag_use_gpu',[]};
	formats(5,4).type = 'check';

	options.Resize = 'off';
	options.AlignControls = 'on';
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
	if IsCancelled == 0
		answer.sim_src_name = answer.sim_src_name{1};
		answer.exp_src_name = answer.exp_src_name{1};

        % record defAns
        defVal = struct;
        defVal.md5 = md5;
        % compute ind
        ind = struct;
        ind.i_sim_src = find(ismember(dir_now.name, answer.sim_src_name));
        ind.i_exp_src = find(ismember(dir_now.array(ind.i_sim_src).name, answer.exp_src_name));
        % save to defVal
        defVal.defInd = ind;
        defVal.defAns = answer;
        if ~exist('_tmp', 'dir');	mkdir('_tmp'); end
        save(defVal_path, 'defVal');
	end
end
