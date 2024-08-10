%
% Generate simulation trials for BER
%
clc;
clear;
close;

%% flags & variables
SHOW_PROGRESS = 1;

SEED = 506;

glob_flag = struct;	% global flags
algo_flag = struct;	% algorithm flags

% %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
glob_flag.add_noise             = 1;
glob_flag.enable_CP_sim         = 0;
glob_flag.assume_known_channel  = 0;
glob_flag.use_correct_dfo       = 0;
glob_flag.use_correct_aoa       = 0;
glob_flag.use_correct_chl       = 0;

% %%%%%%%%%%%%%%%%%%%% flags for BER %%%%%%%%%%%%%%%%%%%%

% ==================== equalizer flags ====================
% [equalizer]: used to estimate data
%   ZF
%   MMSE
% [assume_no_dfo]: consider DFO effect
% [assume_no_aoa]: consider AOA effect
algo_flag.equalizer     = 'MMSE';
algo_flag.assume_no_dfo = 0;
algo_flag.assume_no_aoa = 0;

% ==================== post-processing flags ====================
% post-processing used to correct some error after equalizer
% [post]
%   NONE (use hard decision)
%   ML_SEARCH
algo_flag.post = 'NONE';

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

% -------------------- [post-processing] options --------------------
if strcmpi(algo_flag.post, 'NONE')
    % no options

elseif strcmpi(algo_flag.post, 'ML_SEARCH')
    % the threshold to determine error symbol
    algo_flag.post_opt.thres = 0.01;
    % maximal error symbol would be considerd
    algo_flag.post_opt.max_error_symbol_number = 5;

else
    error("wrong post-processing method");
end

if glob_flag.assume_known_channel == 1
	glob_flag.use_correct_dfo = 1;
	glob_flag.use_correct_aoa = 1;
	glob_flag.use_correct_chl = 1;
end

% convert to one flag:
glob_flag.gen_chl_mtx_flags = uint8(0);
if (glob_flag.use_correct_dfo == 1)
    glob_flag.gen_chl_mtx_flags = bitset(glob_flag.gen_chl_mtx_flags, 1, 1);
end
if (glob_flag.use_correct_aoa == 1)
    glob_flag.gen_chl_mtx_flags = bitset(glob_flag.gen_chl_mtx_flags, 2, 1);
end
if (glob_flag.use_correct_chl == 1)
    glob_flag.gen_chl_mtx_flags = bitset(glob_flag.gen_chl_mtx_flags, 3, 1);
end

% check DEBUG:
if glob_flag.add_noise == 0;warning("noise free");end
if glob_flag.assume_known_channel == 1;warning("ideal case simulation");end
if glob_flag.use_correct_dfo == 1;warning("use correct DFO");end
if glob_flag.use_correct_aoa == 1;warning("use correct AoA");end
if glob_flag.use_correct_chl == 1;warning("use correct channel");end
if algo_flag.assume_no_dfo == 1;warning("receiver ignores DFO");end
if algo_flag.assume_no_aoa == 1;warning("receiver ignores AoA");end

%% add path
addpath(genpath('src'));

%% simulation settings
% no. of Monte Carlo
default_sim_Mc  = 5000;
% range of snr
default_sim_snr = -10 : 10 : 40;
% no. of path
default_sim_Q   = 1 : 5;
% # of Rx antenna
default_sim_Nr  = 1;
% default experiment name
default_exp_name = 'first try';

%% user inputl
[answer, IsCancelled] = request_user_input_sim_ber(default_sim_Mc, default_exp_name);

if IsCancelled;	return;	end

% 1. load flag
glob_flag.use_gpu           = answer.flag_use_gpu;
glob_flag.use_parfor        = answer.flag_use_parfor;
% skip simulated i_Mc runs
glob_flag.skip_simulated    = answer.flag_skip_simulated;
% only collect, no simulation
glob_flag.collect_mode      = answer.flag_collect_mode;

% 2. get infomation
sim_name        = answer.sim_src_name;
exp_name        = answer.exp_src_name;
default_sim_Mc  = answer.Mc;
sub_trial       = answer.exp_name;
dir_exp         = fullfile(sim_name, exp_name);
dir_sim         = fullfile(sim_name, sub_trial);

%% load data
% get config
config_name     = sprintf('%s_config.mat', sim_name);
config          = load(fullfile('_gendata', sim_name, config_name)).config;

gen_config_name = sprintf('%s_gen_config.mat', sim_name);
gen_config      = load(fullfile('_gendata', sim_name, gen_config_name));

est_config_name = sprintf('%s_sim_config.est.mat', sim_name);
est_config      = load(fullfile('_sim_result_est', dir_exp, est_config_name));

%   get Mc
Mc = min([default_sim_Mc, gen_config.Mc, est_config.Mc]);

%   get snr
[snr, unmatched_snr] = intersects(default_sim_snr, gen_config.snr, est_config.snr);
if unmatched_snr;	warning("different simulation snr setting");	end
snr_gendata_idx = arrcmp(snr, gen_config.snr, 1);
snr_estchal_idx = arrcmp(snr, est_config.snr, 1);

%   get Q
[Q, unmatched_Q] = intersects(default_sim_Q, gen_config.Q, est_config.Q);
if unmatched_Q;	warning("different simulation Q setting");	end
Q_gendata_idx = arrcmp(Q, gen_config.Q, 1);
Q_estchal_idx = arrcmp(Q, est_config.Q, 1);

%   get Nr
Nr = default_sim_Nr;
if Nr ~= est_config.Nr
    error("different simulation Nr setting");
end

%% save directory
dir_dump = fullfile('_dump', dir_sim);
if ~exist(dir_dump, 'dir');	mkdir(dir_dump); end

%% auxiliary variables
n_snr   = length(snr);
n_Q     = length(Q);

%% pack variables
sim_Mc_info = struct;
% simulation info
sim_Mc_info.snr     = snr;
sim_Mc_info.Q       = Q;
sim_Mc_info.Nr      = Nr;
sim_Mc_info.config  = config;
% simulation flags
sim_Mc_info.glob_flag = glob_flag;
sim_Mc_info.algo_flag = algo_flag;
% directory of data
sim_Mc_info.dir_gendata = fullfile('_gendata', sim_name, sim_name);
sim_Mc_info.dir_estdata = fullfile('_sim_result_est', dir_exp, sim_name);
% indexes mapping
sim_Mc_info.snr_gendata_idx = snr_gendata_idx;
sim_Mc_info.snr_estchal_idx = snr_estchal_idx;
sim_Mc_info.Q_gendata_idx   = Q_gendata_idx;
sim_Mc_info.Q_estchal_idx   = Q_estchal_idx;

%% declare record var.
r_avg_ber_Mc = zeros(n_Q, n_snr, Mc);

%% run simulation
if glob_flag.use_parfor ~= 0
	poolobj     = startparpool();
	parforArg   = poolobj.NumWorkers;
else
	parforArg   = 0;
end

% set random seed: rng(SEED); for each worker
% otherwise, the random number generated for each worker may not be same
prngsc  = parallel.pool.Constant(RandStream('Threefry', 'Seed', SEED));
psimset = parallel.pool.Constant(sim_Mc_info);

% determine the sim_mode: collect, skip simulated
sim_mode = 2 * glob_flag.collect_mode + glob_flag.skip_simulated;
valid_Mc = ones(Mc, 1);

% start simulation
fprintf(1, "Simulation info:\n");
fprintf(1, "\tsim_src:\t%s\n", sim_name);
fprintf(1, "\texp_src:\t%s\n", exp_name);
fprintf(1, "\texp_name:\t%s\n", sub_trial);
fprintf(1, "\tMc\t=\t%d\n", Mc);
fprintf(1, "\tQ\t=\t[%s]\n", join(string(Q), ','));
fprintf(1, "\tsnr\t=\t[%s]\n", join(string(snr), ','));
fprintf(1, "\tNr\t=\t%d\n", Nr);
fprintf(1, "Start simulation at %s\n", string(datetime('now')));
fprintf(1, "--------------------------------------------------\n");

% if SHOW_PROGRESS == 1;  parfor_progress(Mc);    end
if SHOW_PROGRESS == 1;  PB = ProgressBar_cli_parfor(Mc);    end
tic;
% disp("============!!!!!!!!!!!!!!!!!!!!!============")
%parfor (i_Mc = 1 : Mc, parforArg)
 for (i_Mc = 1 : Mc)
% for (i_Mc = 26)
    % set random seed
    stream              = prngsc.Value;
    stream.Substream    = i_Mc;
    RandStream.setGlobalStream(stream);

    % path name for i_Mc
    path_full_name  = sprintf('%s_avg_ber_%03d.ber.mat', sim_name, i_Mc);
	path_full       = fullfile(dir_dump, path_full_name);
	
	% check mode
	if sim_mode >= 2
		if isfile(path_full)
			% load file
            data        = load(path_full);
			r_avg_ber   = data.r_avg_ber;

            %	get snr & Q
            [load_Q, unmatched_Q]       = intersects(Q, data.Q);
            [load_snr, unmatched_snr]   = intersects(snr, data.snr);
        
            % record
            if unmatched_snr || unmatched_Q
                load_Q_idx      = arrcmp(load_Q, data.Q, 1);
                load_snr_idx    = arrcmp(load_snr, data.snr, 1);
                r_avg_ber_Mc(:, :, i_Mc) = r_avg_ber(load_Q_idx, load_snr_idx);
            else
                r_avg_ber_Mc(:, :, i_Mc) = r_avg_ber;
            end			
		else
			valid_Mc(i_Mc) = 0;
		end

    elseif sim_mode >= 1 && isfile(path_full)
		% load file
        data        = load(path_full);
		r_avg_ber   = data.r_avg_ber;

        %	get snr & Q
        [load_snr, unmatched_snr]   = intersects(snr, data.snr);
        [load_Q, unmatched_Q]       = intersects(Q, data.Q);

		% record
        if unmatched_snr || unmatched_Q
            load_Q_idx      = arrcmp(load_Q, data.Q, 1);
            load_snr_idx    = arrcmp(load_snr, data.snr, 1);            
            r_avg_ber_Mc(:, :, i_Mc) = r_avg_ber(load_Q_idx, load_snr_idx);
        else
            r_avg_ber_Mc(:, :, i_Mc) = r_avg_ber;
        end	
        
	else
		% run simulation
		r_avg_ber = one_Mc_ber(i_Mc, psimset.Value);

		% record
		r_avg_ber_Mc(:, :, i_Mc) = r_avg_ber;

		% dump record
		parsave(path_full, r_avg_ber, Q, snr, Nr);
	end
	
% 	if SHOW_PROGRESS == 1;   parfor_progress;    end
	if SHOW_PROGRESS == 1;   advance(PB);    end
end

% if SHOW_PROGRESS == 1;  parfor_progress(0); end
fprintf(1, "done!\n");

%% show info
sim_exe_time = toc;
fprintf(1, "\nexecution time: %s\n", show_duration(sim_exe_time));

%% average over Mc
% note: average over r_avg_ber_Mc automatically skip miss_Mc since it
% initialized with all zeros
avg_ber = squeeze(sum(r_avg_ber_Mc, 3) ./ nnz(valid_Mc));
Mc      = nnz(valid_Mc);

%% save files
dir_result = fullfile('results', dir_sim);
if ~exist(dir_result, 'dir');	mkdir(dir_result); end

path_full_name  = sprintf('%s_channel_config.ber.mat', sim_name);
path_full       = fullfile(dir_result, path_full_name);
save(path_full, 'config');

path_full_name  = sprintf('%s_sim_config.ber.mat', sim_name);
path_full       = fullfile(dir_result, path_full_name);
save(path_full, 'Mc', 'Q', 'snr', 'Nr', 'sim_exe_time');

path_full_name  = sprintf('%s_avg_ber.ber.mat', sim_name);
path_full       = fullfile(dir_result, path_full_name);
save(path_full, 'avg_ber');
