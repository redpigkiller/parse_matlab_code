%
% Generate simulation trials for estimating channel DFO
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
glob_flag.add_noise     = 1;
glob_flag.enable_CP_sim = 0;

% %%%%%%%%%%%%%%%%%%%% flags for estimating channel %%%%%%%%%%%%%%%%%%%%
% ==================== simulation flags ====================
% set to 0, if you do not want to estimate DFO, and directly estimate
% the channel
algo_flag.estimate_dfo = 1;

% set to 0, if you do not want to estimate channel
algo_flag.estimate_chl = 0;

% ==================== algorithm flags ====================
% algorithm used to estimate DFOs & AOAs:
%   ESPRIT
%   UNITARY_ESPRIT
%   ML_Estimation
%   UNITARY_ESPRIT_2D
%
% enable_sep: separate too close angles
%
% algorithm used to estimate Channel:
%   stage 1:
%       PINV, TIKHONOV, PINV_KRON, TIKHONOV_KRON
%   stage 2:
%       PINV, MMSE, [CS]
algo_flag.use_dfo_aoa_algo      = 'ML_Estimation';
algo_flag.enable_sep            = true;
algo_flag.use_chl_algo_stage_1  = 'TIKHONOV';
algo_flag.use_chl_algo_stage_2  = 'PINV';

% -------------------- [DFO] options --------------------
if strcmpi(algo_flag.use_dfo_aoa_algo, 'ESPRIT')
	% the following option are activated only when estimate_dfo set to 1
	% minus_mean: subtract mean of estimated vector before ESPRIT
	algo_flag.dfo_opt.minus_mean = 0;
    % sig_subspace:
    %   SVD: the direct data approach
    %   EVD: the covariance approach
	% method:
	%   LS: (Least-Squares Solution)
	%   TLS: (Total-Least-Squares Optimization)
    algo_flag.dfo_opt.sig_subspace  = 'SVD';
	algo_flag.dfo_opt.method        = 'TLS';

elseif strcmpi(algo_flag.use_dfo_aoa_algo, 'UNITARY_ESPRIT')
	% the following option are activated only when estimate_dfo set to 1
    % sig_subspace:
    %   SVD: square root approach
    %   EVD: covariance approach
    % method:
	%   LS: (Least-Squares Solution)
	%   TLS: (Total-Least-Squares Optimization)
    algo_flag.dfo_opt.sig_subspace  = 'SVD';
    algo_flag.dfo_opt.method        = 'TLS';

elseif strcmpi(algo_flag.use_dfo_aoa_algo, 'ML_Estimation')
	% [TL;DR] please refers to the help of the ML_Estimation (press F1)
    % algorithms worth to try:
    %   1. PSO with swarn size set to 100
    %   2. PSO (swarn size 20) + HybridFcn (fmincon)
    %   3. fmincon
    %
    % if the MSE of channel is too high, consider set
    % algo_flag.use_chl_algo_stage_1 to 'TIKHONOV'
	algo_flag.dfo_opt.algo = 'PSO';
% 	algo_flag.dfo_opt.algo = 'GA';
    % algoOpt:    options for optimization algorithm
    %   use '' for default option
    %   options states for b(bound), o(order)
	algo_flag.dfo_opt.algoOpt   = @(Q) iif(Q==1, 'b', 'b');
	% method: choose the methods to estimation DFO
	%  1: estimate eps_q * M
	%  2: directly estimate eps_q
	algo_flag.dfo_opt.method    = @(Q) iif(Q==1, 1, 1);

elseif strcmpi(algo_flag.use_dfo_aoa_algo, 'UNITARY_ESPRIT_2D')
	% the following option are activated only when estimate_dfo set to 1
    % sig_subspace:
    %   SVD: square root approach
    %   EVD: covariance approach
    % method:
	%   LS: (Least-Squares Solution)
	%   TLS: (Total-Least-Squares Optimization)
    algo_flag.dfo_aoa_opt.sig_subspace  = 'SVD';
    algo_flag.dfo_aoa_opt.method        = 'TLS';

else
    error("wrong DFO algorithm");
end

% -------------------- [Channel] options for stage 1 --------------------
if strcmpi(algo_flag.use_chl_algo_stage_1, 'PINV')
	% no options

elseif strcmpi(algo_flag.use_chl_algo_stage_1, 'TIKHONOV')
    algo_flag.chl_opt_stage_1.reg = 0.001;

elseif strcmpi(algo_flag.use_chl_algo_stage_1, 'PINV_KRON')
    % no options

elseif strcmpi(algo_flag.use_chl_algo_stage_1, 'TIKHONOV_KRON')
    algo_flag.chl_opt_stage_1.reg = 0.001;

else
    error("wrong Channel stage 1 algorithm");
end

% -------------------- [Channel] options for stage 2 --------------------
if strcmpi(algo_flag.use_chl_algo_stage_2, 'PINV')
	% no options

elseif strcmpi(algo_flag.use_chl_algo_stage_2, 'MMSE')
    algo_flag.chl_opt_stage_2.chlpwr = 0.3;

elseif strcmpi(algo_flag.use_chl_algo_stage_2, 'CS')
    algo_flag.chl_opt_stage_2.chllen = 5;
else
    error("wrong Channel stage 2 algorithm");
end

% check DEBUG:
if glob_flag.add_noise == 0;warning("noise free");end
if algo_flag.estimate_dfo == 0;warning("no DFO estimation");end
if algo_flag.estimate_chl == 0;warning("no channel estimation");end

%% add path
addpath(genpath('src'));            

%% simulation settings
% no. of Monte Carlo
default_sim_Mc  = 5000;
% range of snr
default_sim_snr = -10 : 10 : 40;
% no. of path
default_sim_Q   = 1 : 5;
% no. of subvector (periodic training sequence)
default_sim_P   = 8;
% # of Rx antenna
default_sim_Nr  = 1;

% the multiple of the N
N_ratio = 8;

% default experiment name
default_exp_name = 'first try';

%% user input
[answer, IsCancelled] = request_user_input_sim_estimate(default_sim_Mc, default_exp_name);

if IsCancelled;	return;	end

% 1. load flag
glob_flag.use_parfor    = answer.flag_use_parfor;
glob_flag.use_gpu       = answer.flag_use_gpu;
if strcmpi(algo_flag.use_dfo_aoa_algo, 'ML_Estimation') && glob_flag.use_gpu
	glob_flag.use_gpu = 0;
	warning("ML_Estimation does NOT support GPU! Continue using CPU...");
end
% skip simulated i_Mc runs
glob_flag.skip_simulated    = answer.flag_skip_simulated;
% only collect, no simulation
glob_flag.collect_mode      = answer.flag_collect_mode;

% 2. get infomation
sim_name        = answer.sim_src_name;
default_sim_Mc  = answer.Mc;
sub_trial       = answer.exp_name;
dir_sim         = fullfile(sim_name, sub_trial);

%% load data
% get config
config_name     = sprintf('%s_config.mat', sim_name);
config          = load(fullfile('_gendata', sim_name, config_name)).config;

gen_config_name = sprintf('%s_gen_config.mat', sim_name);
gen_config      = load(fullfile('_gendata', sim_name, gen_config_name));

% unpack variables
sim_Mc  = gen_config.Mc;
sim_snr = gen_config.snr;
sim_Q   = gen_config.Q;

%   get Mc
Mc = min([default_sim_Mc, sim_Mc]);

%   get snr
[snr, unmatched_snr] = intersects(default_sim_snr, sim_snr);
if unmatched_snr;	warning("different simulation snr setting");	end
snr_idx = arrcmp(snr, sim_snr, 1);

%   get Q
[Q, unmatched_Q] = intersects(default_sim_Q, sim_Q);
if unmatched_Q;	warning("different simulation Q setting");	end
Q_idx = arrcmp(Q, sim_Q, 1);

%   get Nr
Nr = default_sim_Nr;
if Nr > 1 && ~strcmpi(algo_flag.use_dfo_aoa_algo, 'UNITARY_ESPRIT_2D')
    error("Can not use 1D algorithm to estimate a 2D problem!");
end

%   get P
P = default_sim_P;

%% save directory
dir_dump = fullfile('_dump', dir_sim);
if ~exist(dir_dump, 'dir');	mkdir(dir_dump); end

dir_save = fullfile('_sim_result_est', dir_sim);
if ~exist(dir_save, 'dir');	mkdir(dir_save); end

%% auxiliary variables
n_snr   = length(snr);
n_Q     = length(Q);

%% pack variables
sim_Mc_info = struct;
sim_Mc_info.snr         = snr;
sim_Mc_info.Q           = Q;
sim_Mc_info.P           = P;
sim_Mc_info.Nr          = Nr;

sim_Mc_info.config      = config;
if N_ratio ~= 1
    sim_Mc_info.config.N = sim_Mc_info.config.N * N_ratio;
    warning("N_ratio != 1");
end

sim_Mc_info.glob_flag   = glob_flag;
sim_Mc_info.algo_flag   = algo_flag;
sim_Mc_info.dir_gendata = fullfile('_gendata', sim_name, sim_name);
sim_Mc_info.snr_idx     = snr_idx;
sim_Mc_info.Q_idx       = Q_idx;

%% declare record var.
r_MSE_dfo_Mc = zeros(n_Q, n_snr, Mc);
r_MSE_aoa_Mc = zeros(n_Q, n_snr, Mc);
r_MSE_chl_Mc = zeros(n_Q, n_snr, Mc);
s_est_dfo_Mc = cell(n_Q, n_snr, Mc);
s_est_aoa_Mc = cell(n_Q, n_snr, Mc);
s_est_chl_Mc = cell(n_Q, n_snr, Mc);

%% run simulation
if glob_flag.use_parfor ~= 0
	poolobj = startparpool();
	parforArg = poolobj.NumWorkers;
else
	parforArg = 0;
end

% set random seed: rng(SEED); for each worker
% otherwise, the random number generated for each worker may not be same
prngsc  = parallel.pool.Constant(RandStream('Threefry', 'Seed', SEED));
psimset = parallel.pool.Constant(sim_Mc_info);

% determine the sim_mode: collect, skip simulated
sim_mode = 2 * glob_flag.collect_mode + glob_flag.skip_simulated;
valid_Q_Mc      = zeros(n_Q, Mc);
valid_snr_Mc    = zeros(n_snr, Mc);

% start simulation
fprintf(1, "Simulation info:\n");
fprintf(1, "\tsim_src:\t%s\n", sim_name);
fprintf(1, "\texp_name:\t%s\n", sub_trial);
fprintf(1, "\tMc\t=\t%d\n", Mc);
fprintf(1, "\tQ\t=\t[%s]\n", join(string(Q), ','));
fprintf(1, "\tsnr\t=\t[%s]\n", join(string(snr), ','));
fprintf(1, "\tP\t=\t%d\n", P);
fprintf(1, "\tNr\t=\t%d\n", Nr);
fprintf(1, "\tN_ratio\t=\t%d\n", N_ratio);
fprintf(1, "Start simulation at %s\n", string(datetime('now')));
fprintf(1, "--------------------------------------------------\n");

% if SHOW_PROGRESS == 1;  parfor_progress(Mc);    end
if SHOW_PROGRESS == 1;  PB = ProgressBar_cli_parfor(Mc);    end
tic;

%parfor (i_Mc = 1 : Mc, parforArg)
 for (i_Mc = 1 : Mc)
    % set random seed
    stream              = prngsc.Value;
    stream.Substream    = i_Mc;
    RandStream.setGlobalStream(stream);

    % path name for i_Mc
    path_full_name  = sprintf('%s_mse_%03d.est.mat', sim_name, i_Mc);
	r_path_full     = fullfile(dir_dump, path_full_name);
    
    % save estimated data
    path_full_name  = sprintf('%s_est_%03d.est.mat', sim_name, i_Mc);
    s_path_full     = fullfile(dir_save, path_full_name);

    % check mode
    if sim_mode >= 2
        if isfile(r_path_full)
            % load file
            r_data = load(r_path_full);
		    est_MSE_dfo = r_data.est_MSE_dfo;
		    est_MSE_aoa = r_data.est_MSE_aoa;
		    est_MSE_chl = r_data.est_MSE_chl;

            s_data = load(s_path_full);
            est_dfoVal = s_data.est_dfoVal;
            est_aoaVal = s_data.est_aoaVal;
            est_chlVal = s_data.est_chlVal;

            %   get snr & Q
            [l_Q, unmatched_Q]       = intersects(Q, r_data.Q);
            [l_snr, unmatched_snr]   = intersects(snr, r_data.snr);

            f_Q_idx      = arrcmp(l_Q, Q, 1);
            f_snr_idx    = arrcmp(l_snr, snr, 1);

		    % record
            if unmatched_snr || unmatched_Q
                l_Q_idx      = arrcmp(l_Q, r_data.Q, 1);
                l_snr_idx    = arrcmp(l_snr, r_data.snr, 1);

                % load DFO
                fill_r_data = zeros(n_Q, n_snr);
                fill_s_data = cell(n_Q, n_snr);
                fill_r_data(f_Q_idx, f_snr_idx) = est_MSE_dfo(l_Q_idx, l_snr_idx);
                fill_s_data(f_Q_idx, f_snr_idx) = est_dfoVal(l_Q_idx, l_snr_idx);
                r_MSE_dfo_Mc(:, :, i_Mc) = fill_r_data;
                s_est_dfo_Mc(:, :, i_Mc) = fill_s_data;

                % load AOA
                fill_r_data = zeros(n_Q, n_snr);
                fill_s_data = cell(n_Q, n_snr);
                fill_r_data(f_Q_idx, f_snr_idx) = est_MSE_aoa(l_Q_idx, l_snr_idx);
                fill_s_data(f_Q_idx, f_snr_idx) = est_aoaVal(l_Q_idx, l_snr_idx);
                r_MSE_aoa_Mc(:, :, i_Mc) = fill_r_data;
                s_est_aoa_Mc(:, :, i_Mc) = fill_s_data;

                % load CHL
                fill_r_data = zeros(n_Q, n_snr);
                fill_s_data = cell(n_Q, n_snr);
                fill_r_data(f_Q_idx, f_snr_idx) = est_MSE_chl(l_Q_idx, l_snr_idx);
                fill_s_data(f_Q_idx, f_snr_idx) = est_chlVal(l_Q_idx, l_snr_idx);
                r_MSE_chl_Mc(:, :, i_Mc) = fill_r_data;
                s_est_chl_Mc(:, :, i_Mc) = fill_s_data;

            else
                r_MSE_dfo_Mc(:, :, i_Mc) = est_MSE_dfo;
	            r_MSE_aoa_Mc(:, :, i_Mc) = est_MSE_aoa;
	            r_MSE_chl_Mc(:, :, i_Mc) = est_MSE_chl;
	            s_est_dfo_Mc(:, :, i_Mc) = est_dfoVal;
	            s_est_aoa_Mc(:, :, i_Mc) = est_aoaVal;
	            s_est_chl_Mc(:, :, i_Mc) = est_chlVal;
            end

            % add valid Mc
            % Convert the indices to a one-hot vector
            valid_Q_Mc(:, i_Mc)     = full(sparse(f_Q_idx, ones(length(f_Q_idx), 1), 1, n_Q, 1));
            valid_snr_Mc(:, i_Mc)   = full(sparse(f_snr_idx, ones(length(f_snr_idx), 1), 1, n_snr, 1));

        end

    elseif sim_mode >= 1 && isfile(r_path_full) && isfile(s_path_full)
        % load file
        r_data = load(r_path_full);
        s_data = load(s_path_full);
        
        %	get snr & Q
        [l_Q, r_unmatched_Q]        = intersects(Q, r_data.Q);
        [l_snr, r_unmatched_snr]    = intersects(snr, r_data.snr);
        [l_Q, s_unmatched_Q]        = intersects(l_Q, s_data.Q);
        [l_snr, s_unmatched_snr]    = intersects(l_snr, s_data.snr);

		% record
        if r_unmatched_snr || s_unmatched_snr || r_unmatched_Q || s_unmatched_Q
            % first, load the existed simulated data
            l_Q_idx      = arrcmp(l_Q, r_data.Q, 1);
            l_snr_idx    = arrcmp(l_snr, r_data.snr, 1);

            % second, run simulation for miss indices
		    [est_MSE_dfo, est_MSE_aoa, est_MSE_chl, est_dfoVal, est_aoaVal, est_chlVal] = ...
                one_Mc_est(i_Mc, psimset.Value, l_Q, l_snr);

            % record MSE
            f_Q_idx      = arrcmp(l_Q, Q, 1);
            f_snr_idx    = arrcmp(l_snr, snr, 1);

            est_MSE_dfo(f_Q_idx, f_snr_idx) = r_data.est_MSE_dfo(l_Q_idx, l_snr_idx);
            est_MSE_aoa(f_Q_idx, f_snr_idx) = r_data.est_MSE_aoa(l_Q_idx, l_snr_idx);
            est_MSE_chl(f_Q_idx, f_snr_idx) = r_data.est_MSE_chl(l_Q_idx, l_snr_idx);
        
            % record estimated data
            s_data = load(s_path_full);
            est_dfoVal(f_Q_idx, f_snr_idx) = s_data.est_dfoVal(l_Q_idx, l_snr_idx);
            est_aoaVal(f_Q_idx, f_snr_idx) = s_data.est_aoaVal(l_Q_idx, l_snr_idx);
            est_chlVal(f_Q_idx, f_snr_idx) = s_data.est_chlVal(l_Q_idx, l_snr_idx);

            % record
            r_MSE_dfo_Mc(:, :, i_Mc) = est_MSE_dfo;
            r_MSE_aoa_Mc(:, :, i_Mc) = est_MSE_aoa;
            r_MSE_chl_Mc(:, :, i_Mc) = est_MSE_chl;
	        s_est_dfo_Mc(:, :, i_Mc) = est_dfoVal;
	        s_est_aoa_Mc(:, :, i_Mc) = est_aoaVal;
	        s_est_chl_Mc(:, :, i_Mc) = est_chlVal;

            % dump save
            dump_name = sprintf('%s_est_%03d.est.mat', sim_name, i_Mc);
	        r_path_full = fullfile(dir_save, dump_name);
	        parsave(r_path_full, est_dfoVal, est_aoaVal, est_chlVal, Q, snr, P, Nr);
        
	        % dump record
            dump_name = sprintf('%s_mse_%03d.est.mat', sim_name, i_Mc);
	        r_path_full = fullfile(dir_dump, dump_name);
	        parsave(r_path_full, est_MSE_dfo, est_MSE_aoa, est_MSE_chl, Q, snr, P, Nr);

        else
            r_MSE_dfo_Mc(:, :, i_Mc) = r_data.est_MSE_dfo;
            r_MSE_aoa_Mc(:, :, i_Mc) = r_data.est_MSE_aoa;
            r_MSE_chl_Mc(:, :, i_Mc) = r_data.est_MSE_chl;
	        s_est_dfo_Mc(:, :, i_Mc) = s_data.est_dfoVal;
	        s_est_aoa_Mc(:, :, i_Mc) = s_data.est_aoaVal;
	        s_est_chl_Mc(:, :, i_Mc) = s_data.est_chlVal;

        end

        % add valid Mc
        valid_Q_Mc(:, i_Mc)     = 1;
        valid_snr_Mc(:, i_Mc)   = 1;
        
	else
		% run simulation
		[est_MSE_dfo, est_MSE_aoa, est_MSE_chl, est_dfoVal, est_aoaVal, est_chlVal] = ...
            one_Mc_est(i_Mc, psimset.Value);

		% record
	    r_MSE_dfo_Mc(:, :, i_Mc) = est_MSE_dfo;
	    r_MSE_aoa_Mc(:, :, i_Mc) = est_MSE_aoa;
	    r_MSE_chl_Mc(:, :, i_Mc) = est_MSE_chl;
	    s_est_dfo_Mc(:, :, i_Mc) = est_dfoVal;
	    s_est_aoa_Mc(:, :, i_Mc) = est_aoaVal;
	    s_est_chl_Mc(:, :, i_Mc) = est_chlVal;

        % add valid Mc
        valid_Q_Mc(:, i_Mc)     = 1;
        valid_snr_Mc(:, i_Mc)   = 1;

		% dump save
        dump_name = sprintf('%s_est_%03d.est.mat', sim_name, i_Mc);
	    r_path_full = fullfile(dir_save, dump_name);
	    parsave(r_path_full, est_dfoVal, est_aoaVal, est_chlVal, Q, snr, P, Nr);
    
	    % dump record
        dump_name = sprintf('%s_mse_%03d.est.mat', sim_name, i_Mc);
	    r_path_full = fullfile(dir_dump, dump_name);
	    parsave(r_path_full, est_MSE_dfo, est_MSE_aoa, est_MSE_chl, Q, snr, P, Nr);
    end

% 	if SHOW_PROGRESS == 1;   parfor_progress;    end
	if SHOW_PROGRESS == 1;   advance(PB);    end
end

% if SHOW_PROGRESS == 1;  parfor_progress(0); end
fprintf(1, "done!\n");

%% show info
sim_exe_time = toc;
fprintf(1, "--------------------------------------------------\n");
fprintf(1, "Simulation ends at %s\n", string(datetime('now')));
fprintf(1, "Execution time: %s\n", show_duration(sim_exe_time));
fprintf(1, "--------------------------------------------------\n");

%% average over Mc
% note: average over r_MSE_dfo_Mc & r_MSE_chl_Mc automatically skip miss_Mc
% since it initialized with all zeros
count_Mc = zeros(n_Q, n_snr);
for i_Mc = 1 : Mc
    count_Mc = count_Mc + valid_Q_Mc(:, i_Mc) * valid_snr_Mc(:, i_Mc).';
end

uni_count_Mc = unique(count_Mc);
uni_count_Mc(uni_count_Mc == 0) = [];

if length(uni_count_Mc) ~= 1
    warning("Unequal number of simulation runs!");
end

Mc = max(uni_count_Mc, [], "all");
MSE_dfo = squeeze(sum(r_MSE_dfo_Mc, 3) ./ count_Mc);
MSE_aoa = squeeze(sum(r_MSE_aoa_Mc, 3) ./ count_Mc);
MSE_chl = squeeze(sum(r_MSE_chl_Mc, 3) ./ count_Mc);

fprintf(1, "Collecting data... (valid Mc = %d)\n", Mc);

%% save files
% 1. for save:
save_name = sprintf('%s_channel_config.est.mat', sim_name);
r_path_full = fullfile(dir_save, save_name);
save(r_path_full, 'config');

save_name = sprintf('%s_sim_config.est.mat', sim_name);
r_path_full = fullfile(dir_save, save_name);
save(r_path_full, 'Mc', 'Q', 'snr', 'P', 'Nr', 'N_ratio');

% 2. for result:
dir_result = fullfile('results', dir_sim);
if ~exist(dir_result, 'dir');	mkdir(dir_result); end

result_name = sprintf('%s_channel_config.est.mat', sim_name);
r_path_full = fullfile(dir_result, result_name);
save(r_path_full, 'config');

result_name = sprintf('%s_sim_config.est.mat', sim_name);
r_path_full = fullfile(dir_result, result_name);
save(r_path_full, 'Mc', 'Q', 'snr', 'P', 'Nr', 'N_ratio', 'sim_exe_time');

result_name = sprintf('%s_mse_dfo_aoa_channel.est.mat', sim_name);
r_path_full = fullfile(dir_result, result_name);
save(r_path_full, 'MSE_dfo', 'MSE_aoa', 'MSE_chl');
