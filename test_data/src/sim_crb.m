%
% Compute CRB
%
clc;
clear;
close;

%% flags & variables
SHOW_PROGRESS   = 1;

SEED = 506;

%% add path
addpath(genpath('src'));

%% simulation settings
% no. of Monte Carlo
default_sim_Mc  = 5000;
% range of snr
default_sim_snr = -10 : 10 : 40;
% no. of path
default_sim_Q   = 1 : 4;
% no. of subvector (periodic training sequence)
default_sim_P   = 4;
% # of Rx antenna
default_sim_Nr  = 1;

% default experiment name
default_exp_name = 'first try';

%% user input
[answer, IsCancelled] = request_user_input_sim_crb(default_sim_Mc, default_exp_name);

if IsCancelled;	return;	end

sim_name        = answer.sim_src_name;
default_sim_Mc  = answer.Mc;
sub_trial       = answer.exp_name;
dir_sim         = fullfile(sim_name, sub_trial);
flag_use_parfor = answer.flag_use_parfor;

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

%   get P
P = default_sim_P;

%   get Nr
Nr = default_sim_Nr;

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
sim_Mc_info.dir_gendata = fullfile('_gendata', sim_name, sim_name);
sim_Mc_info.snr_idx     = snr_idx;
sim_Mc_info.Q_idx       = Q_idx;

%% declare record var.
r_dfo_crb_Mc = zeros(n_Q, n_snr, Mc);
r_aoa_crb_Mc = zeros(n_Q, n_snr, Mc);

%% run simulation
if flag_use_parfor ~= 0
	poolobj = startparpool();
	parforArg = poolobj.NumWorkers;
else
	parforArg = 0;
end

parfevalOnAll(@warning,0,'off','all')

% set random seed: rng(SEED); for each worker
% otherwise, the random number generated for each worker may not be same
prngsc  = parallel.pool.Constant(RandStream('Threefry', 'Seed', SEED));
psimset = parallel.pool.Constant(sim_Mc_info);

% start simulation
fprintf(1, "Simulation info:\n");
fprintf(1, "\tsim_src:\t%s\n", sim_name);
fprintf(1, "\texp_name:\t%s\n", sub_trial);
fprintf(1, "\tMc\t=\t%d\n", Mc);
fprintf(1, "\tQ\t=\t[%s]\n", join(string(Q), ','));
fprintf(1, "\tsnr\t=\t[%s]\n", join(string(snr), ','));
fprintf(1, "\tP\t=\t%d\n", P);
fprintf(1, "\tNr\t=\t%d\n", Nr);
fprintf(1, "Start simulation at %s\n", string(datetime('now')));
fprintf(1, "--------------------------------------------------\n");

% if SHOW_PROGRESS == 1;  parfor_progress(Mc);    end
if SHOW_PROGRESS == 1;  PB = ProgressBar_cli_parfor(Mc);    end
tic;

%parfor (i_Mc = 1 : Mc, parforArg)
 for i_Mc = 1 : Mc
    [r_dfo_crb, r_aoa_crb] = one_Mc_crb(i_Mc, psimset.Value);
    r_dfo_crb_Mc(:, :, i_Mc) = r_dfo_crb;
    r_aoa_crb_Mc(:, :, i_Mc) = r_aoa_crb;

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
avg_dfo_crb = squeeze(sum(r_dfo_crb_Mc, 3) ./ Mc);
avg_aoa_crb = squeeze(sum(r_aoa_crb_Mc, 3) ./ Mc);

%% save files
dir_result = fullfile('results_crb', dir_sim);
if ~exist(dir_result, 'dir');	mkdir(dir_result); end

path_full_name  = sprintf('%s_channel_config.crb.mat', sim_name);
path_full       = fullfile(dir_result, path_full_name);
save(path_full, 'config');

path_full_name  = sprintf('%s_sim_config.crb.mat', sim_name);
path_full       = fullfile(dir_result, path_full_name);
save(path_full, 'Mc', 'Q', 'snr', 'P', 'Nr', 'sim_exe_time');

path_full_name  = sprintf('%s_avg_crb.crb.mat', sim_name);
path_full       = fullfile(dir_result, path_full_name);
save(path_full, 'avg_dfo_crb', 'avg_aoa_crb');
