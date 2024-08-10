%
% Generate simulation trials
%
clc;
clear;
close;

%% default settings
% number of trials
default_Mc  = 5000;
% range of snr
snr         = -10 : 10 : 40;
% no. of path
Q           = 1 : 1 : 5;

%% flags
SEED = 506;

flag_use_parfor = 1;

% 0: call gen_dfo_aoa_channel
% 1: call gen_dfo_aoa_channel_simple 
flag_simple_channel = 0;

% 0: do not change channel, fix the channel all the way
% 1: change the channel with period = period_change
flag_change_ch_Mc   = 10;

% fix DFOs, AoAs, and Channel taps, but regenerate with random delay and
% power
regen_ch = 1;

%% user input
prompt      = {'Enter simulation data name:', 'Number of data'};
dlgtitle    = 'Input';
dims        = [1 35];
definput    = {'sim_01', num2str(default_Mc)};
answer      = inputdlg(prompt, dlgtitle, dims, definput);
sim_name    = answer{1};
Mc          = str2double(answer{2});

%% add path
addpath(genpath('src'));

%% load default config
run('gen_config.m');
config = load(fullfile('_config', 'default_config.mat')).config;
disp(config);

%% unpack variables
N       = config.N;
n_snr   = length(snr);
n_Q     = length(Q);

%% parfor setup
if flag_use_parfor == 1
	poolobj     = startparpool();
	parforArg   = poolobj.NumWorkers;
else
	parforArg = 0;
end

%% start generating
% Set the seed value and use the "twister" generator
rng(SEED, 'twister');

r_dfoVal = cell(n_Q, n_snr, Mc);
r_aoaVal = cell(n_Q, n_snr, Mc);
r_chlVal = cell(n_Q, n_snr, Mc);

fprintf(1, 'Generation Progress: %3d%%\n', 0);
for i_Mc = 1 : Mc
    % print out progress
    fprintf(1, '\b\b\b\b\b\b%5.2f%%', 100 * (i_Mc-1)/Mc);

    % check change channel
    if (i_Mc == 1) || (flag_change_ch_Mc > 0 && (mod(i_Mc - 1, flag_change_ch_Mc) == 0))
        if flag_simple_channel == 0
		    [dfo_cur, aoa_cur, chl_cur, preset_var] = ...
                gen_dfo_aoa_channel(Q, config);
        else
		    [dfo_cur, aoa_cur, chl_cur, preset_var] = ...
                gen_dfo_aoa_channel_simple(Q, config);
        end

    elseif regen_ch == 1
        if flag_simple_channel == 0
		    [dfo_cur, aoa_cur, chl_cur, preset_var] = ...
                gen_dfo_aoa_channel(Q, config, preset_var);
        else
		    [dfo_cur, aoa_cur, chl_cur, preset_var] = ...
                gen_dfo_aoa_channel_simple(Q, config, preset_var);
        end
    end

    for i_snr = 1 : n_snr
        r_dfoVal(:, i_snr, i_Mc) = dfo_cur;
        r_aoaVal(:, i_snr, i_Mc) = aoa_cur;
        r_chlVal(:, i_snr, i_Mc) = chl_cur;		
    end

end
fprintf(1, '\b\b\b\b\b\b%5.2f%%\n', 100);

generated_time = string(datetime('now'));

%% save files
destdirectory = fullfile('_gendata', sim_name);
if ~exist(destdirectory, 'dir')
	mkdir(destdirectory);
else
	delete(strcat(destdirectory, "\*"));
end

% 1. save config
fname = sprintf('%s_config.mat', sim_name);
fulldestination = fullfile(destdirectory, fname);
save(fulldestination, 'config');

fname = sprintf('%s_gen_config.mat', sim_name);
fulldestination = fullfile(destdirectory, fname);
save(fulldestination, 'Mc', 'snr', 'Q', "flag_simple_channel", ...
    "flag_change_ch_Mc", "generated_time");

% 2. save generate channel
fprintf(1, "Start saving ...\n");

PB = ProgressBar_cli_parfor(Mc);
tic;
parfor (i_Mc = 1 : Mc, parforArg)
	fname = sprintf('%s_%03d.mat', sim_name, i_Mc);
	fulldestination = fullfile(destdirectory, fname);

	dfoVal = r_dfoVal(:, :, i_Mc);
	aoaVal = r_aoaVal(:, :, i_Mc);
	chlVal = r_chlVal(:, :, i_Mc);

	parsave(fulldestination, dfoVal, aoaVal, chlVal);

	advance(PB);
end

fprintf(1, "fininshed\n");
