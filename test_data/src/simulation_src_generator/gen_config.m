%% configuration
config = struct;

% system model settings
config.N        = 1024;     % OFDM block size
config.Lcp      = 64;       % CP length
config.beta     = 0.5;      % roll-off ratio of raised-cosine pulse
% config.fc       = 60e+9;    % carrier frequency
config.fc       = 6000e+6;    % carrier frequency
config.fs       = 10e+6;    % sampling frequency
config.mod_type = 'QPSK';   % modulation type [fixed]
config.SP       = 0.5;      % spacing factor for ULA

% channel settings:
config.LOS      = 0;        % line of sight
config.v        = 300;      % velocity of TX (km/h)
config.DS       = 1e-6;     % delay spread
config.pg_low   = -20;      % lower bound of path gain
config.pg_high  = 0;        % upper bound of path gain
config.Lch      = 23;       % channel order
config.dfo_gap  = 0.046;    % smallest gap between two successive angles
config.aoa_gap  = 0.02;     % smallest gap between two successive angles

% simulation settings:
config.H            = 5;    % effetive central length of pulse
config.data_block   = 5;    % # of data_block to simulate BER


%% save config
destdirectory = '_config';
if ~exist(destdirectory, 'dir')
  mkdir(destdirectory);
end
fulldestination = fullfile(destdirectory, 'default_config.mat');
save(fulldestination, 'config');
disp("=====> config saved");
