
import MATLAB.Containers.Vector

flag_use_vec = 0;

addpath(genpath('src'))

sim_name = 'sim_03';
exp_name = 'exp_0322_0.9';

data_path = fullfile('_exp_dump_2', sim_name, exp_name);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mc = 2000;
Q = 1 : 4;
snr = 10 : 5 : 30;
N = 1024;

plot_snr = [0, 10, 20];
subplot_snr = [3, 1];
plot_Q = 1;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_Q = length(Q);
n_snr = length(snr);

ber_Mc_n = zeros(n_Q, n_snr, Mc);
ber_Mc_ml = zeros(n_Q, n_snr, Mc);

skip_collect = 0;

if skip_collect == 0
    for i_Mc = 1 : Mc
    
        disp(i_Mc)
        
        % %%%%%%%% load %%%%%%%%%%%%%%%%%%
        fname = sprintf('%s_r_avg_ber_%03d.ber.mat', sim_name, i_Mc);
        data = load(fullfile(data_path, fname));
        r_ber_Mc_n = data.r_ber_Mc_n;
        r_ber_Mc_ml = data.r_ber_Mc_ml;
    
        ber_Mc_n(:, :, i_Mc) = r_ber_Mc_n;
        ber_Mc_ml(:, :, i_Mc) = r_ber_Mc_ml;
    end
    
    avg_ber_n = mean(ber_Mc_n, 3);
    avg_ber_ml = mean(ber_Mc_ml, 3);
end

% =========================================================================
cmap = ['b', 'k', 'r', 'g', 'c', 'm', 'y'];
makr = ['o', '+', '*', '.', 'x', '-'];
lins = ["-", "--", ":", "-."];

% pos = get(fig, 'pos');
% fig.Position = pos;
for q = 1 : length(Q)
fig = figure(q);
	semilogy(snr, avg_ber_ml(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q)), ' (MLE)'], ...
		'Color', cmap(q), ...
		'Marker', makr(1), ...
		'LineStyle', lins(1));
	hold on;

	semilogy(snr, avg_ber_n(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q)), ' (Original)'], ...
		'Color', cmap(q), ...
		'Marker', makr(2), ...
		'LineStyle', lins(2));
	hold on;

hold off;
xlabel('snr');
ylabel('BER');
title_str = sprintf("Q = %d", (Q(q)));
title(title_str)
grid on;

legend('Location', 'southwest');
end