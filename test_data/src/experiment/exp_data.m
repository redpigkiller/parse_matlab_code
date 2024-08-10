
import MATLAB.Containers.Vector

flag_use_vec = 0;

addpath(genpath('src'))

sim_name = 'sim_release';
exp_name = 'exp_0224';

data_path = fullfile('_exp_dump', sim_name, exp_name);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mc = 500;
Q = 1 : 4;
snr = -10:10:40;
N = 1024;

plot_snr = [0, 10, 20];
subplot_snr = [3, 1];
plot_Q = 1;

thres_partition = 40;
ratio_in_mean = 1/10;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_Q = length(Q);
n_snr = length(snr);

skip_collect = 1;

if skip_collect == 0

% 0: record val
e0_ZF_est = zeros(n_Q, n_snr, Mc);
e0_ZF_idl = zeros(n_Q, n_snr, Mc);

e0_MMSE_est = zeros(n_Q, n_snr, Mc);
e0_MMSE_idl = zeros(n_Q, n_snr, Mc);

% 1, 2: record pts
e1_ZF_err_pts = cell(n_Q, n_snr, Mc);
e1_ZF_cor_pts = cell(n_Q, n_snr, Mc);
e1_MMSE_err_pts = cell(n_Q, n_snr, Mc);
e1_MMSE_cor_pts = cell(n_Q, n_snr, Mc);

e2_ZF_err_pts = cell(n_Q, n_snr, Mc);
e2_ZF_cor_pts = cell(n_Q, n_snr, Mc);
e2_MMSE_err_pts = cell(n_Q, n_snr, Mc);
e2_MMSE_cor_pts = cell(n_Q, n_snr, Mc);

tic;
% parfor i_Mc = 1 : Mc
for i_Mc = 1 : Mc

    disp(i_Mc)
    
    % %%%%%%%% load %%%%%%%%%%%%%%%%%%
    fname = sprintf('%s_r_avg_ber_%03d.ber.mat', sim_name, i_Mc);
    data = load(fullfile(data_path, fname));
    r_eZF_Mc = data.r_eZF_Mc;
    r_eMMSE_Mc = data.r_eMMSE_Mc;
%     r_eCORR_Mc = data.r_eCORR_Mc;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i_Q = 1 : n_Q
        for i_snr = 1 : n_snr
            data_ZF = r_eZF_Mc{i_Q, i_snr};
            data_MMSE = r_eMMSE_Mc{i_Q, i_snr};
%             data_CORR = r_eCORR_Mc{i_Q, i_snr};

            [n_trial_q, n_block] = size(data_ZF);

            % 0
            r_e0_ZF_est = 0;
            r_e0_ZF_idl = 0;
            r_e0_MMSE_est = 0;
            r_e0_MMSE_idl = 0;

            % 1, 2
            if flag_use_vec == 0
                r_e1_ZF_err_pts = [];
                r_e1_ZF_cor_pts = [];
                r_e1_MMSE_err_pts = [];
                r_e1_MMSE_cor_pts = [];

                r_e2_ZF_err_pts = [];
                r_e2_ZF_cor_pts = [];
                r_e2_MMSE_err_pts = [];
                r_e2_MMSE_cor_pts = [];

            else
                pr_e1_ZF_err_pts = Vector();
                pr_e1_ZF_cor_pts = Vector();
                pr_e1_MMSE_err_pts = Vector();
                pr_e1_MMSE_cor_pts = Vector();
    
    
                pr_e2_ZF_err_pts = Vector();
                pr_e2_ZF_cor_pts = Vector();
                pr_e2_MMSE_err_pts = Vector();
                pr_e2_MMSE_cor_pts = Vector();

            end


            for i_trial_q = 1 : n_trial_q
                for i_block = 1 : n_block
                    %%% 0               
                    r_e0_ZF_est = r_e0_ZF_est + ...
                        data_ZF{i_trial_q, i_block}.e0_est;
                    r_e0_ZF_idl = r_e0_ZF_idl + ...
                        data_ZF{i_trial_q, i_block}.e0_idl;
                    r_e0_MMSE_est = r_e0_MMSE_est + ...
                        data_MMSE{i_trial_q, i_block}.e0_est;
                    r_e0_MMSE_idl = r_e0_MMSE_idl + ...
                        data_MMSE{i_trial_q, i_block}.e0_idl;

                    if flag_use_vec == 0
%%% 1
r_e1_ZF_err_pts = [r_e1_ZF_err_pts; data_ZF{i_trial_q, i_block}.e1_err_pts];
r_e1_ZF_cor_pts = [r_e1_ZF_cor_pts; data_ZF{i_trial_q, i_block}.e1_cor_pts];
r_e1_MMSE_err_pts = [r_e1_MMSE_err_pts; data_MMSE{i_trial_q, i_block}.e1_err_pts];
r_e1_MMSE_cor_pts = [r_e1_MMSE_cor_pts; data_MMSE{i_trial_q, i_block}.e1_cor_pts];

%%% 2
r_e2_ZF_err_pts = [r_e2_ZF_err_pts; data_ZF{i_trial_q, i_block}.e2_err_pts];
r_e2_ZF_cor_pts = [r_e2_ZF_cor_pts; data_ZF{i_trial_q, i_block}.e2_cor_pts];
r_e2_MMSE_err_pts = [r_e2_MMSE_err_pts; data_MMSE{i_trial_q, i_block}.e2_err_pts];
r_e2_MMSE_cor_pts = [r_e2_MMSE_cor_pts; data_MMSE{i_trial_q, i_block}.e2_cor_pts];

                    else
%%% 1
pr_e1_ZF_err_pts.PushBack(data_ZF{i_trial_q, i_block}.e1_err_pts);
pr_e1_ZF_cor_pts.PushBack(data_ZF{i_trial_q, i_block}.e1_cor_pts);
pr_e1_MMSE_err_pts.PushBack(data_MMSE{i_trial_q, i_block}.e1_err_pts);
pr_e1_MMSE_cor_pts.PushBack(data_MMSE{i_trial_q, i_block}.e1_cor_pts);

%%% 2
pr_e2_ZF_err_pts.PushBack(data_ZF{i_trial_q, i_block}.e2_err_pts);
pr_e2_ZF_cor_pts.PushBack(data_ZF{i_trial_q, i_block}.e2_cor_pts);
pr_e2_MMSE_err_pts.PushBack(data_MMSE{i_trial_q, i_block}.e2_err_pts);
pr_e2_MMSE_cor_pts.PushBack(data_MMSE{i_trial_q, i_block}.e2_cor_pts);

                    end

                end

            end

%%% 0
e0_ZF_est(i_Q, i_snr, i_Mc) = r_e0_ZF_est / n_trial_q / n_block;
e0_ZF_idl(i_Q, i_snr, i_Mc) = r_e0_ZF_idl / n_trial_q / n_block;
e0_MMSE_est(i_Q, i_snr, i_Mc) = r_e0_MMSE_est / n_trial_q / n_block;
e0_MMSE_idl(i_Q, i_snr, i_Mc) = r_e0_MMSE_idl / n_trial_q / n_block;

if flag_use_vec == 0
%%% 1
e1_ZF_err_pts{i_Q, i_snr, i_Mc} = r_e1_ZF_err_pts;
e1_ZF_cor_pts{i_Q, i_snr, i_Mc} = r_e1_ZF_cor_pts;
e1_MMSE_err_pts{i_Q, i_snr, i_Mc} = r_e1_MMSE_err_pts;
e1_MMSE_cor_pts{i_Q, i_snr, i_Mc} = r_e1_MMSE_cor_pts;
%%% 2
e2_ZF_err_pts{i_Q, i_snr, i_Mc} = r_e2_ZF_err_pts;
e2_ZF_cor_pts{i_Q, i_snr, i_Mc} = r_e2_ZF_cor_pts;
e2_MMSE_err_pts{i_Q, i_snr, i_Mc} = r_e2_MMSE_err_pts;
e2_MMSE_cor_pts{i_Q, i_snr, i_Mc} = r_e2_MMSE_cor_pts;
else
%%% 1
e1_ZF_err_pts{i_Q, i_snr, i_Mc} = pr_e1_ZF_err_pts.Data;
e1_ZF_cor_pts{i_Q, i_snr, i_Mc} = pr_e1_ZF_cor_pts.Data;
e1_MMSE_err_pts{i_Q, i_snr, i_Mc} = pr_e1_MMSE_err_pts.Data;
e1_MMSE_cor_pts{i_Q, i_snr, i_Mc} = pr_e1_MMSE_cor_pts.Data;
%%% 2
e2_ZF_err_pts{i_Q, i_snr, i_Mc} = pr_e2_ZF_err_pts.Data;
e2_ZF_cor_pts{i_Q, i_snr, i_Mc} = pr_e2_ZF_cor_pts.Data;
e2_MMSE_err_pts{i_Q, i_snr, i_Mc} = pr_e2_MMSE_err_pts.Data;
e2_MMSE_cor_pts{i_Q, i_snr, i_Mc} = pr_e2_MMSE_cor_pts.Data;
end
        end
    end
end

fprintf('collecting time: %g\n', toc);

%%       post-processing
% 0
e0_ZF_est = sum(e0_ZF_est, 3) ./ Mc;
e0_ZF_idl = sum(e0_ZF_idl, 3) ./ Mc;
e0_MMSE_est = sum(e0_MMSE_est, 3) ./ Mc;
e0_MMSE_idl = sum(e0_MMSE_idl, 3) ./ Mc;

% 1
mean_ZF_err1 = zeros(n_Q, n_snr);
mean_ZF_cor1 = zeros(n_Q, n_snr);
mean_MMSE_err1 = zeros(n_Q, n_snr);
mean_MMSE_cor1 = zeros(n_Q, n_snr);

e1_ZF_err_pts = cellcat(e1_ZF_err_pts);
e1_ZF_cor_pts = cellcat(e1_ZF_cor_pts);
e1_MMSE_err_pts = cellcat(e1_MMSE_err_pts);
e1_MMSE_cor_pts = cellcat(e1_MMSE_cor_pts);
for i_snr = 1 : n_snr
    for i_Q = 1 : n_Q
        mean_ZF_err1(i_Q, i_snr) = ...
            nan2zero(mean(e1_ZF_err_pts{i_Q, i_snr}));
        mean_ZF_cor1(i_Q, i_snr) = ...
            nan2zero(mean(e1_ZF_cor_pts{i_Q, i_snr}));
        mean_MMSE_err1(i_Q, i_snr) = ...
            nan2zero(mean(e1_MMSE_err_pts{i_Q, i_snr}));
        mean_MMSE_cor1(i_Q, i_snr) = ...
            nan2zero(mean(e1_MMSE_cor_pts{i_Q, i_snr}));
    end
end

% 2
mean_ZF_err2 = zeros(n_Q, n_snr);
mean_ZF_cor2 = zeros(n_Q, n_snr);
mean_MMSE_err2 = zeros(n_Q, n_snr);
mean_MMSE_cor2 = zeros(n_Q, n_snr);

e2_ZF_err_pts = cellcat(e2_ZF_err_pts);
e2_ZF_cor_pts = cellcat(e2_ZF_cor_pts);
e2_MMSE_err_pts = cellcat(e2_MMSE_err_pts);
e2_MMSE_cor_pts = cellcat(e2_MMSE_cor_pts);

for i_snr = 1 : n_snr
    for i_Q = 1 : n_Q
        mean_ZF_err2(i_Q, i_snr) = ...
            nan2zero(mean(e2_ZF_err_pts{i_Q, i_snr}));
        mean_ZF_cor2(i_Q, i_snr) = ...
            nan2zero(mean(e2_ZF_cor_pts{i_Q, i_snr}));
        mean_MMSE_err2(i_Q, i_snr) = ...
            nan2zero(mean(e2_MMSE_err_pts{i_Q, i_snr}));
        mean_MMSE_cor2(i_Q, i_snr) = ...
            nan2zero(mean(e2_MMSE_cor_pts{i_Q, i_snr}));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmap = ['b', 'k', 'r', 'g', 'c', 'm', 'y'];
makr = ['o', '+', '*', '.', 'x', '-'];
lins = ["-", "--", ":", "-."];


close all

%% 1
% =========================================================================
snr_plot(1, snr, e0_ZF_est, Q, 'snr', '', 'ZF (estimation) e0');

%% 2
% =========================================================================
snr_plot(2, snr, e0_ZF_idl, Q, 'snr', '', 'ZF (ideal) e0');

%% 3
% =========================================================================
snr_plot(3, snr, e0_MMSE_est, Q, 'snr', '', 'MMSE (estimation) e0');

%% 4
% =========================================================================
snr_plot(4, snr, e0_MMSE_idl, Q, 'snr', '', 'MMSE (ideal) e0');

%% 5
% % =========================================================================
% fig = figure(5);
% pos = get(fig, 'pos');
% fig.Position = pos;
% for q = 1 : length(Q)
% 	semilogy(snr, mean_ZF_err1(q, :), ...
% 		'DisplayName', ['Q=', num2str(Q(q)), ' (error)'], ...
% 		'Color', cmap(q), ...
% 		'Marker', makr(1), ...
% 		'LineStyle', lins(1));
% 	hold on;
% end
% for q = 1 : length(Q)
% 	semilogy(snr, mean_ZF_cor1(q, :), ...
% 		'DisplayName', ['Q=', num2str(Q(q)), ' (correct)'], ...
% 		'Color', cmap(q), ...
% 		'Marker', makr(2), ...
% 		'LineStyle', lins(2));
% 	hold on;
% end
% hold off;
% xlabel('snr');
% ylabel('value');
% title("mean of ZF (e1)")
% grid on;
% 
% legend('Location', 'southwest');

%% 6
% =========================================================================
% fig = figure(6);
% pos = get(fig, 'pos');
% fig.Position = pos;
% for q = 1 : length(Q)
% 	semilogy(snr, mean_MMSE_err1(q, :), ...
% 		'DisplayName', ['Q=', num2str(Q(q)), ' (error)'], ...
% 		'Color', cmap(q), ...
% 		'Marker', makr(1), ...
% 		'LineStyle', lins(1));
% 	hold on;
% end
% for q = 1 : length(Q)
% 	semilogy(snr, mean_MMSE_cor1(q, :), ...
% 		'DisplayName', ['Q=', num2str(Q(q)), ' (correct)'], ...
% 		'Color', cmap(q), ...
% 		'Marker', makr(2), ...
% 		'LineStyle', lins(2));
% 	hold on;
% end
% hold off;
% xlabel('snr');
% ylabel('value');
% title("mean of MMSE (e1)")
% grid on;
% 
% legend('Location', 'southwest');

%% 7
% =========================================================================
fig = figure(7);
pos = get(fig, 'pos');
fig.Position = pos;
for q = 1 : length(Q)
	semilogy(snr, mean_ZF_err2(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q)), ' (error)'], ...
		'Color', cmap(q), ...
		'Marker', makr(1), ...
		'LineStyle', lins(1));
	hold on;
end
for q = 1 : length(Q)
	semilogy(snr, mean_ZF_cor2(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q)), ' (correct)'], ...
		'Color', cmap(q), ...
		'Marker', makr(2), ...
		'LineStyle', lins(2));
	hold on;
end
hold off;
xlabel('snr');
ylabel('value');
title("mean of ZF (e2)")
grid on;

legend('Location', 'southwest');

%% 8
% =========================================================================
fig = figure(8);
pos = get(fig, 'pos');
fig.Position = pos;
for q = 1 : length(Q)
	semilogy(snr, mean_MMSE_err2(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q)), ' (error)'], ...
		'Color', cmap(q), ...
		'Marker', makr(1), ...
		'LineStyle', lins(1));
	hold on;
end
for q = 1 : length(Q)
	semilogy(snr, mean_MMSE_cor2(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q)), ' (correct)'], ...
		'Color', cmap(q), ...
		'Marker', makr(2), ...
		'LineStyle', lins(2));
	hold on;
end
hold off;
xlabel('snr');
ylabel('value');
title("mean of MMSE (e2)")
grid on;

legend('Location', 'southwest');

%% 13
% % =========================================================================
% % ZF e1
% % =========================================================================
% fig = figure(9);
% pos = get(fig, 'pos');
% fig.Position = [pos(1) pos(2)-420 560 840];
% 
% for idx = 1 : length(plot_snr)
%     subplot(subplot_snr(1), subplot_snr(2), idx);
% 
%     i_snr = find(snr == plot_snr(idx));
% 
%     count_TP = zeros(n_Q, thres_partition);
%     count_FN = zeros(n_Q, thres_partition);
%     count_FP = zeros(n_Q, thres_partition);
%     count_TN = zeros(n_Q, thres_partition);
% 
%     for i_Q = 1 : n_Q
%         if plot_Q > 0 && plot_Q ~= Q(i_Q)
%             continue;
%         end
% 
%         a = mean_ZF_cor1(i_Q, i_snr);
%         b = mean_ZF_err1(i_Q, i_snr);
%         if a > b
%             [b, a] = deal(a, b);
%         end
% 
%         thres_gap = (b - a) / floor(thres_partition*ratio_in_mean);
%         thres_lo_count = floor((thres_partition - floor(thres_partition*ratio_in_mean)) / 2);
%         thres_up_count = thres_partition - floor(thres_partition*ratio_in_mean) - thres_lo_count;
%         
%         lb = a - thres_gap * thres_lo_count;
%         if lb < 0;    lb = 0;    end
%         ub = b + thres_gap * thres_up_count;
%         
%         thres = linspace(lb, ub, thres_partition);
% 
%         for i_thres = 1 : length(thres)
%             tot = length(e1_ZF_cor_pts{i_Q, i_snr}) + ...
%                     length(e1_ZF_err_pts{i_Q, i_snr});
%             if mod(tot, N) ~= 0
%                 tot
%                 error("11");
%             end
% 
%             count_TP(i_Q, i_thres) = ...
%                 nnz(e1_ZF_err_pts{i_Q, i_snr} >= thres(i_thres));
% 
%             count_FN(i_Q, i_thres) = ...
%                 nnz(e1_ZF_err_pts{i_Q, i_snr} < thres(i_thres));
% 
%             count_FP(i_Q, i_thres) = ...
%                 nnz(e1_ZF_cor_pts{i_Q, i_snr} >= thres(i_thres));
% 
%             count_TN(i_Q, i_thres) = ...
%                 nnz(e1_ZF_cor_pts{i_Q, i_snr} < thres(i_thres));
%         end
%     end
% 
%     for i_Q = 1 : n_Q
%         if plot_Q > 0 && plot_Q ~= Q(i_Q)
%             continue;
%         end
% 
%         FNR = count_FN(i_Q, :) ./ (count_TP(i_Q, :) + count_FN(i_Q, :));
%         FPR = count_FP(i_Q, :) ./ (count_FP(i_Q, :) + count_TN(i_Q, :));
%         
%         semilogy(thres, FNR, ...
% 		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (FNR)'], ...
% 		    'Color', cmap(i_Q), ...
% 		    'Marker', makr(1), ...
% 		    'LineStyle', lins(1));
%         hold on;
% 
%         semilogy(thres, FPR, ...
% 		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (FPR)'], ...
% 		    'Color', cmap(i_Q), ...
% 		    'Marker', makr(3), ...
% 		    'LineStyle', lins(2));
%         hold on;
%     end
% 	hold off;
%     grid on;
% 
%     xlabel('threshold');
%     ylabel('probabilty');
% 
% 	title_str = sprintf('snr = %d ZF (e1)', snr(i_snr));
% 	title(title_str);
%     title_str = sprintf('Q = %d ZF (e1)', Q(i_Q));
% 	sgtitle(title_str);
% 	
% 	legend('Location', 'southwest');
% end


%% 14
% % =========================================================================
% % MMSE e1
% % =========================================================================
% fig = figure(10);
% pos = get(fig, 'pos');
% fig.Position = [pos(1) pos(2)-420 560 840];
% 
% for idx = 1 : length(plot_snr)
%     subplot(subplot_snr(1), subplot_snr(2), idx);
% 
%     i_snr = find(snr == plot_snr(idx));
% 
%     count_TP = zeros(n_Q, thres_partition);
%     count_FN = zeros(n_Q, thres_partition);
%     count_FP = zeros(n_Q, thres_partition);
%     count_TN = zeros(n_Q, thres_partition);
% 
%     for i_Q = 1 : n_Q
%         if plot_Q > 0 && plot_Q ~= Q(i_Q)
%             continue;
%         end
% 
%         a = mean_MMSE_cor1(i_Q, i_snr);
%         b = mean_MMSE_err1(i_Q, i_snr);
%         if a > b
%             [b, a] = deal(a, b);
%         end
% 
%         thres_gap = (b - a) / floor(thres_partition*ratio_in_mean);
%         thres_lo_count = floor((thres_partition - floor(thres_partition*ratio_in_mean)) / 2);
%         thres_up_count = thres_partition - floor(thres_partition*ratio_in_mean) - thres_lo_count;
%         
%         lb = a - thres_gap * thres_lo_count;
%         if lb < 0;    lb = 0;    end
%         ub = b + thres_gap * thres_up_count;
%         
%         thres = linspace(lb, ub, thres_partition);
% 
%         for i_thres = 1 : length(thres)
%             tot = length(e1_MMSE_cor_pts{i_Q, i_snr}) + ...
%                     length(e1_MMSE_err_pts{i_Q, i_snr});
%             if mod(tot, N) ~= 0
%                 tot
%                 error("11");
%             end
% 
%             count_TP(i_Q, i_thres) = ...
%                 nnz(e1_MMSE_err_pts{i_Q, i_snr} >= thres(i_thres));
% 
%             count_FN(i_Q, i_thres) = ...
%                 nnz(e1_MMSE_err_pts{i_Q, i_snr} < thres(i_thres));
% 
%             count_FP(i_Q, i_thres) = ...
%                 nnz(e1_MMSE_cor_pts{i_Q, i_snr} >= thres(i_thres));
% 
%             count_TN(i_Q, i_thres) = ...
%                 nnz(e1_MMSE_cor_pts{i_Q, i_snr} < thres(i_thres));
%         end
%     end
% 
%     for i_Q = 1 : n_Q
%         if plot_Q > 0 && plot_Q ~= Q(i_Q)
%             continue;
%         end
% 
%         FNR = count_FN(i_Q, :) ./ (count_TP(i_Q, :) + count_FN(i_Q, :));
%         FPR = count_FP(i_Q, :) ./ (count_FP(i_Q, :) + count_TN(i_Q, :));
%         
%         semilogy(thres, FNR, ...
% 		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (FNR)'], ...
% 		    'Color', cmap(i_Q), ...
% 		    'Marker', makr(1), ...
% 		    'LineStyle', lins(1));
%         hold on;
% 
%         semilogy(thres, FPR, ...
% 		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (FPR)'], ...
% 		    'Color', cmap(i_Q), ...
% 		    'Marker', makr(3), ...
% 		    'LineStyle', lins(2));
%         hold on;
%     end
% 	hold off;
%     grid on;
% 
%     xlabel('threshold');
%     ylabel('probabilty');
% 
% 	title_str = sprintf('snr = %d MMSE (e1)', snr(i_snr));
% 	title(title_str);
%     title_str = sprintf('Q = %d MMSE (e1)', Q(i_Q));
% 	sgtitle(title_str);
% 	
% 	legend('Location', 'southwest');
% end

%% 15
% =========================================================================
% ZF e2
% =========================================================================
% fig = figure(11);
% pos = get(fig, 'pos');
% fig.Position = [pos(1) pos(2)-420 560 840];
% 
% for idx = 1 : length(plot_snr)
%     subplot(subplot_snr(1), subplot_snr(2), idx);
% 
%     i_snr = find(snr == plot_snr(idx));
% 
%     count_TP = zeros(n_Q, thres_partition);
%     count_FN = zeros(n_Q, thres_partition);
%     count_FP = zeros(n_Q, thres_partition);
%     count_TN = zeros(n_Q, thres_partition);
% 
%     for i_Q = 1 : n_Q
%         if plot_Q > 0 && plot_Q ~= Q(i_Q)
%             continue;
%         end
% 
%         a = mean_ZF_cor2(i_Q, i_snr);
%         b = mean_ZF_err2(i_Q, i_snr);
%         if a > b
%             [b, a] = deal(a, b);
%         end
% 
%         thres_gap = (b - a) / floor(thres_partition*ratio_in_mean);
%         thres_lo_count = floor((thres_partition - floor(thres_partition*ratio_in_mean)) / 2);
%         thres_up_count = thres_partition - floor(thres_partition*ratio_in_mean) - thres_lo_count;
%         
%         lb = a - thres_gap * thres_lo_count;
%         if lb < 0;    lb = 0;    end
%         ub = b + thres_gap * thres_up_count;
%         
%         thres = linspace(lb, ub, thres_partition);
% 
%         for i_thres = 1 : length(thres)
%             tot = length(e2_ZF_cor_pts{i_Q, i_snr}) + ...
%                     length(e2_ZF_err_pts{i_Q, i_snr});
%             if mod(tot, N) ~= 0
%                 tot
%                 error("11");
%             end
% 
%             count_TP(i_Q, i_thres) = ...
%                 nnz(e2_ZF_err_pts{i_Q, i_snr} >= thres(i_thres));
% 
%             count_FN(i_Q, i_thres) = ...
%                 nnz(e2_ZF_err_pts{i_Q, i_snr} < thres(i_thres));
% 
%             count_FP(i_Q, i_thres) = ...
%                 nnz(e2_ZF_cor_pts{i_Q, i_snr} >= thres(i_thres));
% 
%             count_TN(i_Q, i_thres) = ...
%                 nnz(e2_ZF_cor_pts{i_Q, i_snr} < thres(i_thres));
%         end
%     end
% 
%     for i_Q = 1 : n_Q
%         if plot_Q > 0 && plot_Q ~= Q(i_Q)
%             continue;
%         end
% 
%         FNR = count_FN(i_Q, :) ./ (count_TP(i_Q, :) + count_FN(i_Q, :));
%         FPR = count_FP(i_Q, :) ./ (count_FP(i_Q, :) + count_TN(i_Q, :));
%         
%         semilogy(thres, FNR, ...
% 		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (FNR)'], ...
% 		    'Color', cmap(i_Q), ...
% 		    'Marker', makr(1), ...
% 		    'LineStyle', lins(1));
%         hold on;
% 
%         semilogy(thres, FPR, ...
% 		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (FPR)'], ...
% 		    'Color', cmap(i_Q), ...
% 		    'Marker', makr(3), ...
% 		    'LineStyle', lins(2));
%         hold on;
%     end
% 	hold off;
%     grid on;
% 
%     xlabel('threshold');
%     ylabel('probabilty');
% 
% 	title_str = sprintf('snr = %d ZF (e2)', snr(i_snr));
% 	title(title_str);
%     title_str = sprintf('Q = %d ZF (e2)', Q(i_Q));
% 	sgtitle(title_str);
% 	
% 	legend('Location', 'southwest');
% end

%% 16
% =========================================================================
% MMSE e2
% =========================================================================
fig = figure(12);
pos = get(fig, 'pos');
fig.Position = [pos(1) pos(2)-420 560 840];

for idx = 1 : length(plot_snr)
    subplot(subplot_snr(1), subplot_snr(2), idx);

    i_snr = find(snr == plot_snr(idx));

    count_TP = zeros(n_Q, thres_partition);
    count_FN = zeros(n_Q, thres_partition);
    count_FP = zeros(n_Q, thres_partition);
    count_TN = zeros(n_Q, thres_partition);

    for i_Q = 1 : n_Q
        if plot_Q > 0 && plot_Q ~= Q(i_Q)
            continue;
        end

        a = mean_MMSE_cor2(i_Q, i_snr);
        b = mean_MMSE_err2(i_Q, i_snr);
        if a > b
            [b, a] = deal(a, b);
        end

        thres_gap = (b - a) / floor(thres_partition*ratio_in_mean);
        thres_lo_count = floor((thres_partition - floor(thres_partition*ratio_in_mean)) / 2);
        thres_up_count = thres_partition - floor(thres_partition*ratio_in_mean) - thres_lo_count;
        
        lb = a - thres_gap * thres_lo_count;
        if lb < 0;    lb = 0;    end
        ub = b + thres_gap * thres_up_count;
        lb = 0;
        ub = 1;
        thres = linspace(lb, ub, thres_partition);

        for i_thres = 1 : length(thres)
            tot = length(e2_MMSE_cor_pts{i_Q, i_snr}) + ...
                    length(e2_MMSE_err_pts{i_Q, i_snr});
            if mod(tot, N) ~= 0
                tot
                error("11");
            end

            count_TP(i_Q, i_thres) = ...
                nnz(e2_MMSE_err_pts{i_Q, i_snr} >= thres(i_thres));

            count_FN(i_Q, i_thres) = ...
                nnz(e2_MMSE_err_pts{i_Q, i_snr} < thres(i_thres));

            count_FP(i_Q, i_thres) = ...
                nnz(e2_MMSE_cor_pts{i_Q, i_snr} >= thres(i_thres));

            count_TN(i_Q, i_thres) = ...
                nnz(e2_MMSE_cor_pts{i_Q, i_snr} < thres(i_thres));
        end
        
    end

    for i_Q = 1 : n_Q
        if plot_Q > 0 && plot_Q ~= Q(i_Q)
            continue;
        end

        FNR = count_FN(i_Q, :) ./ (count_TP(i_Q, :) + count_FN(i_Q, :));
        FPR = count_FP(i_Q, :) ./ (count_FP(i_Q, :) + count_TN(i_Q, :));
        ACC = (count_TP(i_Q, :) + count_TN(i_Q, :)) ./ ...
            (count_TP(i_Q, :) + count_FN(i_Q, :) + count_FP(i_Q, :) + count_TN(i_Q, :));
        Precision = count_TP(i_Q, :) ./ (count_TP(i_Q, :) + count_FP(i_Q, :));
        Recall = count_TP(i_Q, :) ./ (count_TP(i_Q, :) + count_FN(i_Q, :));
        F1_score = 2*(Precision.*Recall)./(Precision+Recall);


        semilogy(thres, F1_score, ...
		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (F1_score)'], ...
		    'Color', cmap(i_Q), ...
		    'Marker', makr(1), ...
		    'LineStyle', lins(1));
        hold on;

%         semilogy(thres, FNR, ...
% 		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (FNR)'], ...
% 		    'Color', cmap(i_Q), ...
% 		    'Marker', makr(1), ...
% 		    'LineStyle', lins(1));
%         hold on;

%         semilogy(thres, FPR, ...
% 		    'DisplayName', ['Q=', num2str(Q(i_Q)), ' (FPR)'], ...
% 		    'Color', cmap(i_Q), ...
% 		    'Marker', makr(3), ...
% 		    'LineStyle', lins(2));
%         hold on;
    end
	hold off;
    grid on;

    xlabel('threshold');
    ylabel('probabilty');

	title_str = sprintf('snr = %d MMSE (e2)', snr(i_snr));
	title(title_str);
    title_str = sprintf('Q = %d MMSE (e2)', plot_Q);
	sgtitle(title_str);
	
	legend('Location', 'southwest');
end