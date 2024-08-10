sim_name = 'sim_release';
exp_name = 'exp_test_0222';

data_path = fullfile('_exp_dump', sim_name, exp_name);

Mc = 1;
Q = 1 : 4;
snr = -10 : 10 : 40;

n_Q = length(Q);
n_snr = length(snr);


r_ZF_err_pt = cell(n_Q, n_snr);
r_ZF_cor_pt = cell(n_Q, n_snr);
r_ZF_err2_pt = cell(n_Q, n_snr);
r_ZF_cor2_pt = cell(n_Q, n_snr);

r_MMSE_err_pt = cell(n_Q, n_snr);
r_MMSE_cor_pt = cell(n_Q, n_snr);
r_MMSE_err2_pt = cell(n_Q, n_snr);
r_MMSE_cor2_pt = cell(n_Q, n_snr);

e_ZF = zeros(n_Q, n_snr);
e_MMSE = zeros(n_Q, n_snr);
e_CORR = zeros(n_Q, n_snr);

for i_Mc = 1 : Mc
    i_Mc
    fname = sprintf('%s_r_avg_ber_%03d.ber.mat', sim_name, i_Mc);

    data = load(fullfile(data_path, fname));

    r_eZF_Mc = data.r_eZF_Mc;
    r_eMMSE_Mc = data.r_eMMSE_Mc;
    r_eCORR_Mc = data.r_eCORR_Mc;

    for i_Q = 1 : n_Q
        for i_snr = 1 : n_snr
            data_ZF = r_eZF_Mc{i_Q, i_snr};
            data_MMSE = r_eMMSE_Mc{i_Q, i_snr};
            data_CORR = r_eCORR_Mc{i_Q, i_snr};

            [n_trial_q, n_block] = size(data_ZF);
            c1 = n_trial_q * n_block;

            ZFval = 0;
            MMSEval = 0;
            CORRval = 0;

            ZF_err_pt = [];
            ZF_cor_pt = [];
            ZF_err2_pt = [];
            ZF_cor2_pt = [];

            MMSE_err_pt = [];
            MMSE_cor_pt = [];
            MMSE_err2_pt = [];
            MMSE_cor2_pt = [];

            for i_trial_q = 1 : n_trial_q
                for i_block = 1 : n_block
                    ZFval = ZFval + (1 / c1) * data_ZF{n_trial_q, n_block}.val;
                    MMSEval = MMSEval + (1 / c1) * data_MMSE{n_trial_q, n_block}.val;
                    CORRval = CORRval + (1 / c1) * data_CORR{n_trial_q, n_block};

ZF_err_pt = [ZF_err_pt; data_ZF{n_trial_q, n_block}.err];
ZF_cor_pt = [ZF_cor_pt; data_ZF{n_trial_q, n_block}.cor];
ZF_err2_pt = [ZF_err2_pt; data_ZF{n_trial_q, n_block}.err2];
ZF_cor2_pt = [ZF_cor2_pt; data_ZF{n_trial_q, n_block}.cor2];

MMSE_err_pt = [MMSE_err_pt; data_MMSE{n_trial_q, n_block}.err];
MMSE_cor_pt = [MMSE_cor_pt; data_MMSE{n_trial_q, n_block}.cor];
MMSE_err2_pt = [MMSE_err2_pt; data_MMSE{n_trial_q, n_block}.err2];
MMSE_cor2_pt = [MMSE_cor2_pt; data_MMSE{n_trial_q, n_block}.cor2];

                end

            end


            e_ZF(i_Q, i_snr) = e_ZF(i_Q, i_snr) + ZFval;
            e_MMSE(i_Q, i_snr) = e_MMSE(i_Q, i_snr) + MMSEval;
            e_CORR(i_Q, i_snr) = e_CORR(i_Q, i_snr) + CORRval;

r_ZF_err_pt{i_Q, i_snr} = [r_ZF_err_pt{i_Q, i_snr}; ZF_err_pt];
r_ZF_cor_pt{i_Q, i_snr} = [r_ZF_cor_pt{i_Q, i_snr}; ZF_cor_pt];
r_ZF_err2_pt{i_Q, i_snr} = [r_ZF_err2_pt{i_Q, i_snr}; ZF_err2_pt];
r_ZF_cor2_pt{i_Q, i_snr} = [r_ZF_cor2_pt{i_Q, i_snr}; ZF_cor2_pt];
r_MMSE_err_pt{i_Q, i_snr} = [r_MMSE_err_pt{i_Q, i_snr}; MMSE_err_pt];
r_MMSE_cor_pt{i_Q, i_snr} = [r_MMSE_cor_pt{i_Q, i_snr}; MMSE_cor_pt];
r_MMSE_err2_pt{i_Q, i_snr} = [r_MMSE_err2_pt{i_Q, i_snr}; MMSE_err2_pt];
r_MMSE_cor2_pt{i_Q, i_snr} = [r_MMSE_cor2_pt{i_Q, i_snr}; MMSE_cor2_pt];

        end
    end

end

e_ZF = e_ZF / Mc;
e_MMSE = e_MMSE / Mc;
e_CORR = e_CORR / Mc;

%% plot

cmap = ['b', 'k', 'r', 'g', 'c', 'm', 'y'];
makr = ['o', '+', '*', '.', 'x', '-'];
lins = ["-", "--", ":", "-."];

snr_plot(1, snr, e_ZF, Q, 'snr', '', '');

fig = figure(2);
for q = 1 : length(Q)
	semilogy(snr, e_ZF(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(3), ...
		'LineStyle', lins(1));
	hold on;
end
hold off;
xlabel('snr');
% ylabel('MSE DFO');
grid on;

% title_str = sprintf('N=%d, M=%d, Mc=%d', config.N, config.M, sim_config.Mc);
% title(title_str);

legend('Location', 'southwest');

return
%% 
fig = figure(2);
for q = 1 : length(Q)
	semilogy(snr, e_MMSE(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(3), ...
		'LineStyle', lins(1));
	hold on;
end
hold off;
xlabel('snr');
% ylabel('MSE DFO');
grid on;

% title_str = sprintf('N=%d, M=%d, Mc=%d', config.N, config.M, sim_config.Mc);
% title(title_str);

legend('Location', 'southwest');

%%
fig = figure(3);
for q = 1 : length(Q)
	semilogy(snr, e_CORR(q, :), ...
		'DisplayName', ['Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(3), ...
		'LineStyle', lins(1));
	hold on;
end
hold off;
xlabel('snr');
% ylabel('MSE DFO');
grid on;

% title_str = sprintf('N=%d, M=%d, Mc=%d', config.N, config.M, sim_config.Mc);
% title(title_str);

legend('Location', 'southwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ZF_err = zeros(n_Q, n_snr);
ZF_cor = zeros(n_Q, n_snr);
MMSE_err = zeros(n_Q, n_snr);
MMSE_cor = zeros(n_Q, n_snr);
for i_snr = 1 : n_snr
    for i_Q = 1 : n_Q
        ZF_err(i_Q, i_snr) = nan2zero(mean(r_ZF_err_pt{i_Q, i_snr}));
        ZF_cor(i_Q, i_snr) = nan2zero(mean(r_ZF_cor_pt{i_Q, i_snr}));
        MMSE_err(i_Q, i_snr) = nan2zero(mean(r_MMSE_err_pt{i_Q, i_snr}));
        MMSE_cor(i_Q, i_snr) = nan2zero(mean(r_MMSE_cor_pt{i_Q, i_snr}));
    end
end

fig = figure(4);
for q = 1 : length(Q)
	semilogy(snr, ZF_err(q, :), ...
		'DisplayName', ['error Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(1), ...
		'LineStyle', lins(1));
    semilogy(snr, ZF_cor(q, :), ...
		'DisplayName', ['correct Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(2), ...
		'LineStyle', lins(2));
	hold on;
end
hold off;
xlabel('snr');
grid on;
title('ZF');
legend('Location', 'southwest');

fig = figure(5);
for q = 1 : length(Q)
	semilogy(snr, MMSE_err(q, :), ...
		'DisplayName', ['error Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(1), ...
		'LineStyle', lins(1));
    semilogy(snr, MMSE_cor(q, :), ...
		'DisplayName', ['correct Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(2), ...
		'LineStyle', lins(2));
	hold on;
end
hold off;
xlabel('snr');
grid on;
title('MMSE');
legend('Location', 'southwest');





%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ZF_err2 = zeros(n_Q, n_snr);
ZF_cor2 = zeros(n_Q, n_snr);
MMSE_err2 = zeros(n_Q, n_snr);
MMSE_cor2 = zeros(n_Q, n_snr);
for i_snr = 1 : n_snr
    for i_Q = 1 : n_Q
        ZF_err2(i_Q, i_snr) = nan2zero(mean(r_ZF_err2_pt{i_Q, i_snr}));
        ZF_cor2(i_Q, i_snr) = nan2zero(mean(r_ZF_cor2_pt{i_Q, i_snr}));
        MMSE_err2(i_Q, i_snr) = nan2zero(mean(r_MMSE_err2_pt{i_Q, i_snr}));
        MMSE_cor2(i_Q, i_snr) = nan2zero(mean(r_MMSE_cor2_pt{i_Q, i_snr}));
    end
end

fig = figure(6);
for q = 1 : length(Q)
	semilogy(snr, ZF_err2(q, :), ...
		'DisplayName', ['error Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(1), ...
		'LineStyle', lins(1));
    semilogy(snr, ZF_cor2(q, :), ...
		'DisplayName', ['correct Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(2), ...
		'LineStyle', lins(2));
	hold on;
end
hold off;
xlabel('snr');
grid on;
title('ZF');
legend('Location', 'southwest');

fig = figure(7);
for q = 1 : length(Q)
	semilogy(snr, MMSE_err2(q, :), ...
		'DisplayName', ['error Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(1), ...
		'LineStyle', lins(1));
    semilogy(snr, MMSE_cor2(q, :), ...
		'DisplayName', ['correct Q=', num2str(Q(q))], ...
		'Color', cmap(q), ...
		'Marker', makr(2), ...
		'LineStyle', lins(2));
	hold on;
end
hold off;
xlabel('snr');
grid on;
title('MMSE');
legend('Location', 'southwest');


return;
%%
x_ZF_err_pt = {};
ZF_err_pt = [];
for i_snr = 1 : n_snr
    for i_Q = 1 : n_Q
        nlen = length(r_ZF_cor_pt{i_Q, i_snr});
        if nlen == 0
            continue;
        end
        Qlabel = repmat(num2str( (i_snr-1)*n_Q + i_Q ), nlen, 1);
        x_ZF_err_pt = vertcat(x_ZF_err_pt , cellstr(Qlabel));
        ZF_err_pt = [ZF_err_pt ; r_ZF_cor_pt{i_Q, i_snr}]; 
    end
end

x_ZF_cor_pt = {};
ZF_cor_pt = [];
for i_snr = 1 : n_snr
    for i_Q = 1 : n_Q
        nlen = length(r_ZF_cor_pt{i_Q, i_snr});
        if nlen == 0
            continue;
        end
        Qlabel = repmat(num2str( (i_snr-1)*n_Q + i_Q ), nlen, 1);
        x_ZF_cor_pt = vertcat(x_ZF_cor_pt , cellstr(Qlabel));
        ZF_cor_pt = [ZF_cor_pt ; r_ZF_cor_pt{i_Q, i_snr}]; 
    end
end

x1 = categorical(x_ZF_err_pt);
x2 = categorical(x_ZF_cor_pt);

figure(8)
scatter(x1, ZF_err_pt, 10, 'red', 'filled');
hold on;
scatter(x2, ZF_cor_pt, 10, 'blue', 'filled');
grid on;
hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_MMSE_err_pt = {};
MMSE_err_pt = [];
for i_snr = 1 : n_snr
    for i_Q = 1 : n_Q
        nlen = length(r_MMSE_err_pt{i_Q, i_snr});
        if nlen == 0
            continue;
        end
        Qlabel = repmat(num2str( (i_snr-1)*n_Q + i_Q ), nlen, 1);
        x_MMSE_err_pt = vertcat(x_MMSE_err_pt , cellstr(Qlabel));
        MMSE_err_pt = [MMSE_err_pt ; r_MMSE_err_pt{i_Q, i_snr}]; 
    end
end
%%
x_MMSE_cor_pt = {};
MMSE_cor_pt = [];
for i_snr = 1 : n_snr
    for i_Q = 1 : n_Q
        nlen = length(r_MMSE_cor_pt{i_Q, i_snr});
        if nlen == 0
            continue;
        end
        Qlabel = repmat(num2str( (i_snr-1)*n_Q + i_Q ), nlen, 1);
        x_MMSE_cor_pt = vertcat(x_MMSE_cor_pt , cellstr(Qlabel));
        MMSE_cor_pt = [MMSE_cor_pt ; r_MMSE_cor_pt{i_Q, i_snr}]; 
    end
end

x1 = categorical(x_MMSE_err_pt);
x2 = categorical(x_MMSE_cor_pt);

figure(9)
scatter(x1, MMSE_err_pt, 10, 'red', 'filled');
hold on;
scatter(x2, MMSE_cor_pt, 10, 'blue', 'filled');
grid on;
hold off;
