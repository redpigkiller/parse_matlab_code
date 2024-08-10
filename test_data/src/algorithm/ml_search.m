function [s_hat, n_err] = ml_search(yp, ss_hat, chl_mtx, arg)
%ML_SEARCH Use ML to correct some erro in hard decision
%   yp:         recv signal
%   ss_hat:     equalizer output
%   chl_mtx:    channel matrix
%   arg:        options     
%
% Try to improve using DP...
%

% hard decision
s_hat = nrSymbolDemodulate(ss_hat, 'QPSK', 'DecisionType', 'Hard');

% recover symbol
ss_hat_d = nrSymbolModulate(s_hat, 'QPSK');

% compute distance:
eVal = abs(ss_hat - ss_hat_d).^2;

% eVal(924)
% find indices of symbol whose eVal >= threshold
eVal(eVal < arg.thres) = -1;

% keep largest max_error_symbol_number symbols
[~, e_idx] = maxk(eVal, arg.max_error_symbol_number);
% e_idx.'
% eVal(e_idx)
n_err = length(e_idx);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try all possible ss_hat_d:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. calculate correct part of H_hat:
h_cor = chl_mtx(:, setdiff(1 : size(chl_mtx, 2), e_idx)) * ...
    ss_hat_d(setdiff(1 : size(chl_mtx, 2), e_idx));

% 2. iterate through all possible ss_hat_d
yp = yp - h_cor;
try_value = nrSymbolModulate([0 0 0 1 1 0 1 1].', 'QPSK');
% try_value = pskmod(0:3, 4, pi/4).';

%% brute-force
min_val = inf;

for idx = 1 : length(try_value)^n_err
    % Convert i to a base-N number to get the indices of the elements to use from the values array
    indices = baseN(idx - 1, length(try_value), n_err) + 1;
    % Use the indices to get the corresponding elements from the values array
    try_s = try_value(indices);
    % Calculate the value of f for this array
    val = norm(yp - chl_mtx(:, e_idx) * try_s)^2;
    % Update the minimum value and corresponding array if necessary
    if val < min_val
        min_val = val;
        s_hat = try_s;
    end
end

% get s_hat
ss_hat_d(e_idx) = s_hat;

% demodulation
s_hat = nrSymbolDemodulate(ss_hat_d, 'QPSK', 'DecisionType', 'Hard');

end