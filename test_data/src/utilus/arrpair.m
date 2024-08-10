function [arr_hat, indices, cost] = arrpair(arr_ref, arr_hat, cirbnd, method)
%ARRPAIR This function performs pairing between two arrays with same size, 
% it minimizes the "circular" distance and returns the paired second array.
% Make sure the input range is [-cirbnd, cirbnd], otherwise, the distance
% between two points may be negative. (cirbnd = inf by default)
% 
% Algorithm: Hungarian Algorithm / Munkres algorithm
% 
% Usage:
%   [arr_hat, indices, cost] = arrpair(arr_ref, arr_hat, cir_bnd, method);
% 
%   + input args:
%         eps:      input vector of existing angle
%         eps_hat:  input vector of measured angle
%         cir_bnd:  the angle is periodic. By default, cir_bnd = inf, i.e.
%                   the value of angle is in [-0.5, 0.5). For cir_bnd < 0,
%                   then there is no circular bound.
%         method:   0 for assignmentoptimal (default)
%                   1 for assignmentsuboptimal1
%                   2 for assignmentsuboptimal2
%     + output args:
%         eps_hat_paired:   output vector of measured angle paired with eps such that 2-norm of eps-eps_hat_paired is minimized
%         indices:          indices mapping from eps_hat to eps_hat_paired, i.e., eps_hat_paired = eps_hat(indices) in MATLAB
%         cost:             2-norm of eps-eps_hat_paired
%     
%   ref:
%         https://github.com/mcximing/hungarian-algorithm-cpp
%         https://github.com/phoemur/hungarian_algorithm
%         https://www.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem
%

error('can not find the mex compiled file.');

end