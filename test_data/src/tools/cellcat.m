function Y = cellcat(X)
%CELLCAT Concatenate vertically along the 3rd dimension of the cell array
%with each entry constains a 1D array
%   Consider X = cell(4, 6, 100), and X{i, j, k} = [ 1.23, 4.56, ... ] is a
%   1D array. The output Y is a cell array with size (4, 6), and Y{i, j} =
%   [ X{i, j, 1} ; X{i, j, 2} ; ... ; X{i, j, 100} ].
%   Note that the 1D array in X should be a row vector.
%
% Reference: https://www.mathworks.com/matlabcentral/answers/447469-looping-through-a-matrix-of-unknown-dimensions
%

nv   = ndims(X);
v    = num2cell(ones(1, nv));
vLim = size(X);

Y = cell(vLim(1:end-1));

ready = false;
while ~ready
    % Do what you need with X and the index vector v
    Y{v{1:end-1}} = [Y{v{1:end-1}} ; X{v{:}}];
    
    % Update the index vector:
    ready = true;           % Assume that the WHILE loop is ready
    for k = nv : -1 : 1
        v{k} = v{k} + 1;
        if v{k} <= vLim(k)
            ready = false;  % No, WHILE loop is not ready now
            break;          % v(k) increased successfully, leave "for k" loop
        end
        v{k} = 1;           % Reset v(k), proceed to next k
    end
end

end

