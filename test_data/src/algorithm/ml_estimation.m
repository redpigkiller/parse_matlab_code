function eps_hat = ml_estimation(Y, Q, flags, arg)
% ML_ESTIMATION Use ML estimation to estimate DFO (extend to ML estimation in paper1)
%   Y:  recv signal matrix
%   Q:  # of paths
%   algo:   choose the optimization algorithm (function name but not function pointer)
%       - Global optimization algorithm
%           + GA (Genetic Algorithm)
%           + PSO (Particle Swarm Optimization)
%           + SA (Simulated Annealing)
%           + PATTERNSEARCH (Pattern Search)
%       - Global optimization using local optimization
%           + MULTI_FMINUNC (Multistart with @fminunc)
%           + MULTI_FMINCON (Multistart with @fmincon)
%           + GLOBAL_FMINCON (GlobalSearch with @fmincon)
%       - Local optimization
%           + FMINUNC (@fminunc)
%           + FMINCON (@fmincon)
%           + FMINSEARCH (@fminsearch)
%           + MINFUNC
%   algoOpt:    options for optimization algorithm
%       use '' for default option
%       options states for b(bound), o(order)
%           + GA accept 'b', 'o'
%           + PSO accept 'b'
%           + SA accept 'b'
%           + PS accept 'b', 'o'
%           + MULTI_FMINUNC accept ''
%           + MULTI_FMINCON accept 'b', 'o'
%           + GLOBAL_FMINCON accept 'b', 'o'
%           + FMINUNC accept ''
%           + FMINCON accept 'b', 'o'
%           + FMINSEARCH accept ''
%           + MINFUNC accept ''
%   method: choose the methods to estimation DFO
%       1: estimate eps_q * M
%       2: directly estimate eps_q
%   arg:    needs to include information: fc, fs, v to calculate maximum
%           value of DFO
%
%
% NOTE1: for minFunc, please check out the following url:
%       https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
%
% NOTE2: reference paper: "Maximum Likelihood Timing and Carrier Frequency
%       Offset Estimation for OFDM Systems With Periodic Preambles"
%
% NOTE3: in this function, directly estimate autocorrelation matrix
%       of Y (set delta=0, lam=1)
%

%% check
if ~ismember(flags.method, [1 2])
	error("wrong method");
end

P = size(Y, 1);
M = size(Y, 2);

%% calculate estimated R
% R_hat = zeros(P, P);
% for m = 1 : M
% 	R_hat = R_hat + (Y(:, m) * Y(:, m)');
% end
% R_hat = 1/M * R_hat;
R_hat = 1/M * Y * (Y');

%% set up objective function
if isa(flags.method, 'function_handle')
	method = flags.method(Q);
else
    method = flags.method;
end

if method == 1
	cal_A = @(x) exp(1j * 2 * pi * (0:1:P-1).' * x);
else
	cal_A = @(x) exp(1j * 2 * pi * (0:1:P-1).' * x * M);
end

ML_func = @(y) -real(trace(R_hat * cal_A(y) * pinv(cal_A(y))));
% ML_func = @(y) -real(trace(R_hat * cal_A(y) * invvander(exp(2j*pi*y), P, true)));

%% constraints if needed
if isa(flags.algoOpt, 'function_handle')
	algoOpt = flags.algoOpt(Q);
elseif ~isa(flags.algoOpt, 'char')
	error("wrong algoOpt");
end

cstr_type = 0;
if contains(algoOpt, 'b')
	cstr_type = cstr_type + 2;
end

if contains(algoOpt, 'o')
	cstr_type = cstr_type + 1;
end

if bitand(cstr_type, 2)
	if method == 1
		lb = repmat(-.5, Q, 1);
		ub = repmat(.5, Q, 1);
	else
		lb = repmat(-.5/M, Q, 1);
		ub = repmat(.5/M, Q, 1);
	end
else
	lb = [];
	ub = [];
end

if bitand(cstr_type, 1) && Q > 1
	IJ = nchoosek(1:Q, 2);
	I = IJ(:, 1);
	J = IJ(:, 2);
	
	cstr_A = zeros(size(IJ, 1), Q);
	
	for q = 1 : size(IJ, 1)
		cstr_A(q, I(q)) = 1;
		cstr_A(q, J(q)) = -1;
	end
	cstr_b = zeros(size(cstr_A, 1), 1);
else
	cstr_A = [];
	cstr_b = [];
end

%% initial points
eps0 = linspace(-0.5, 0.5, Q);

%% setup
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(flags.algo, 'GA')
	options = optimoptions('ga', 'Display', 'off');
	f = ga(ML_func, Q, cstr_A, cstr_b, [], [], lb, ub, [], options);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'PSO')
% 	hybridopts = optimoptions('fmincon', 'Display', 'off', ...
%         'Algorithm', 'interior-point', ...
%         'OptimalityTolerance', 1e-15, ...
%         'StepTolerance', 1e-15);
%     options = optimoptions('particleswarm', 'Display', 'off',...
%         'SwarmSize', min(100, 10*Q), ...
%         'HybridFcn', { 'fmincon', hybridopts });
    options = optimoptions('particleswarm', 'Display', 'off', ...
        'SwarmSize', 250);
    f = particleswarm(ML_func, Q, lb, ub, options);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'SA')
	options = optimoptions('simulannealbnd', 'Display', 'off');
	f = simulannealbnd(ML_func, eps0, lb, ub, options);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'PATTERNSEARCH')
	options = optimoptions('patternsearch', 'Display', 'off');
	f = patternsearch(ML_func, eps0, cstr_A, cstr_b, [], [], lb, ub, [], options);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'MULTI_FMINUNC')
	options = optimoptions('fminunc', 'Display', 'off', ...
        'Algorithm', 'quasi-newton', ...
        'HessianApproximation', 'lbfgs', ...
        'OptimalityTolerance', 1e-15, ...
        'StepTolerance', 1e-15);
	problem = createOptimProblem('fminunc', 'objective', ML_func, ...
		'x0', eps0, 'options', options);
	ms = MultiStart('Display', 'off');
	f = run(ms, problem, 20);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'MULTI_FMINCON')
	options = optimoptions('fmincon', 'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'OptimalityTolerance', 1e-15, ...
        'StepTolerance', 1e-15);
	problem = createOptimProblem('fmincon', 'objective', ML_func, ...
		'x0', eps0, 'Aineq', cstr_A, 'bineq', cstr_b, ...
		'lb', lb, 'ub', ub, 'options', options);
	ms = MultiStart('Display', 'off');
	f = run(ms, problem, 20);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'GLOBAL_FMINCON')
	options = optimoptions('fmincon', 'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'OptimalityTolerance', 1e-15, ...
        'StepTolerance', 1e-15);
	problem = createOptimProblem('fmincon', 'objective', ML_func, ...
		'x0', eps0, 'Aineq', cstr_A, 'bineq', cstr_b, ...
		'lb', lb, 'ub', ub, 'options', options);
	gs = GlobalSearch('Display', 'off');
	f = run(gs, problem);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'FMINUNC')
	options = optimoptions('fminunc', 'Display', 'off', ...
        'Algorithm', 'quasi-newton', ...
        'HessianApproximation', 'lbfgs', ...
        'OptimalityTolerance', 1e-15, ...
        'StepTolerance', 1e-15);
	f = fminunc(ML_func, eps0, options);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'FMINCON')
	options = optimoptions('fmincon', 'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'OptimalityTolerance', 1e-15, ...
        'StepTolerance', 1e-15);
	f = fmincon(ML_func, eps0, cstr_A, cstr_b, [], [], lb, ub, [], options);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'FMINSEARCH')
	options = optimoptions('fminsearch', 'Display', 'off');
	f = fminsearch(ML_func, eps0, options);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(flags.algo, 'MINFUNC')
    options = struct;
    options.display = 'none';
    options.Method = 'lbfgs';
    options.numDiff = 2;
    % transpose the input of ML_func
    ML_func = @(y) ML_func(y.');
    f = minFunc(ML_func, eps0.', options).';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
	error("Wrong input algorithm")

end

%% output
% determine the multiples of 2*pi
if method == 2
	f = f * M;
end
k = round(f);

% obtain eps_hat
eps_hat = (f - k) ./ M;

% in case of eps_hat out of range
eps_q_max = arg.fc * (arg.v / 3.6) / physconst('LightSpeed') / arg.fs;

ind = abs(eps_hat) > eps_q_max;
eps_hat(ind) = sign(eps_hat(ind)) * eps_q_max;
eps_hat = eps_hat.';

end

