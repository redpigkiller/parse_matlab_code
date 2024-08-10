function parsave(savefile, varargin)
%PARSAVE The function allows user to save variables within parfor.
%   In parfor loop, MATLAB can NOT determine which variables the user wants
%   to save. Pass the variables into the PARSAVE function can resolve this
%   problem.
% 
%   Note that you should pass variables to parsave rather than the string
%   of the name of the variables. For example,
%   use 
%      PARSAVE(path, x);
%   instead of
%      PARSAVE(path, 'x');
%
%   See also PARFOR.

	name = cell(nargin - 1, 1);
	for ii = 1 : nargin - 1
		name{ii} = inputname(ii + 1);
		eval([name{ii} '=varargin{' num2str(ii) '};']);
		name{ii} = [',''',name{ii},''''];
	end
	varnames = [name{:}];
	comstring = ['save(''',savefile,'''',varnames,');'];
	eval(comstring);
end