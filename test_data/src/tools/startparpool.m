function poolobj = startparpool(n_workers, idletimeout)
%STARTPARPOOL Start the parpool with specified number of workers.
%   STARTPARPOOL() returns pool object with maximum number of workers.
%
%   STARTPARPOOL(n_workers) returns pool object with n_workers workers.
%
%   STARTPARPOOL(n_workers, idletimeout) returns pool object with n_workers
%   workers and with idle timeout being idletimeout.
%
%   See also PARPOOL.

	narginchk(0, 2);

	if nargin < 1
		n_workers = 0;
	end
	if nargin < 2
		idletimeout = 0;
	end

	% get the pool object
	poolobj = gcp("nocreate");

	% check if there is actually a pool
	if ~isempty(poolobj)
		if n_workers > 0
			if poolobj.NumWorkers ~= n_workers
				delete(poolobj);
			end
		end
	end

	% finally, if there is not a pool
	if isempty(gcp('nocreate'))
		if n_workers > 0
			if idletimeout > 0
				poolobj = parpool(n_workers, 'IdleTimeout', idletimeout);
			else
				poolobj = parpool(n_workers);
			end
		else
			if idletimeout > 0
				poolobj = parpool('IdleTimeout', idletimeout);
			else
				poolobj = parpool();
			end
		end
	end
end